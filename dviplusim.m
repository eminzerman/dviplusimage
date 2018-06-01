function dviplusim(fileFolderName, outFolder)
%DVIPLUSIM Creates DVI Plus image for SIM2 HDR47 screen.
%   DVIPLUSIM(FILEFOLDERNAME) Creates DVI Plus image and writes three  
%   images as output: 
%       * a PNG image, which is used to show the image in DVI+ mode
%       * an EXR image with estimated emitted luminance information
%       * an EXR image with estimated backlight luminance information
%   The function requires FILEFOLDERNAME as input. This input can both 
%   be a file or directory.
%
% Examples:
%   dviplusim('hdrImages\cloudyDay.hdr');
%   dviplusim('hdrImages\cloudyDay.hdr', 'outFolder');
%   dviplusim([hdrPath filesep 'hdrImagesFolder'], 'outFolder');
%
%   Check the following paper for details:
%     E. Zerman, G. Valenzise, and F. Dufaux, "A Dual Modulation Algorithm 
%     for Accurate Reproduction of High Dynamic Range Video", IEEE 12th 
%     Image, Video, and Multidimensional Signal Processing Workshop (IVMSP), 
%     Bordeaux, France, July 2016.
%
% ---------------------
% - Emin Zerman / emin.zerman@telecom-paristech.fr
% - Created:  23/03/2015 -- Updated: 16/09/2015
% - Telecom ParisTech - TSI - MM
% ---------------------
%  Copyright (C) 2018 - Emin Zerman
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------

% Check if flags exist and initialize
if(~exist('outFolder', 'var')), 
    outFolder ='dvipOut'; 
    
    % Create outFolder either in target folder or code directory
    if isdir(fileFolderName)
        if ~isdir([fileFolderName filesep outFolder])
            mkdir([fileFolderName filesep outFolder]);
        end
        outFolder = [fileFolderName filesep outFolder];
    else
        if ~isdir(outFolder)
            mkdir(outFolder);
        end
    end
else
    if ~isdir(outFolder)
        mkdir(outFolder);
    end
end
Params.outFolder = outFolder;
Params.fileFolderName = fileFolderName;

% CONSTANT values
Params.LED_MAX_LUM = 180;
Params.LUMINANCE_HDR_CONSTANT = 179;
Params.MAX_LUMINANCE = 4250;
Params.EXTRAPOLATE_MARGIN = 80;
Params.LCD_LEAKAGE_CONSTANT = 1/0.005;
Params.fileFlag = 0;

% Get point spread function, resize PSF to reduce complexity
load('led_psf.mat');
Params.led_psf = led_psf;
clear led_psf;
%Params.led_psf = lum(hdrimread('psf_model.pfm'));
Params.led_psfSm = imresize(Params.led_psf, 1/4); 

% Find locations of LEDs, resize LED maps to reduce complexity
[~, Params.ledLabels, Params.ledIndices] = findledpositions([1080 1920], 0); 
[Params.ledMapSm, Params.ledLabelsSm, Params.ledIndicesSm] = findledpositions([1080 1920]/4, 0);

% Check if the first input is a file name or folder name
[Params, DbContents] = determineFileFolder(Params);


%% Process the indicated image or all images in the indicated folder
for imNo = Params.fileIndices
    
    % Find fileName for naming
    if Params.fileFlag
        [~, Ims.fileName, Ims.fileXt] = fileparts(Params.fileFolderName);
    else
        [~, Ims.fileName, Ims.fileXt] = fileparts(DbContents(imNo).name);
    end
    
    %% Acquire Image
    % Acquire HDR image by utilizing Matlab's hdrread function, and tonemap
    % function for user to see the image.
    [Ims, rFlag] = readHdrIm(DbContents(imNo).name, Params, Ims);
    if rFlag == 1
        continue;
    end

    %% Find desired backlight
    Ims = desiredBacklight(Params, Ims);

    %% Find LED values
    Ims = findLedValues(Params, Ims);

    %% Finding LCD values
    % After LED values and backlight illumination caused by LEDs have been
    % found, calculate and set LCD values for DVI Plus image.
    Ims = findLcdValues(Params, Ims);

    %% Save DVI+ output image
        
    % Save DVI+ Image
    getdviplusim(round(255*Ims.ledVals)',uint8(Ims.lcdGamCorrdUint8),...
                           [Params.outFolder filesep Ims.fileName '.png']);

    % Find the reconstructed HDR image
    [dvipHdr, ledLum3d] = dviplus2hdrim(round(255*Ims.ledVals)',...
                            uint8(Ims.lcdGamCorrdUint8), Params.led_psf);
    exrwrite(dvipHdr, [Params.outFolder filesep 'lum_' Ims.fileName '.exr' ]);
    exrwrite(ledLum3d, [Params.outFolder filesep 'ledLum_' Ims.fileName '.exr' ]);
    
end

end

function [Params, DbContents] = determineFileFolder(Params)

    if ~isempty(strfind(Params.fileFolderName(end-4:end),'.hdr')) || ...
            ~isempty(strfind(Params.fileFolderName(end-4:end),'.exr'))
        % ========== It is an HDR or EXR file ==========
        Params.fileFlag = 1;
        Params.fileIndices = 1;
        DbContents.name = Params.fileFolderName;
    elseif isdir(Params.fileFolderName)
        % ========== It is a folder ==========
        Params.fileFlag = 0;
        % Find folder contents
        countHdrFiles = length(dir([Params.fileFolderName '\*.hdr']));
        countExrFiles = length(dir([Params.fileFolderName '\*.exr']));
        if countHdrFiles == 0 && countExrFiles == 0
            error('No .hdr or .ext files present in specified path!');
        else
            if countHdrFiles > countExrFiles
                % Read HDR
                DbContents = dir([Params.fileFolderName '\*.hdr']);
            else
                %Read EXR
                DbContents = dir([Params.fileFolderName '\*.exr']);
            end
        end
        Params.fileIndices = 1:length(DbContents);
    else
        error('This is not an HDR image file or a folder!');
    end

end

function [Ims, rFlag] = readHdrIm(imName, Params, Ims)
    rFlag = 0;
    if Params.fileFlag
        if strcmp(Ims.fileXt,'.hdr')
            try
                hdrIm = hdrread(imName); 
            catch
                try
                    hdrIm = hdrimread(imName);
                catch
                    error('Problem acquiring the HDR image/frame!!');
                end
            end
        elseif strcmp(Ims.fileXt,'.exr')
            hdrIm = exrread(imName);
        elseif strcmp(Ims.fileXt,'.pfm')
            hdrIm = hdrimread(imName);
        else
            rFlag = 1;
        end
    else
        if strcmp(Ims.fileXt,'.hdr')
            try
                hdrIm = hdrread([Params.fileFolderName filesep imName]); 
            catch
                try
                    hdrIm = hdrimread([Params.fileFolderName filesep imName]);
                catch
                    error('Problem acquiring the HDR image/frame!!');
                end
            end
        elseif strcmp(Ims.fileXt,'.exr')
            hdrIm = exrread([Params.fileFolderName filesep imName]);
        elseif strcmp(Ims.fileXt,'.pfm')
            hdrIm = hdrimread([Params.fileFolderName filesep imName]);
        else
            rFlag = 1;
        end
    end
    
    % Adjust image to be 1920/1080 by resizing and/or padding with black
    [Ims.hdrIm, Ims.lastIm] = makeitfullhd(hdrIm);
    
    % Find target luminance, multiply w/ constant if ".hdr"
    if strcmp(Ims.fileXt,'.hdr')
        hdrImL = Ims.hdrIm*Params.LUMINANCE_HDR_CONSTANT;
        %lastImL = Ims.lastIm*Params.LUMINANCE_HDR_CONSTANT;
    elseif strcmp(Ims.fileXt,'.exr')
        hdrImL = Ims.hdrIm;
        %lastImL = Ims.lastIm;
    end
    hdrLum = max(hdrImL, [], 3);

    if max(max(hdrLum))>Params.MAX_LUMINANCE, 
        hdrLum(hdrLum>Params.MAX_LUMINANCE) = Params.MAX_LUMINANCE; 
    end

    Ims.hdrLum = hdrLum;
end

function [hdrIm, lastIm, typeFlag, diff] = makeitfullhd(hdrIm)

ratioHD = 1920/1080;
    ratioIm = size(hdrIm,2)/size(hdrIm,1);
    
    if (ratioHD == ratioIm) && (size(hdrIm,2) ~= 1920)
        % Same resolution with different aspect ratio
        hdrIm = imresize(hdrIm, [1080 1920]); 
        hdrIm(hdrIm<0) = 0; % Avoid having negative values
        lastIm = hdrIm;
        typeFlag = 1;       % Type 1: Same resolution w/ diff. aspect ratio
        diff = 0;
        
    elseif (ratioHD == ratioIm) && (size(hdrIm,2) == 1920)
        % Same resolution
        lastIm = hdrIm;
        typeFlag = 2;       % Type 2: Same resolution
        diff = 0;
        
    elseif (ratioHD > ratioIm)
        % Different ratio, resize the image
        tempIm = imresize(hdrIm, [1080 ratioIm*1080]); 
        tempIm(tempIm<0) = 0;   % Avoid having negative values
        
        % Find the difference between 1920 and the column pixels
        diff = 1920 - size(tempIm,2);
        lastIm = zeros(1080,1920,size(tempIm,3));
        lastIm(:,floor(diff/2)+1:(floor(diff/2)+size(tempIm,2)),:) = tempIm;
        hdrIm = lastIm;     % lastIm is the black padded image
        
        hdrIm(:,1:floor(diff/2),:) = repmat(hdrIm(:,floor(diff/2)+1,:), [1 floor(diff/2) 1]);
        hdrIm(:,(floor(diff/2)+size(tempIm,2))+1:end,:) = repmat(hdrIm(:,(floor(diff/2)+size(tempIm,2)),:), [1 diff-floor(diff/2) 1]);
        typeFlag = 3;       % Type 3: Diff ratio, tight frame
        
    elseif (ratioHD < ratioIm)
        % Different ratio, resize the image
        tempIm = imresize(hdrIm, [1920/ratioIm 1920]); 
        tempIm(tempIm<0) = 0;   % Avoid having negative values
        
        % Find the difference between 1080 and the row pixels
        diff = 1080 - size(tempIm,1);
        lastIm = zeros(1080,1920,size(tempIm,3));
        lastIm(floor(diff/2)+1:(floor(diff/2)+size(tempIm,1)),:,:) = tempIm;
        hdrIm = lastIm;     % lastIm is the black padded image
        
        hdrIm(1:floor(diff/2),:,:) = repmat(hdrIm(floor(diff/2)+1,:,:), [floor(diff/2) 1 1]);
        hdrIm((floor(diff/2)+size(tempIm,1))+1:end,:,:) = repmat(hdrIm((floor(diff/2)+size(tempIm,1)),:,:), [diff-floor(diff/2) 1 1]);
        typeFlag = 4;       % Type 4: Diff ratio, wide frame
        
    else
        error('Something wrong with the spatial resolution!!');
    end
end

function [ledImage, ledLabels, ledIndices, ledPos] = findledpositions(imSize, debugFlag)

    if size(imSize)<2,               error('Wrong imsize!');   end
    if(~exist('debugFlag', 'var')),  debugFlag = 0;            end

    %% Determine the locations of LEDs
    % Since LEDs of SIM2 HDR47S4 display is put in a hexagonal order, the odd
    % and even rows have different distance values from the edge of the screen.
    % The calculations below find the positions of LEDs from origin, and they
    % can be used for plotting.

    dispXOdd = 0.64:17.6:(17.6*59+0.64*2);
    dispXEven = (0.64+17.6/2):17.6:(17.6*58+0.64*2+17.6);
    dispYOdd = -8.01:-31.6:-(31.6*18+8.01*2);
    dispYEven = -(8.01+31.6/2):-31.6:-(31.6*17+8.01*2+31.6);

    % scale if necessary (i.e. if resolution is different - N/A for SIM2)
    scale = 1080/imSize(1);
    dispXOdd = dispXOdd/scale;
    dispXEven = dispXEven/scale;
    dispYOdd = dispYOdd/scale;
    dispYEven = dispYEven/scale;

    % Create a grid for odd and even rows
    [ledsOddY, ledsOddX] = meshgrid(dispYOdd, dispXOdd);
    ledsOddX = reshape(ledsOddX, [], 1);
    ledsOddY = reshape(ledsOddY, [], 1);
    [ledsEvenY, ledsEvenX] = meshgrid(dispYEven, dispXEven);
    ledsEvenX = reshape(ledsEvenX, [], 1);
    ledsEvenY = reshape(ledsEvenY, [], 1);

    % Find LED positions by alternately concatanating odd and even rows
    ledPosX = [];
    ledPosY = [];
    for ind = 1:37
        if mod(ind,2)==1
            ledPosX = cat(1, ledPosX, ledsOddX((ceil(ind/2)-1)*60+1:ceil(ind/2)*60));
            ledPosY = cat(1, ledPosY, ledsOddY((ceil(ind/2)-1)*60+1:ceil(ind/2)*60));
        else
            ledPosX = cat(1, ledPosX, ledsEvenX((ceil(ind/2)-1)*59+1:ceil(ind/2)*59));
            ledPosY = cat(1, ledPosY, ledsEvenY((ceil(ind/2)-1)*59+1:ceil(ind/2)*59));
        end
    end
    ledPos = [ledPosX ledPosY (1:2202)'];

    if debugFlag, figure, plot(ledPosX, ledPosY, 'r*'); axis equal, axis([0 1039.68 -584.82 0]); end;
    if debugFlag, figure, voronoi(ledPosX, ledPosY); axis equal, axis([0 1039.68 -584.82 0]); end;

    %% LEDs Pixel Positions
    % In order to handle other image related operations, LED pixel locations
    % have to be known. Hence, pixel locations are estimated considering the
    % true locations of LEDs which is found above.

    tempLedImage = zeros(imSize(1:2));
    tempLedLabels = zeros(imSize(1:2));

    % Find row and col pixels
    pixRows = ceil(ledPos(:,2)/-0.5415);
    pixCols = ceil(ledPos(:,1)/0.5415);
    ledIndices = sub2ind(imSize(1:2), pixRows, pixCols);

    % Assign values
    tempLedImage(ledIndices) = 1;
    tempLedLabels(ledIndices) = ledPos(:,3);

    if debugFlag, figure, imshow(tempLedImage); end;
    if debugFlag, figure, imagesc(tempLedLabels); end;

    %% Output
    ledImage = tempLedImage;
    ledLabels = tempLedLabels;
end

function Ims = desiredBacklight(Params, Ims)

    % Extrapolate (by inpainting) the image to avoid information loss 
    % while filtering. Downsample the hdrLum by 4 to reduce complexity.
    hdrLumBar = zeros( (Params.EXTRAPOLATE_MARGIN+1080)/4,...
                       (Params.EXTRAPOLATE_MARGIN+1920)/4 ); 
    hdrLumBar(:) = NaN;
    hdrLumBar(Params.EXTRAPOLATE_MARGIN/8 +1 : 1080/4 +Params.EXTRAPOLATE_MARGIN/8,...
              Params.EXTRAPOLATE_MARGIN/8 +1 : 1920/4 +Params.EXTRAPOLATE_MARGIN/8) ...
              = imresize(Ims.hdrLum, 1/4);
    hdrLumBarN = inpaintn(hdrLumBar);
    hdrLumExtp = imresize(hdrLumBarN, 4);
    hdrLumExtp(Params.EXTRAPOLATE_MARGIN/2 +1 : 1080 +Params.EXTRAPOLATE_MARGIN/2,...
                  Params.EXTRAPOLATE_MARGIN/2 +1 : 1920 +Params.EXTRAPOLATE_MARGIN/2) ...
                  = Ims.hdrLum;
    
    % -- Find Target Luminance value --
    % Find the minimum required luminance for fidelity by filtering the max
    % value image with a Gaussian filter (i.e. fspecial('disk',N) )
    maxValImg = imdilate(hdrLumExtp*1.1, strel('disk', 30));
    targetLumMin = conv2(maxValImg, fspecial('disk',30),'same');
    
    % Find the maximum allowed luminance to avoid any LCD leakage and light
    % halos that may occur. Median filter in order to reduce the artifacts
    % created due to the high frequency characteristics (e.g. stars, pixel
    % defects caused during image acquisition, etc.)
    targetLumMax = hdrLumExtp*Params.LCD_LEAKAGE_CONSTANT; 
    targetLumMax(targetLumMax>Params.MAX_LUMINANCE) = Params.MAX_LUMINANCE; 
    targetLumMax = medfilt2(targetLumMax, [5 5]);
    
    % Find the minumum value (or the constraining value) of these minumum 
    % and maximum images.
    tempArr = cat(3, targetLumMin, targetLumMax);
    targetLum_pre = min(tempArr, [], 3); 
    
    % Smooth out the preliminary result to avoid direct faults of 
    % individual LEDs. Median filter in order to reduce the artifacts
    % created due to the high frequency characteristics (e.g. stars, pixel
    % defects caused during image acquisition, etc.)
    targetLum_filt = imfilter(targetLum_pre, fspecial('disk', Params.EXTRAPOLATE_MARGIN));
    targetLum_temp = targetLum_filt(Params.EXTRAPOLATE_MARGIN/2 +1 : 1080 +Params.EXTRAPOLATE_MARGIN/2,...
                                    Params.EXTRAPOLATE_MARGIN/2 +1 : 1920 +Params.EXTRAPOLATE_MARGIN/2);
    targetLum = medfilt2(max(Ims.hdrLum,targetLum_temp), [5 5]); 
    targetLum(targetLum==0) = targetLum_temp(targetLum==0);
    targetLumSm = imresize(targetLum, 1/4); 
    targetLumSm(targetLumSm < 0) = 0;

    % Pass the pictures
    Ims.targetLum = targetLum;
    Ims.targetLumSm = targetLumSm;
end    

function y = inpaintn(x,n,y0,m)

    if nargin==0&&nargout==0, RunTheExample, return, end

    x = double(x);
    if nargin==1 || isempty(n), n = 100; end

    sizx = size(x);
    d = ndims(x);
    Lambda = zeros(sizx);
    for i = 1:d
        siz0 = ones(1,d);
        siz0(i) = sizx(i);
        Lambda = bsxfun(@plus,Lambda,...
            cos(pi*(reshape(1:sizx(i),siz0)-1)/sizx(i)));
    end
    Lambda = 2*(d-Lambda);

    % Initial condition
    W = isfinite(x);
    if nargin==3 && ~isempty(y0)
        y = y0;
        s0 = 3; % note: s = 10^s0
    else
        if any(~W(:))
            [y,s0] = InitialGuess(x,isfinite(x));
        else
            y = x;
            return
        end
    end
    x(~W) = 0;

    if isempty(n) || n<=0, n = 100; end

    % Smoothness parameters: from high to negligible values
    s = logspace(s0,-6,n); 

    RF = 2; % relaxation factor

    if nargin<4 || isempty(m), m = 2; end
    Lambda = Lambda.^m;

    % h = waitbar(0,'Inpainting...');
    for i = 1:n
            Gamma = 1./(1+s(i)*Lambda);
            y = RF*idctn(Gamma.*dctn(W.*(x-y)+y)) + (1-RF)*y;
            % waitbar(i/n,h)
    end
    % close(h)

    y(W) = x(W);

end

%% Initial Guess
function [z,s0] = InitialGuess(y,I)

    if license('test','image_toolbox')
        %-- nearest neighbor interpolation
        [~,L] = bwdist(I);
        z = y;
        z(~I) = y(L(~I));
        s0 = 3; % note: s = 10^s0
    else
        warning('MATLAB:inpaintn:InitialGuess',...
            ['BWDIST (Image Processing Toolbox) does not exist. ',...
            'The initial guess may not be optimal; additional',...
            ' iterations can thus be required to ensure complete',...
            ' convergence. Increase N value if necessary.'])
        z = y;
        z(~I) = mean(y(I));
        s0 = 6; % note: s = 10^s0
    end

end

%% Example (3-D)
function RunTheExample
          load wind
          xmin = min(x(:)); xmax = max(x(:)); %#ok
          zmin = min(z(:)); ymax = max(y(:)); %#ok
          %-- wind velocity
          vel0 = interp3(sqrt(u.^2+v.^2+w.^2),1,'cubic');
          x = interp3(x,1); y = interp3(y,1); z = interp3(z,1);
          %-- remove randomly 90% of the data
          I = randperm(numel(vel0));
          velNaN = vel0;
          velNaN(I(1:round(numel(I)*.9))) = NaN;
          %-- inpaint using INPAINTN
          vel = inpaintn(velNaN);
          %-- display the results
          subplot(221), imagesc(velNaN(:,:,15)), axis equal off
          title('Corrupt plane, z = 15')
          subplot(222), imagesc(vel(:,:,15)), axis equal off
          title('Reconstructed plane, z = 15')    
          subplot(223)
          hsurfaces = slice(x,y,z,vel0,[xmin,100,xmax],ymax,zmin);
          set(hsurfaces,'FaceColor','interp','EdgeColor','none')
          hcont = contourslice(x,y,z,vel0,[xmin,100,xmax],ymax,zmin);
          set(hcont,'EdgeColor',[.7,.7,.7],'LineWidth',.5)
          view(3), daspect([2,2,1]), axis tight
          title('Actual data compared with...')
          subplot(224)
          hsurfaces = slice(x,y,z,vel,[xmin,100,xmax],ymax,zmin);
          set(hsurfaces,'FaceColor','interp','EdgeColor','none')
          hcont = contourslice(x,y,z,vel,[xmin,100,xmax],ymax,zmin);
          set(hcont,'EdgeColor',[.7,.7,.7],'LineWidth',.5)
          view(3), daspect([2,2,1]), axis tight
          title('... reconstructed data')
end

%% DCTN
function y = dctn(y)


    y = double(y);
    sizy = size(y);
    y = squeeze(y);
    dimy = ndims(y);

    % Some modifications are required if Y is a vector
    if isvector(y)
        dimy = 1;
        if size(y,1)==1, y = y.'; end
    end

    % Weighting vectors
    w = cell(1,dimy);
    for dim = 1:dimy
        n = (dimy==1)*numel(y) + (dimy>1)*sizy(dim);
        w{dim} = exp(1i*(0:n-1)'*pi/2/n);
    end

    % --- DCT algorithm ---
    if ~isreal(y)
        y = complex(dctn(real(y)),dctn(imag(y)));
    else
        for dim = 1:dimy
            siz = size(y);
            n = siz(1);
            y = y([1:2:n 2*floor(n/2):-2:2],:);
            y = reshape(y,n,[]);
            y = y*sqrt(2*n);
            y = ifft(y,[],1);
            y = bsxfun(@times,y,w{dim});
            y = real(y);
            y(1,:) = y(1,:)/sqrt(2);
            y = reshape(y,siz);
            y = shiftdim(y,1);
        end
    end
            
    y = reshape(y,sizy);

end

%% IDCTN
function y = idctn(y)

    y = double(y);
    sizy = size(y);
    y = squeeze(y);
    dimy = ndims(y);

    % Some modifications are required if Y is a vector
    if isvector(y)
        dimy = 1;
        if size(y,1)==1
            y = y.';
        end
    end

    % Weighing vectors
    w = cell(1,dimy);
    for dim = 1:dimy
        n = (dimy==1)*numel(y) + (dimy>1)*sizy(dim);
        w{dim} = exp(1i*(0:n-1)'*pi/2/n);
    end

    % --- IDCT algorithm ---
    if ~isreal(y)
        y = complex(idctn(real(y)),idctn(imag(y)));
    else
        for dim = 1:dimy
            siz = size(y);
            n = siz(1);
            y = reshape(y,n,[]);
            y = bsxfun(@times,y,w{dim});
            y(1,:) = y(1,:)/sqrt(2);
            y = ifft(y,[],1);
            y = real(y*sqrt(2*n));
            I = (1:n)*0.5+0.5;
            I(2:2:end) = n-I(1:2:end-1)+1;
            y = y(I,:);
            y = reshape(y,siz);
            y = shiftdim(y,1);            
        end
    end
            
    y = reshape(y,sizy);

end

function Ims = findLedValues(Params, Ims)

    %if imNo == 1

    % Sample the Target Luminance
    leds = 3*sqrt(Ims.targetLumSm).*Params.ledMapSm;

    % Clip if the LED values exceed maximum LED value
    if max(max(leds))>Params.LED_MAX_LUM, 
        leds(leds>Params.LED_MAX_LUM)=Params.LED_MAX_LUM; 
    end

    % Find the luminance created by LEDs by taking convolution of LED
    % values and Point Spread Function (PSF). To reduce complexity both LED
    % values and PSF are their small version.
    ledLumSm = conv2(leds, Params.led_psfSm, 'same');

    % Update LED values by finding the ratio between the expected 
    % illumination created by L vector and the desired back illumination. 
    % Then, this scale has been multiplied with the LED values in order to 
    % get a much more closer approximation to the desired backlight values.

    % Scale the LEDs 
    Params.SAFETY_MARGIN = 30;
    scale = imdilate(Ims.targetLumSm + Params.SAFETY_MARGIN, ones(25,25))./ledLumSm;
    leds2 = (leds + 0.1).*scale.*spreadmap(size(Ims.targetLumSm),Params.led_psfSm);

    % Clip if the LED values exceed maximum LED value
    if max(max(leds2))>Params.LED_MAX_LUM, 
        leds2(leds2>Params.LED_MAX_LUM)=Params.LED_MAX_LUM; 
    end

    % Find the luminance created by LEDs
    ledLum2 = conv2(leds2,Params.led_psfSm, 'same');
    Ims.ledLumPrevSm = imfilter(ledLum2, fspecial('disk',5));
    
    tempVec = findledvals(leds2./Params.LED_MAX_LUM, Params.ledLabelsSm);
    Ims.ledValsPrev = tempVec(:,2);
    
    %end

    %% Iterative Update
    % Find the final LED values after an iterative scaling step
    
    % Initialize for iteration
    ledLumNew = Ims.ledLumPrevSm;
    ledLumOld = ledLumSm;
    ledsNew = zeros(size(Ims.ledLumPrevSm));
    ledsNew(Params.ledIndicesSm) = Ims.ledValsPrev;
    iterNum = 0;

    % Iterate until the error between consecutive luminance value maps are
    % decreased a certain threshold, decided as 100.
    while sum(sum((pu_encode(ledLumNew) - pu_encode(ledLumOld)).^2)) > 10
        scaleIt = Ims.targetLumSm./ledLumNew;
        
        % Scale the values and clip values under 0 and above MAX
        ledsNew = clipim(ledsNew.*scaleIt, Params.LED_MAX_LUM, 0); 
        
        % Find the luminance values wrt newly found LED values
        ledLumOld = ledLumNew;
        ledLumNew = conv2(ledsNew, Params.led_psfSm, 'same');
        iterNum = iterNum + 1;
    end
    
    % Find LED values after iteration
    tempVec = findledvals(ledsNew./Params.LED_MAX_LUM, Params.ledLabelsSm);
    ledsValsNew = ceil(tempVec(:,2).*255)./255;
    
    %% Check Power Constraint
    % Assume each LED consumes 1.4 Watts of power. Although the presentation
    % mentions 1500W power supply, undercalculate the power to be on the safe
    % side always.
    totalPower = sum(ledsValsNew*1.44);
    if(totalPower > 1400)
        disp('-!-Maximum power has been reached. LED values will be revised...');
        scalingFactor = totalPower/1400;
        
        % If maximum power is reached, scale it down
        ledsValsNewPow = ledsValsNew/scalingFactor;
        ledsValsNewPow = ceil(ledsValsNewPow.*255)./255;
    else
        ledsValsNewPow = ledsValsNew;
    end

    % Replace the LED values after scaling, and find the backlight
    % luminance
    ledsValsNewIm = placeleds(ledsValsNewPow, [1080 1920], Params.ledIndices);
    baseLum = conv2(ledsValsNewIm.*Params.LED_MAX_LUM, Params.led_psf, 'same');

    % Pass values
    %Ims.ledLumPrevSm = ledLumNew;
    Ims.ledVals = ledsValsNewPow;
    Ims.baseLum = baseLum;
end

function Ims = findLcdValues(Params, Ims)

    % If the extension is ".hdr" then multiply w/ constant
    if strcmp(Ims.fileXt,'.hdr')
        lcdLinear = cat(3, Params.LUMINANCE_HDR_CONSTANT*Ims.lastIm(:,:,1)./Ims.baseLum,...
                    cat(3, Params.LUMINANCE_HDR_CONSTANT*Ims.lastIm(:,:,2)./Ims.baseLum,...
                           Params.LUMINANCE_HDR_CONSTANT*Ims.lastIm(:,:,3)./Ims.baseLum));
    elseif strcmp(Ims.fileXt,'.exr')
        lcdLinear = cat(3, Ims.lastIm(:,:,1)./Ims.baseLum,... 
                    cat(3, Ims.lastIm(:,:,2)./Ims.baseLum,...
                           Ims.lastIm(:,:,3)./Ims.baseLum));
    end
    
    % Clip LCD values that exceed 1, and apply Gamma correction
    lcdLinear(lcdLinear>1) = 1;
    lcdLinear(lcdLinear<0) = 0;
    lcdGamCorrd = lcdgammacorr(lcdLinear);

    % Apply contrast enhancement step - Find histogram of LCD values, when
    % the cumulative value of these histogram bins reach pre-determined
    % threshold of << mean(histVals(:))/2 >> take that bin and stretch that
    % value to 0
    histVals = hist(lcdGamCorrd(:),500);
    cumHistVals = cumsum(histVals);
    minInd = find( cumHistVals > mean(histVals(:))/2 ,1)/500;
    lcdGamCorrd = imadjust(lcdGamCorrd, [minInd, minInd, minInd; 1 1 1], []);
    
    % Convert it to unit8 data type
    Ims.lcdGamCorrdUint8 = uint8(lcdGamCorrd*255);
end

function ledVals = findledvals(image, ledLabels)

    % Concatanate the image and ledLabels
    togetherArray = cat(3, ledLabels, image);
    togetherVector = reshape(togetherArray, [], 2);

    % Sort the rows to have LED IDs at the bottom
    ledValsTemp = sortrows(togetherVector,1);
    ledVals = ledValsTemp(end-2201:end,:);

end

function resArray = lcdgammacorr(inpVal, chanNum, valsArray)

    % If not manually entered, load predetermined values
    if ~exist('valsArray','var') 
        load('lcdGammaCorr.mat');
    end

    % Gamma correct image according to channel number
    if size(inpVal, 3) == 3
        resArray(:,:,1) = fnredval(inpVal(:,:,1), valsArray(1,:));
        resArray(:,:,2) = fngrnval(inpVal(:,:,2), valsArray(2,:));
        resArray(:,:,3) = fnbluval(inpVal(:,:,3), valsArray(3,:));
    elseif size(inpVal, 3) == 1
        if chanNum == 1
            resArray = fnredval(inpVal,valsArray(1,:));
        elseif chanNum == 2
            resArray = fngrnval(inpVal,valsArray(2,:));
        elseif chanNum == 3
            resArray = fnbluval(inpVal,valsArray(3,:));
        end
    else
        error('wrong operation in lcdGammaVals...');
    end

end

% Process Red Channel
function res = fnredval(inpX, valsRed)

    % Interpolate for intermediate values
    lowVal = valsRed(floor(inpX*255)+1); 
    upVal = valsRed(ceil(inpX*255)+1);
    powRed = lowVal + (inpX-floor(inpX)).*(upVal-lowVal);

    res = inpX.^(1./powRed);
end

% Process Green Channel
function res = fngrnval(inpX, valsGrn)

    % Interpolate for intermediate values
    lowVal = valsGrn(floor(inpX*255)+1); 
    upVal = valsGrn(ceil(inpX*255)+1);
    powGrn = lowVal + (inpX-floor(inpX)).*(upVal-lowVal);

    res = inpX.^(1./powGrn);
end

% Process Blue Channel
function res = fnbluval(inpX, valsBlu)

    % Interpolate for intermediate values
    lowVal = valsBlu(floor(inpX*255)+1); 
    upVal = valsBlu(ceil(inpX*255)+1);
    powBlu = lowVal + (inpX-floor(inpX)).*(upVal-lowVal);

    res = inpX.^(1./powBlu);
end

function resArray = lcdgammacorrinv(inpVal, chanNum, valsArray)

    if ~exist('valsArray','var') 
        load('lcdGammaCorr.mat');
        valsArray = ones(size(valsArray))./valsArray;
    end

    % Inverse Gamma correct image according to channel number
    if size(inpVal, 3) == 3
        resArray(:,:,1) = fnRedVal(inpVal(:,:,1),valsArray(1,:));
        resArray(:,:,2) = fnGrnVal(inpVal(:,:,2),valsArray(2,:));
        resArray(:,:,3) = fnBluVal(inpVal(:,:,3),valsArray(3,:));
    elseif size(inpVal, 3) == 1
        if chanNum == 1
            resArray = fnRedVal(inpVal,valsArray(1,:));
        elseif chanNum == 2
            resArray = fnGrnVal(inpVal,valsArray(2,:));
        elseif chanNum == 3
            resArray = fnBluVal(inpVal,valsArray(3,:));
        end
    else
        error('wrong operation in lcdGammaVals...');
    end

end

% Process Red Channel
function res = fnRedVal(inpX, valsRed)

    % Interpolate for intermediate values
    lowVal = valsRed(floor(inpX*255)+1); 
    upVal = valsRed(ceil(inpX*255)+1);
    powRed = lowVal + (inpX-floor(inpX)).*(upVal-lowVal);

    res = inpX.^(1./powRed);
end

% Process Green Channel
function res = fnGrnVal(inpX, valsGrn)

    % Interpolate for intermediate values
    lowVal = valsGrn(floor(inpX*255)+1); 
    upVal = valsGrn(ceil(inpX*255)+1);
    powGrn = lowVal + (inpX-floor(inpX)).*(upVal-lowVal);

    res = inpX.^(1./powGrn);
end

% Process Blue Channel
function res = fnBluVal(inpX, valsBlu)

    % Interpolate for intermediate values
    lowVal = valsBlu(floor(inpX*255)+1); 
    upVal = valsBlu(ceil(inpX*255)+1);
    powBlu = lowVal + (inpX-floor(inpX)).*(upVal-lowVal);

    res = inpX.^(1./powBlu);
end

function [ledsPlaced, ledIndices] = placeleds(ledArray, imSize, ledIndices)

    if(~exist('ledIndices', 'var'))
        % Find locations of LEDs
        [~, ~, ledIndices] = findLEDpositions(imSize);
    end

    % Create an empty image and place LED values
    ledsPlaced = zeros(imSize(1:2));
    ledsPlaced(ledIndices) = ledArray;

end

function impactMap = spreadmap(imSize, ledPsf, ledMax)

    if(~exist('ledMax', 'var')), ledMax = 235;   end;

    % Find luminance values if all the LEDs are turned on
    ledMap = findledpositions(imSize(1:2), 0);
    ledLum = conv2(ledMap*ledMax,ledPsf, 'same');

    % Divide all the luminance values found by the center pixel's value
    midVal = round(ledLum(imSize(1)/2,imSize(2)/2));
    impactMap = repmat(midVal, imSize(1:2))./ledLum;

end

function clippedIm = clipim(image, upLim, lowLim)

    % Check if upLim >= lowLim
    if upLim < lowLim, 
        error('Upper limit cannot be smaller than lower limit!!');
    end
        
    clippedIm = image;

    clippedIm(image>upLim) = upLim;
    clippedIm(image<lowLim) = lowLim;

end

function [hdrImDviP, ledLum3d] = dviplus2hdrim(ledVals, lcdVals, led_psf)

    if ~exist('led_psf','var'), 
        led_psf = lum(hdrimread('psf_model.pfm')); 
    end
    LED_MAX_LUM = 180;

    % Get LED Values
    [~, ~, ledIndices] = findledpositions(size(lcdVals), 0);

    % Find luminance map
    tempIm = zeros(size(lcdVals,1),size(lcdVals,2));
    tempIm(ledIndices) = double(ledVals)./255;
    ledLum = conv2(tempIm, led_psf*LED_MAX_LUM, 'same');
    ledLum3d = repmat(ledLum,[1 1 3]);

    % Take RGB Values
    tempIm = double(lcdVals);
    tempIm = tempIm/255;

    % Remove gamma-correction
    lcdRatio = lcdgammacorrinv(tempIm);

    % Normalize the image between 0.005-1 to simulate the light leakage caused 
    % by the LCD non-ideality.
    lcdRatio = normim(lcdRatio, 1, 0.005);

    % Find HDR Image by multiplying the luminance and LCD values
    hdrImDviP = ledLum3d.*lcdRatio;

end

function normalizedIm = normim(image, upLim, lowLim)
%NORMIM Normalize the given image within the range defined.
%
%   NORMALIZEDIM = NORMIM(IMAGE, UPLIM, LOWLIM) Normalizes the input IMAGE
%   within the range defined by UPLIM and LOWLIM. So that, the minimum
%   value of the NORMALIZEDIM becomes LOWLIM and the maximum value becomes
%   UPLIM.
%
% Examples:
%   normalizedIm = normim(image, 1, 0)
%   normalizedIm = normim(image, 1, 0.005)
%   normalizedIm = normim(image, 10, -10)
%
% ---------------------
% - Emin Zerman / emin.zerman@telecom-paristech.fr
% - Created:  04/05/2015
% - Telecom ParisTech - TSI - MM
% ---------------------

% If lowLim and upLim are not supplied, normalize between 0-1
if ~exist('upLim','var'), 
    upLim = 1; 
end
if ~exist('lowLim','var'), 
    lowLim = 0; 
end

% Find the scale and remove any offset
targetScale = upLim - lowLim;
imScale = max(image(:)) - min(image(:));

% Scale the image and normalize within range
normalizedIm = image*targetScale/imScale;
normalizedIm = normalizedIm - (min(normalizedIm(:)) - lowLim);

end

function dviPlusIm = getdviplusim(ledVals, lcdVals, fileName)

    %% Create DVI Plus header
    % Since the DVI+ image uses first two lines to pass the information about
    % LEDs, we have to write this information in the first two lines of the
    % image.
    newImage = zeros(size(lcdVals));

    % Write DVI+ header for display to understand this image is a special image
    dviPheader = [...
        158 158 158;
        23  23  23;
        211 211 211;
        36  36  36;
        198 198 198;
        53  53  53;
        187 187 187;
        79  79  79;
        0   0   0;
        0   0   0;
        59  36  0;
        127 71  55;
        0   0   0;
        8   153 7];

    newImage(1,1:14,:) = dviPheader;

    %% Get LED header
    ledVals = repmat(ledVals', 1, 3); % Write same LCD values on R, G, and B

    % Write LED values in the first two lines of the image
    newImage(1,15:end,:) = ledVals(1:1906,:);
    newImage(2, 1:296, :) = ledVals(1907:end,:);

    %% Set LCD values

    % Write LCD values to the rest of the image
    newImage(2, 297:end, :) = lcdVals(2, 297:end, :);
    newImage(3:end, :, :) = lcdVals(3:end, :, :);

    %% Pass the result (and save the Image)

    dviPlusIm = uint8(newImage);

    % If a fileName is specified, write the image to that file
    if exist('fileName', 'var')
        imwrite(dviPlusIm, fileName);
    end

end



















