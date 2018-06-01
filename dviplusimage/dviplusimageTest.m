function dviPlusIm = dviplusimageTest(fileFolderName, debugFlag, testPhrase, outFolder)
%DVIPLUSIMAGE Creates DVI Plus image for SIM2 HDR47 screen.
%   DVIPLUSIM = DVIPLUSIMAGE(FILEFOLDERNAME) Creates DVI Plus image and 
%   returns the image as DVIPLUSIM. The function requires FILEFOLDERNAME
%   as input. This input can both be a file or directory.
%
%   DVIPLUSIMAGE(FILEFOLDERNAME, DEBUGFLAG) To see different level of 
%   information, DEBUGFLAG can be used. It may take a value between 0-4. 
%   DEBUGFLAG 0 gives no information, 4 give all information available.
%
%   DVIPLUSIMAGE(FILEFOLDERNAME, DEBUGFLAG, TESTPHRASE, OUTFOLDER) The 
%   input TESTPHRASE can be used to differentiate between different test 
%   batches. OUTFOLDER specifies the output folder of DVI Plus images.
%
% Examples:
%   dvipImage = dviplusimage('hdrImages\cloudyDay.hdr',3);
%   dviplusimage('hdrTestFolder', 1, 'lumVersion2', 'hdrOutFolder');
%   dviplusimage('hdrImages', 0, '', 'outFolder');
%   dviplusimage('hdrImages');
%
% Debug Levels:
%   0: No debug messages/images
%   1: Save luminance values in ".exr" format
%   2: Save images in ".exr" format at multiple levels
%   3: Both save and show images for some steps (only for files)
%   4: Both save and show images for each step (only for files)
%
% ---------------------
% - Emin Zerman / emin.zerman@telecom-paristech.fr
% - Created:  23/03/2015 -- Updated: 16/09/2015
% - Telecom ParisTech - TSI - MM
% ---------------------

% Check if flags exist and initialize
if(~exist('debugFlag', 'var')),     debugFlag = 0;  end
if(~exist('testPhrase', 'var')),    testPhrase =''; end
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
end

% CONSTANT values
LED_MAX_LUM = 180;
LUMINANCE_HDR_CONSTANT = 179;
MAX_LUMINANCE = 4250;
EXTRAPOLATE_MARGIN = 80;
LCD_LEAKAGE_CONSTANT = 1/0.005;
fileFlag = 0;

% Check if the first input is a file name or folder name
if ~isempty(strfind(fileFolderName(end-4:end),'.hdr')) || ...
        ~isempty(strfind(fileFolderName(end-4:end),'.exr'))
    % ========== It is an HDR or EXR file ==========
    fileFlag = 1;
    fileIndices = 1;
    DbContents.name = fileFolderName;
elseif isdir(fileFolderName)
    % ========== It is a folder ==========
    fileFlag = 0;
    % Find folder contents
    DbContents = dir(fileFolderName);
    fileIndices = 3:length(DbContents);
else
    error('This is not an HDR image file or a folder!');
end

% Process the indicated image or all images in the indicated folder
for imNo = fileIndices
    close all;
    
    % Find fileName for naming
    if fileFlag
        [~, fileName, fileXt] = fileparts(fileFolderName);
    else
        [~, fileName, fileXt] = fileparts(DbContents(imNo).name);
    end
    
    %% Acquire Image
    % Acquire HDR image by utilizing Matlab's hdrread function, and tonemap
    % function for user to see the image.
    if fileFlag
        if strcmp(fileXt,'.hdr')
            try
                hdrIm = hdrread(DbContents(imNo).name); 
            catch
                try
                    hdrIm = hdrimread(DbContents(imNo).name);
                catch
                    error('Problem acquiring the HDR image/frame!!');
                end
            end
        elseif strcmp(fileXt,'.exr')
            hdrIm = exrread(DbContents(imNo).name);
        else
            continue;
        end
    else
        if strcmp(fileXt,'.hdr')
            try
                hdrIm = hdrread([fileFolderName filesep DbContents(imNo).name]); 
            catch
                try
                    hdrIm = hdrimread([fileFolderName filesep DbContents(imNo).name]);
                catch
                    error('Problem acquiring the HDR image/frame!!');
                end
            end
        elseif strcmp(fileXt,'.exr')
            hdrIm = exrread([fileFolderName filesep DbContents(imNo).name]);
        else
            continue;
        end
    end
    
    % Adjust image to be 1920/1080 by resizing and/or padding with black
    [hdrIm, lastIm] = makeitfullhd(hdrIm);
%     hdrIm = hdrIm./2;
    
    % Find target luminance, multiply w/ constant if ".hdr"
    if strcmp(fileXt,'.hdr')
        hdrImL = hdrIm*LUMINANCE_HDR_CONSTANT;
        lastImL = lastIm*LUMINANCE_HDR_CONSTANT;
    elseif strcmp(fileXt,'.exr')
        hdrImL = hdrIm;
        lastImL = lastIm;
    end
    hdrLum = max(hdrImL, [], 3);
    
    disp(['==== ' DbContents(imNo).name ' started ====']);

    if debugFlag>=1, exrwrite(lastImL, [outFolder filesep fileName '_00_hdrImage.exr' ]);  end;
    if debugFlag>=1, imwrite(ReinhardBilTMO(double(lastImL)),[outFolder filesep fileName '_01_tonemap.png' ]);  end;
    if debugFlag>=3 && fileFlag, figure, imshow(tonemap(hdrImL, 'AdjustSaturation', 3)); title(['Tonemapped version of ' fileName]); end;

    % Get point spread function, resize PSF to reduce complexity
    led_psf = lum(hdrimread('psf_model.pfm'));
    led_psfSm = imresize(led_psf, 1/4); 

    % Find locations of LEDs, resize LED maps to reduce complexity
    [~, ledLabels, ledIndices] = findledpositions(size(hdrImL), 0); 
    [ledMapSm, ledLabelsSm, ledIndicesSm] = findledpositions(size(hdrImL)/4, 0);

    %% Find desired backlight
    
    % If the luminance is too high, clip
    if max(max(hdrLum))>MAX_LUMINANCE, 
        hdrLum(hdrLum>MAX_LUMINANCE) = MAX_LUMINANCE; 
        disp('-!-hdrLum values are over MAX! Clipping...'); 
    end

    % Extrapolate (by inpainting) the image to avoid information loss 
    % while filtering. Downsample the hdrLum by 4 to reduce complexity.
    hdrLumBar = zeros( (EXTRAPOLATE_MARGIN+1080)/4,...
                       (EXTRAPOLATE_MARGIN+1920)/4 ); 
    hdrLumBar(:) = NaN;
    hdrLumBar(EXTRAPOLATE_MARGIN/8 +1 : 1080/4 +EXTRAPOLATE_MARGIN/8,...
              EXTRAPOLATE_MARGIN/8 +1 : 1920/4 +EXTRAPOLATE_MARGIN/8) ...
              = imresize(hdrLum, 1/4, 'nearest');
    hdrLumBarN = inpaintn(hdrLumBar);
    hdrLumExtp = imresize(hdrLumBarN, 4, 'nearest');
    hdrLumExtp(EXTRAPOLATE_MARGIN/2 +1 : 1080 +EXTRAPOLATE_MARGIN/2,...
                  EXTRAPOLATE_MARGIN/2 +1 : 1920 +EXTRAPOLATE_MARGIN/2) ...
                  = hdrLum;
    
    % -- Find Target Luminance value --
    % Find the minimum required luminance for fidelity by filtering the max
    % value image with a Gaussian filter (i.e. fspecial('disk',N) )
    maxValImg = imdilate(hdrLumExtp*1.1, strel('disk', 30));
    targetLumMin = conv2(maxValImg, fspecial('disk',30),'same');
    
    % Find the maximum allowed luminance to avoid any LCD leakage and light
    % halos that may occur. Median filter in order to reduce the artifacts
    % created due to the high frequency characteristics (e.g. stars, pixel
    % defects caused during image acquisition, etc.)
    targetLumMax = hdrLumExtp*LCD_LEAKAGE_CONSTANT; 
    targetLumMax(targetLumMax>MAX_LUMINANCE) = MAX_LUMINANCE; 
    targetLumMax = medfilt2(targetLumMax, [5 5]);
    
    % Find the minumum value (or the constraining value) of these minumum 
    % and maximum images.
    tempArr = cat(3, targetLumMin, targetLumMax);
    targetLum_pre = min(tempArr, [], 3); 
    
    % Smooth out the preliminary result to avoid direct faults of 
    % individual LEDs. Median filter in order to reduce the artifacts
    % created due to the high frequency characteristics (e.g. stars, pixel
    % defects caused during image acquisition, etc.)
    targetLum_filt = imfilter(targetLum_pre, fspecial('disk', EXTRAPOLATE_MARGIN));
    targetLum_temp = targetLum_filt(EXTRAPOLATE_MARGIN/2 +1 : 1080 +EXTRAPOLATE_MARGIN/2,...
                                    EXTRAPOLATE_MARGIN/2 +1 : 1920 +EXTRAPOLATE_MARGIN/2);
	targetLum = medfilt2(max(hdrLum,targetLum_temp), [5 5]); 
    targetLum(targetLum==0) = targetLum_temp(targetLum==0);
    targetLumSm2 = imresize(targetLum, 1/4, 'nearest'); 
    targetLumSm2(targetLumSm2 < 0) = 0;
    targetLumSm3 = imfilter(targetLumSm2,fspecial('disk',17)); %fimagesc(targetLumSm)
    %IMDILATE
%     targetLumSm = imdilate(targetLumSm2,strel('disk',5));   %fimagesc(targetLumSm)
    %NORM
    targetLumSm = normim(targetLumSm3, max(targetLumSm2(:)), min(targetLumSm2(:)));
    

    if debugFlag>=1, exrwrite(targetLumMin, [outFolder filesep fileName '_02_targetLumMin.exr' ]);  end;
    if debugFlag>=1, exrwrite(targetLumMax, [outFolder filesep fileName '_03_targetLumMax.exr' ]);  end;
    if debugFlag>=1, exrwrite(targetLum, [outFolder filesep fileName '_04_targetLuminance.exr' ]);  end;
    if debugFlag>=4 && fileFlag, figure, imagesc(targetLum); title('Target Luminance'); axis image; end;
    if debugFlag>=4 && fileFlag, figure, imagesc(targetLum - hdrLum); title('Target Luminance - hdrLum'); axis image; end;

    % Sample the Target Luminance
    leds = 3*sqrt(targetLumSm).*ledMapSm;

    % Clip if the LED values exceed maximum LED value
    if max(max(leds))>LED_MAX_LUM, 
        leds(leds>LED_MAX_LUM)=LED_MAX_LUM; 
        disp('-!-Led values are over MAX! Clipping...'); 
    end

    % Find the luminance created by LEDs by taking convolution of LED
    % values and Point Spread Function (PSF). To reduce complexity both LED
    % values and PSF are their small version.
    ledLumSm = conv2(leds, led_psfSm, 'same');

    if debugFlag>=1, exrwrite(imresize(ledLumSm, 4), [outFolder filesep fileName '_05_luminanceAfterSampling.exr' ]);  end;
    if debugFlag>=3 && fileFlag, figure, imagesc(hdrLum); title('hdrLum'); axis image; end;
    if debugFlag>=3 && fileFlag, figure, imagesc(ledLumSm); title('ledLum'); axis image; end;
    if debugFlag>=3 && fileFlag, figure, imagesc(ledLumSm - hdrLum); title('ledLum - hdrLum'); axis image; end;
    if debugFlag>=3 && fileFlag, figure, imagesc(ledLumSm - targetLum); title('ledLum - targetLum'); axis image; end;

    %% Update LED values
    % Update LED values by finding the ratio between the expected 
    % illumination created by L vector and the desired back illumination. 
    % Then, this scale has been multiplied with the LED values in order to 
    % get a much more closer approximation to the desired backlight values.

    % Scale the LEDs 
%     SAFETY_MARGIN = 30;
%     scale = imdilate(targetLumSm + SAFETY_MARGIN, ones(25,25))./(ledLumSm+ 1e-8);
%     leds2 = (leds + 0.1).*scale.*spreadmap(size(targetLumSm),led_psfSm);
    leds2 = leds;

    % Clip if the LED values exceed maximum LED value
    if max(max(leds2))>LED_MAX_LUM, 
        leds2(leds2>LED_MAX_LUM)=LED_MAX_LUM; 
        disp('-!-Led2 values are over MAX! Clipping...'); 
    end

    % Find the luminance created by LEDs
    ledLum2 = conv2(leds2,led_psfSm, 'same');
    ledLum3 = imfilter(ledLum2, fspecial('disk',5));

    if debugFlag>=1, exrwrite(imresize(scale, 4), [outFolder filesep fileName '_06_scaleAfterSampling.exr' ]);  end;
    if debugFlag>=1, exrwrite(imresize(ledLum3, 4), [outFolder filesep fileName '_07_luminanceFirstScaling.exr' ]);  end;
    if debugFlag>=1 && fileFlag, figure, imagesc(ledLum3); title('ledLum2'); axis image; end;
    if debugFlag>=3 && fileFlag, figure, imagesc(ledLum3 - hdrLum); title('ledLum2 - hdrLum'); axis image; end;
    if debugFlag>=3 && fileFlag, figure, imagesc(ledLum3 - targetLum); title('ledLum2 - targetLum'); axis image; end    
    
    tempVec = findledvals(leds2./LED_MAX_LUM, ledLabelsSm);
    ledVals2 = tempVec(:,2);

    %% Iterative Update
    % Find the final LED values after an iterative scaling step
    
    % Initialize for iteration
    ledLumNew = ledLum3;
    ledLumOld = ledLumSm;
    ledsNew = zeros(size(ledLum3));
    ledsNew(ledIndicesSm) = ledVals2;
    iterNum = 0;

    % Iterate until the error between consecutive luminance value maps are
    % decreased a certain threshold, decided as 100.
    while sum(sum((pu_encode(ledLumNew) - pu_encode(ledLumOld)).^2)) > 10
        scaleIt = targetLumSm./(ledLumNew + 1e-8);
        
        % Scale the values and clip values under 0 and above MAX
        ledsNew = clipim(ledsNew.*scaleIt, LED_MAX_LUM, 0); 
        
        % Find the luminance values wrt newly found LED values
        ledLumOld = ledLumNew;
        ledLumNew = conv2(ledsNew, led_psfSm, 'same');
        iterNum = iterNum + 1;
        
        if debugFlag>=2 && iterNum <= 3, exrwrite(imresize(scaleIt, 4), [outFolder filesep fileName '_08_scaleIteration' num2str(iterNum, '%02d') '.exr' ]);  end;
        if debugFlag>=2 && iterNum <= 3, exrwrite(imresize(ledLumNew, 4), [outFolder filesep fileName '_09_luminanceIteration' num2str(iterNum, '%02d') '.exr' ]);  end;
        if debugFlag>=2, disp(['iter ' num2str(iterNum) ' - error: ' num2str(sum(sum((pu_encode(ledLumNew) - pu_encode(ledLumOld)).^2))) ' - errorLum: ' num2str(sum(sum((ledLumNew - targetLumSm).^2)))]); end
    end
    
    % Find LED values after iteration
    tempVec = findledvals(ledsNew./LED_MAX_LUM, ledLabelsSm);
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
    ledsValsNewIm = placeleds(ledsValsNewPow, [1080 1920], ledIndices);
    baseLum = conv2(ledsValsNewIm.*LED_MAX_LUM, led_psf, 'same');
    ledLumNew = (baseLum + 1e-8);
    
    if debugFlag>=2, exrwrite(ledLumNew, [outFolder filesep '10_' fileName '_luminanceBacklight.exr' ]);  end;
    if debugFlag>=2, imwrite(uint8(normim(ledLumNew)*255), hot(256), [outFolder filesep fileName '_10_luminanceBacklight.png' ]);  end;

    %% Draw LED values
%     if debugFlag>=2, 
        ledVis = visualizeleds(ledLabels, ledsValsNewPow);
        imwrite(uint8(ledVis*255), hot(256), [outFolder filesep fileName '_11_ledValues.png' ]);  
        imwrite(imresize(imfuse(targetLumSm,ledMapSm),4), [outFolder filesep fileName '_11_ledCorrespondance.png' ]);
%     end;
    if debugFlag>=4 && fileFlag,   fimshow(ledVis);   end;

    %% Finding LCD values
    % After LED values and backlight illumination caused by LEDs have been
    % found, calculate and set LCD values for DVI Plus image.

    % If the extension is ".hdr" then multiply w/ constant
    if strcmp(fileXt,'.hdr')
        lcdLinear = cat(3, LUMINANCE_HDR_CONSTANT*lastIm(:,:,1)./ledLumNew,...
                    cat(3, LUMINANCE_HDR_CONSTANT*lastIm(:,:,2)./ledLumNew,...
                           LUMINANCE_HDR_CONSTANT*lastIm(:,:,3)./ledLumNew));
    elseif strcmp(fileXt,'.exr')
        lcdLinear = cat(3, lastIm(:,:,1)./ledLumNew,... 
                    cat(3, lastIm(:,:,2)./ledLumNew,...
                           lastIm(:,:,3)./ledLumNew));
    end
    
    % Clip LCD values that exceed 1, and apply Gamma correction
    lcdLinear(lcdLinear>1) = 1;
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
    lcdGamCorrdUint8 = uint8(lcdGamCorrd*255);

    %% Save DVI+ output image
        
    % Save DVI+ Image
    if debugFlag>=2, 
        dviPlusIm = getdviplusim(round(255*ledsValsNewPow)',...
                                 uint8(lcdGamCorrdUint8),...
                                 [outFolder filesep fileName '_' date testPhrase '.png']);
        
        % --Find and show the error in both EXR and PNG format--
        % Find the reconstructed HDR image
        dvipHdr = dviplus2hdrim(round(255*ledsValsNewPow)',...
                                uint8(lcdGamCorrdUint8), led_psf);
        
        hdrErrorIm = hdrImL - dvipHdr;
        exrwrite(dvipHdr, [outFolder filesep fileName '_12_hdrReconstructed.exr' ]);
        exrwrite(hdrErrorIm, [outFolder filesep fileName '_13_hdrErrorImage.exr' ]);
        
        % Find maximum amount of error
        devErr = max(abs(min(hdrErrorIm(:))), abs(max(hdrErrorIm(:))));
        hdrErrorIm(1,1) = -devErr; hdrErrorIm(1,2) = devErr;
        lcdErrorIm = normim(hdrErrorIm, 255, 0);
        imwrite(uint8(lcdErrorIm), [outFolder filesep fileName '_14_lcdErrorImage_' date testPhrase '.png']); 
    else
        dviPlusIm = getdviplusim(round(255*ledsValsNewPow)',...
                                 uint8(lcdGamCorrdUint8),...
                                 [outFolder filesep fileName '.png']);
                                                  
        % Find the reconstructed HDR image
%         dvipHdr = dviplus2hdrim(round(255*ledsValsNewPow)',...
%                                 uint8(lcdGamCorrdUint8), led_psf);
%         exrwrite(dvipHdr, [outFolder filesep fileName '_12_hdrReconstructed.exr' ]);
    end
    
end





















