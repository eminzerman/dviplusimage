function [hdrIm, lastIm, typeFlag, diff] = makeitfullhd(hdrIm)
%MAKEITFULLHD Changes the resolution of the input image to Full HD.
%
%   [HDRIM, LASTIM, TYPEFLAG, DIFF] = MAKEITFULLHD(HDRIM) Changes the
%   resolution of input image HDRIM and makes it Full HD resolution
%   (1920 x 1080) by image resizing and/or padding. The function returns
%   HDRIM, the image with repeated rows (or columns) to preserve color;
%   LASTIM, the image which is padded with black row (or column) pixels;
%   TYPEFLAG, the flag that indicates the the type of process done; DIFF is
%   the number of difference pixels between the expected row number and the
%   actual row number (or column).
%
% Examples:
%   [hdrOutIm, lastIm, typeFlag, diff] = makeitfullhd(hdrImage)
%
% ---------------------
% - Emin Zerman / emin.zerman@telecom-paristech.fr
% - Created:  01/06/2015
% - Telecom ParisTech - TSI - MM
% ---------------------

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