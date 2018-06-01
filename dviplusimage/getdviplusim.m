function dviPlusIm = getdviplusim(ledVals, lcdVals, fileName)
%GETDVIPLUSIM Creates a DVI+ image givne parameters.
%
%   GETDVIPLUSIM(LEDVALS, LCDVALS, FILENAME) Creates a DVI+ image with the
%   given FILENAME using LEDVALS and LCDVALS. LEDVALS should be a 2202x1
%   matrix and LCDVALS should be 1920x1080 uint8 RGB image.
%
%   DVIPLUSIM = GETDVIPLUSIM(LEDVALS, LCDVALS) Creates a DVI+ image called
%   DVIPLUSIM in MATLAB workspace, and does not create a file.
%
%   DVIPLUSIM = GETDVIPLUSIM(LEDVALS, LCDVALS, FILENAME) Creates both an
%   image in MATLAB workspace and writes it as a file on the computer.
%
% Examples:
%   getdviplusim(ledVals, lcdVals, 'myDviPlusImage.png')
%   dviPlusIm = getdviplusim(ledVals, lcdVals)
%   dviPlusIm = getdviplusim(ledVals, lcdVals, 'myDviPlusImage.png')
%
% ---------------------
% - Emin Zerman / emin.zerman@telecom-paristech.fr
% - Created:  20/02/2015
% - Telecom ParisTech - TSI - MM
% ---------------------

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



