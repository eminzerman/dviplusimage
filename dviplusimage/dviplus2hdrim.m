function hdrImDviP = dviplus2hdrim(ledVals, lcdVals, led_psf)
%DVIPLUS2HDRIM Converts DVI+ images to HDR images.
%
%   HDRIMDVIP = DVIPLUS2HDRIM(LEDVALS, LCDVALS) Converts DVI+ images to HDR
%   images. This function requires DVI+ image parameters such as LEDVALS
%   and LCDVALS, and outputs the HDR image as HDRIMDVIP.
%
%   HDRIMDVIP = DVIPLUS2HDRIM(LEDVALS, LCDVALS, LED_PSF) Does the same
%   operation. Instead of default PSF, it uses the given input LED_PSF.
%
% Examples:
%   hdrImDviP = dviplus2hdrim(ledValsNew, lcdGammCorrd)
%
% ---------------------
% - Emin Zerman / emin.zerman@telecom-paristech.fr
% - Created:  05/03/2015
% - Telecom ParisTech - TSI - MM
% ---------------------

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

% Take RGB Values
tempIm = double(lcdVals);
tempIm = tempIm/255;

% Remove gamma-correction
lcdRatio = lcdgammacorrinv(tempIm);

% Normalize the image between 0.005-1 to simulate the light leakage caused 
% by the LCD non-ideality.
lcdRatio = normim(lcdRatio, 1, 0.005);

% Find HDR Image by multiplying the luminance and LCD values
hdrImDviP = repmat(ledLum,[1 1 3]).*lcdRatio;

end






