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

