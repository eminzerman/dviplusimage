function clippedIm = clipim(image, upLim, lowLim)
%CLIPIM Clips the values of the image in the given threshold.
%
%   CLIPPEDIM = CLIPIM(IMAGE, UPLIM, LOWLIM) Clips the IMAGE values between
%   UPLIM and LOWLIM.
%
% Examples:
%   clippedIm = clipim(image, 255, 0)
%   clippedIm = clipim(image, 5, -5)
%
% ---------------------
% - Emin Zerman / emin.zerman@telecom-paristech.fr
% - Created:  20/03/2015
% - Telecom ParisTech - TSI - MM
% ---------------------

% Check if upLim >= lowLim
if upLim < lowLim, 
    error('Upper limit cannot be smaller than lower limit!!');
end
    
clippedIm = image;

clippedIm(image>upLim) = upLim;
clippedIm(image<lowLim) = lowLim;