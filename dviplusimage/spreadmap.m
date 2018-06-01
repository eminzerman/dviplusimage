function impactMap = spreadmap(imSize, ledPsf, ledMax)
%SPREADMAP Finds the spread effect on the SIM2 display.
%
%   IMPACTMAP = SPREADMAP(IMSIZE, LEDPSF, LEDMAX) Finds the light spread
%   effect on the SIM2 display. This effect reduces the amount of light at
%   the edges and the corners of the display. This function returns
%   IMPACTMAP which is the weight of how the LED values should be
%   multiplied to reduce this light spread effect. This function requires
%   IMSIZE, the image resolution, LEDPSF, the point spread function,
%   LEDMAX, the maximum luminance levels for each LED.
%
% Examples:
%   impactMap = spreadmap(size(hdrIm), ledPsf)
%
% ---------------------
% - Emin Zerman / emin.zerman@telecom-paristech.fr
% - Created:  05/03/2015
% - Telecom ParisTech - TSI - MM
% ---------------------

if(~exist('ledMax', 'var')), ledMax = 235;   end;

% Find luminance values if all the LEDs are turned on
ledMap = findledpositions(imSize(1:2), 0);
ledLum = conv2(ledMap*ledMax,ledPsf, 'same');

% Divide all the luminance values found by the center pixel's value
midVal = round(ledLum(imSize(1)/2,imSize(2)/2));
impactMap = repmat(midVal, imSize(1:2))./ledLum;

