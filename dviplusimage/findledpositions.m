function [ledImage, ledLabels, ledIndices, ledPos] = findledpositions(imSize, debugFlag)
%FINDLEDPOSITIONS Finds LED positions for SIM2 HDR47S4 display.
% 
%   [LEDIMAGE, LEDLABELS, LEDINDICES, LEDPOS] = FINDLEDPOSITIONS(IMSIZE, 
%   DEBUGFLAG) Finds the LED positions for SIM2 HDR47S4 display. This
%   function returns LEDIMAGE, a black image with white pixels on the
%   locations of the LEDs; LEDLABELS, corresponding labels for those LEDs;
%   LEDINDICES, the indice numbers for LEDs; LEDPOS, the array that
%   indicates the positions of LEDs wrt origin. This function requires
%   IMSIZE, the image resolution, and DEBUGFLAG, the flag which will be
%   used to show or hide debug messages/figures.
%
% Examples:
%   [ledImage, ledLabels] = findledpositions([1080 1920 3])
%   [ledImage, ledLabels, ~, ledPos] = findledpositions([1080 1920 3], 1)
%
% ---------------------
% - Emin Zerman / emin.zerman@telecom-paristech.fr
% - Created:  19/02/2015
% - Telecom ParisTech - TSI - MM
% ---------------------

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