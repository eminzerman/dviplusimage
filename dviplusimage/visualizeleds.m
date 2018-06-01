function ledVisImage = visualizeleds(ledLabels, ledArray, debugFlag, strWidth)
%VISUALIZELEDS Creates an RGB image to visualize the LED values.
%
%   LEDVISIMAGE = VISUALIZELEDS(LEDLABELS, LEDARRAY, DEBUGFLAG, STRWIDTH)
%   Creates an RGB image to visualize the LED values. This function returns
%   an RGB image which shows each LED as brightened by the indicated
%   amount. This function requires LEDLABELS, the image that indicates LED
%   IDs, LEDARRAY, the LED values represented in an array of 2202x1,
%   DEBUGFLAG, the flag to show or hide the debug messages/figures,
%   STRWIDTH, the width of structuring element.
%
% Examples:
%   ledVisImage = visualizeleds(ledLabels, ledArray)
%   ledVisImage = visualizeleds(ledLabels, ledArray, 0, 12)
%
% ---------------------
% - Emin Zerman / emin.zerman@telecom-paristech.fr
% - Created:  02/03/2015
% - Telecom ParisTech - TSI - MM
% ---------------------

if(~exist('debugFlag', 'var')),  debugFlag = 0;   end
if(~exist('strWidth', 'var')),   strWidth = 12;   end

% Assign zero LEDs to a very small value
ledArray(ledArray==0) = 1e-6;

% Assign LED values to their locations
tempX = zeros(size(ledLabels))';
tempX(ledLabels'~=0) = ledArray;
tempX = tempX';

% Dilate the LED values by an octagon structuring element
ledVisImage = imdilate(tempX, strel('octagon',strWidth));
ledVisImage(ledVisImage == 0) = -0.05;

% Create a border to differentiate the black LEDs from the frame edges
ledVisImage = ledVisImage - min(ledVisImage(:));
ledVisImage = ledVisImage./max(ledVisImage(:));
if debugFlag, figure, imshow(ledVisImage); end;

end