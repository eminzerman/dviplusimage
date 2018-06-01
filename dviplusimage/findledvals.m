function ledVals = findledvals(image, ledLabels)
%FINDLEDVALS Finds the led values from the image of LEDs given.
%
%   LEDVALS = FINDLEDVALS(IMAGE, LEDLABELS) Finds the led values from the 
%   IMAGE of LEDs and LEDLABELS given. This function returns LED values as
%   a Nx2 array which contains LED IDs and LED values on each column of
%   LEDVALS.
%
% Examples:
%   ledVals = findledvals(image, ledLabels)
%
% ---------------------
% - Emin Zerman / emin.zerman@telecom-paristech.fr
% - Created:  23/02/2015
% - Telecom ParisTech - TSI - MM
% ---------------------

% Concatanate the image and ledLabels
togetherArray = cat(3, ledLabels, image);
togetherVector = reshape(togetherArray, [], 2);

% Sort the rows to have LED IDs at the bottom
ledValsTemp = sortrows(togetherVector,1);
ledVals = ledValsTemp(end-2201:end,:);

end