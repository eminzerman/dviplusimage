function [ledsPlaced, ledIndices] = placeleds(ledArray, imSize, ledIndices)
%PLACELEDS Replaces the found LED values in image for backlight creation.
%
%   [LEDSPLACED, LEDINDICES] = PLACELEDS(LEDARRAY, IMSIZE) Replaces the LED
%   values from LEDARRAY to the image LEDSPLACED. IMSIZE is required to
%   find the LEDINDICES and to replace LED values correctly. The found
%   LEDINDICES values are also given as an output of this function.
%
%   LEDSPLACED = PLACELEDS(LEDARRAY, IMSIZE, LEDINDICES) Does the same
%   operation. If known, the LEDINDICES values can be given to this 
%   function to reduce computation time.
%
% Examples:
%   ledsPlaced = placeleds(ledArray, imSize, ledIndices)
%   [ledsPlaced, ledIndices] = placeleds(ledArray, imSize)
%
% ---------------------
% - Emin Zerman / emin.zerman@telecom-paristech.fr
% - Created:  23/03/2015
% - Telecom ParisTech - TSI - MM
% ---------------------

	if(~exist('ledIndices', 'var'))
	    % Find locations of LEDs
		[~, ~, ledIndices] = findLEDpositions(imSize);
    end

    % Create an empty image and place LED values
	ledsPlaced = zeros(imSize(1:2));
	ledsPlaced(ledIndices) = ledArray;

end