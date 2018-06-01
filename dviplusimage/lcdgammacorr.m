function resArray = lcdgammacorr(inpVal, chanNum, valsArray)
%LCDGAMMACORR Corrects LCD values according to Gamma correction curves.
%
%   RESARRAY = LCDGAMMACORR(INPVAL, CHANNUM, VALSARRAY) Takes the LCD image
%   and applies a Gamma correction process. INPVAL is the input LCD image
%   where CHANNUM is the channel number an it is relevant only if INPVAL
%   dimensions are less than 3. VALSARRAY has the parameters for Gamma
%   correction. If known, these parameters can be supplied.
%
% Examples:
%   resArray = lcdgammacorr(hdrIm)
%
% ---------------------
% - Emin Zerman / emin.zerman@telecom-paristech.fr
% - Created:  20/04/2015  -- Updated: 22/04/2015
% - Telecom ParisTech - TSI - MM
% ---------------------

% If not manually entered, load predetermined values
if ~exist('valsArray','var') 
    load('lcdGammaCorr.mat');
end

% Gamma correct image according to channel number
if size(inpVal, 3) == 3
    resArray(:,:,1) = fnredval(inpVal(:,:,1), valsArray(1,:));
    resArray(:,:,2) = fngrnval(inpVal(:,:,2), valsArray(2,:));
    resArray(:,:,3) = fnbluval(inpVal(:,:,3), valsArray(3,:));
elseif size(inpVal, 3) == 1
    if chanNum == 1
        resArray = fnredval(inpVal,valsArray(1,:));
    elseif chanNum == 2
        resArray = fngrnval(inpVal,valsArray(2,:));
    elseif chanNum == 3
        resArray = fnbluval(inpVal,valsArray(3,:));
    end
else
    error('wrong operation in lcdGammaVals...');
end

end

% Process Red Channel
function res = fnredval(inpX, valsRed)

    % Interpolate for intermediate values
    lowVal = valsRed(floor(inpX*255)+1); 
    upVal = valsRed(ceil(inpX*255)+1);
    powRed = lowVal + (inpX-floor(inpX)).*(upVal-lowVal);

    res = inpX.^(1./powRed);
end

% Process Green Channel
function res = fngrnval(inpX, valsGrn)

    % Interpolate for intermediate values
    lowVal = valsGrn(floor(inpX*255)+1); 
    upVal = valsGrn(ceil(inpX*255)+1);
    powGrn = lowVal + (inpX-floor(inpX)).*(upVal-lowVal);

    res = inpX.^(1./powGrn);
end

% Process Blue Channel
function res = fnbluval(inpX, valsBlu)

    % Interpolate for intermediate values
    lowVal = valsBlu(floor(inpX*255)+1); 
    upVal = valsBlu(ceil(inpX*255)+1);
    powBlu = lowVal + (inpX-floor(inpX)).*(upVal-lowVal);

    res = inpX.^(1./powBlu);
end