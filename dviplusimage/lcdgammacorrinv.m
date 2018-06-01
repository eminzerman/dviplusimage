function resArray = lcdgammacorrinv(inpVal, chanNum, valsArray)
%LCDGAMMACORRINV Inverts Gamma correction operation.
%
%   RESARRAY = LCDGAMMACORRINV(INPVAL, CHANNUM, VALSARRAY) Takes Gamma
%   corrected LCD image and applies an inverse Gamma correction process. 
%   INPVAL is the input gamma-corrected LCD image where CHANNUM is the 
%   channel number an it is relevant only if INPVAL dimensions are less 
%   than 3. VALSARRAY has the parameters for Gamma correction. If known, 
%   these parameters can be supplied.
%
% Examples:
%   resArray = lcdgammacorrinv(lcdGammaCorr)
%
% ---------------------
% - Emin Zerman / emin.zerman@telecom-paristech.fr
% - Created:  20/04/2015  -- Updated: 22/04/2015
% - Telecom ParisTech - TSI - MM
% ---------------------

if ~exist('valsArray','var') 
    load('lcdGammaCorr.mat');
    valsArray = ones(size(valsArray))./valsArray;
end

% Inverse Gamma correct image according to channel number
if size(inpVal, 3) == 3
    resArray(:,:,1) = fnRedVal(inpVal(:,:,1),valsArray(1,:));
    resArray(:,:,2) = fnGrnVal(inpVal(:,:,2),valsArray(2,:));
    resArray(:,:,3) = fnBluVal(inpVal(:,:,3),valsArray(3,:));
elseif size(inpVal, 3) == 1
    if chanNum == 1
        resArray = fnRedVal(inpVal,valsArray(1,:));
    elseif chanNum == 2
        resArray = fnGrnVal(inpVal,valsArray(2,:));
    elseif chanNum == 3
        resArray = fnBluVal(inpVal,valsArray(3,:));
    end
else
    error('wrong operation in lcdGammaVals...');
end

end

% Process Red Channel
function res = fnRedVal(inpX, valsRed)

    % Interpolate for intermediate values
    lowVal = valsRed(floor(inpX*255)+1); 
    upVal = valsRed(ceil(inpX*255)+1);
    powRed = lowVal + (inpX-floor(inpX)).*(upVal-lowVal);

    res = inpX.^(1./powRed);
end

% Process Green Channel
function res = fnGrnVal(inpX, valsGrn)

    % Interpolate for intermediate values
    lowVal = valsGrn(floor(inpX*255)+1); 
    upVal = valsGrn(ceil(inpX*255)+1);
    powGrn = lowVal + (inpX-floor(inpX)).*(upVal-lowVal);

    res = inpX.^(1./powGrn);
end

% Process Blue Channel
function res = fnBluVal(inpX, valsBlu)

    % Interpolate for intermediate values
    lowVal = valsBlu(floor(inpX*255)+1); 
    upVal = valsBlu(ceil(inpX*255)+1);
    powBlu = lowVal + (inpX-floor(inpX)).*(upVal-lowVal);

    res = inpX.^(1./powBlu);
end