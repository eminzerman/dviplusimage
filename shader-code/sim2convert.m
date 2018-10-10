function outIm = sim2convert(im)
%SIM2CONVERT Creates a LogLuv image for built-in rendering of SIM2 HDR47.
%   OUTIM = SIM2CONVERT(IM) Creates a LogLuv image for built-in rendering
%   of SIM2 HDR47 display and returns OUTIM. This OUTIM can be saved as 
%   a bitmap or PNG image and displayed on a SIM2 HDR display.
% 
%   This code is the Matlab implementation of the shader presented in
%   EasyHDRPlayer software, provided with SIM2 display.
%
% Example:
%   hdrIm = hdrimread('hdrImages\cloudyDay.hdr');
%   outIm = sim2convert(hdrInpIm);
%   imwrite(outIm, 'hdrImages\cloudyDay.png');
%
% ---------------------
% - Emin Zerman / emin.zerman@telecom-paristech.fr
% - Created:  27/04/2016
% - Telecom ParisTech - TSI - MM
% ---------------------
%  Copyright (C) 2018 - Emin Zerman
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------

% Check if image size is Full HD
if size(im,1) ~= 1080 || size(im,2) ~= 1920
	error('Image should be Full HD!!');
end

% Init some parameters
im = double(im);
alphaVal = 0.0376;
Lscale = 32;

% Color conversion
imXyz = rgb2xyz_sim2(im);
imYuv = xyz2yuv_m(imXyz);

% Take Y channel
imY = imYuv(:,:,1);

% Process for u and v channel
imYuvFilt = imfilter(imYuv, [1/3 1/3 1/3]);
imYuvFilt(:,1,:) = (imYuv(:,1,:) + imYuv(:,2,:))/2;
imYuvFilt(:,1920,:) = (imYuv(:,1919,:) + imYuv(:,1920,:))/2;
imU = imYuvFilt(:,:,2);
imV = imYuvFilt(:,:,3);

% Find L channel
idxMinY = imY < 0.00001;
idxMaxY = imY >= 0.00001;
imLtemp = zeros(size(imY));
imLtemp(idxMinY) = 0;
imLtemp(idxMaxY) = (alphaVal.*log2(imY(idxMaxY))) + 0.5;

imL = (253.*imLtemp + 1).*Lscale;
imL = floor(imL + 0.5);
imL = imClamp(imL, 32, 8159);
imLh= floor(imL./32);
imLl= imL - (imLh.*32);

% oddCase
imV_odd = floor(imV + 0.5);
imV_odd = imClamp(imV_odd, 4, 1019);
imCh_odd = floor(imV_odd./4);
imCl_odd = imV_odd - (imCh_odd.*4);

outR_odd = (imLl.*8) + (imCl_odd.*2);
outR_odd = imClamp(outR_odd, 1, 254);
outG_odd = imLh;
outB_odd = imCh_odd;

% evenCase
imU_even = floor(imU + 0.5);
imU_even = imClamp(imU_even, 4, 1019);
imCh_even = floor(imU_even./4);
imCl_even = imU_even - (imCh_even.*4);

outR_even = imCh_even;
outG_even = imLh;
outB_even = (imLl.*8) + (imCl_even.*2);
outB_even = imClamp(outB_even, 1, 254);

% Prepare RGB image
outRGB = zeros(size(im));
outRGB(:,1:2:end,1) = outR_odd(:,1:2:end); outRGB(:,2:2:end,1) = outR_even(:,2:2:end);
outRGB(:,1:2:end,2) = outG_odd(:,1:2:end); outRGB(:,2:2:end,2) = outG_even(:,2:2:end);
outRGB(:,1:2:end,3) = outB_odd(:,1:2:end); outRGB(:,2:2:end,3) = outB_even(:,2:2:end);

% Type casting
outRGB = outRGB./255;
outIm = uint8(255*outRGB);

end

% Colorspace transfer function from RGB to XYZ
%   Caution! : You may need to change the mRGB2XYZ matrix
%   depending on the SIM2 HDR display characteristics.
function out = rgb2xyz_sim2(inp)
    % For 709
	mRGB2XYZ = [0.334785, 0.340199, 0.172095; 
				0.191003, 0.717471, 0.091526; 
				0.004775, 0.032711, 0.380929]';
            
    size1 = size(inp, 1);
    size2 = size(inp, 2);
    out = reshape(reshape(inp, [], 3)*mRGB2XYZ, size1, size2, 3);
end

% Colorspace transfer function from XYZ to Yuv
function out = xyz2yuv_m(inp)

	X = inp(:,:,1);
	Y = inp(:,:,2);
	Z = inp(:,:,3);

	% Take Y
	out(:,:,1) = Y;

	% For u channel,
	out(:,:,2) = (((1626.6875 * X) ./ (X + 15.0 * Y + 3.0 * Z)) + 0.546875) * 4.0;

	% For v channel,
	out(:,:,3) = (((3660.046875 * Y) ./ (X + 15.0 * Y + 3.0 * Z)) + 0.546875) * 4.0;
end

% Clamp function
function out = imClamp(inp, minVal, maxVal)
	temp = inp;
	temp(temp<minVal) = minVal;
	temp(temp>maxVal) = maxVal;
	out = temp;
end