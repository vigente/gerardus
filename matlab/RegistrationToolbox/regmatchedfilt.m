function [tElx, cmax, imm] = regmatchedfilt(imf, imm, alpha)
% REGMATCHEDFILT  Matched filter registration for translation and rotation.
%
% REGMATCHEDFILT uses a matched filter approach to find a global optimum
% for the rigid (translation and rotation) registration of two images.
%
% Assuming white noise, if we have signal s(t) with duration T, its matched
% filter is h(t)=s(T-t). If we convolve s(t)*h(t), the convolution has a
% peak at t=T, when both signals best overlap. We apply this idea to image
% registration, with the fixed image as s(t) and the moving image as h(t).
% This function uses FFTs to implement the convolution.
%
% An alternative is phase correlation registration (imregcorr). Phase
% correlation is similar to a matched filter, but normalizes the magnitude
% of the transforms in Fourier space. This may work better or worse than
% the matched filter, depending on the characteristics of the images.
%
% Matched filters and phase correlation in principle use brute force for
% translations. Matlab's function imregcorr uses Keller (2005) to find the
% angle and scale faster with a domain transform. However, when applied to
% our histology images we find that the angle estimation fails
% substantially. That's why we have implemented REGMATCHEDFILT with an
% angle brute-force approach. It's slower, but it produces good matches.
%
% TELX = REGMATCHEDFILT(IMF, IMM)
%
%   IMF, IMM are two grayscale images, not necessarily of the same size.
%   IMF is the fixed image, and IMM is the moving image. The images can be
%   provided as arrays or, if physical units are required, as SCIMAT
%   structs (see "help scimat"). When arrays are used, then we assume pixel
%   size = [1 1], offset = [0 0].
%
%   TELX is the transform in Elastix format that registers IMM onto IMF.
%
% TELX = REGMATCHEDFILT(..., ALPHA)
%
%   ALPHA is a vector of rotation angles. The centre of rotation is
%   computed internally as the centre of mass of IMM. The translation
%   optimum is found for each ALPHA(i) value. The ALPHA(i) with the best
%   match is chosen as the optimal angle.
%
% [TELX, CMAX, IMMREG] = REGMATCHEDFILT(...)
%
%   CMAX is a scalar with the value of the convolution peak.
%
%   IMMREG is the image resulting from applying the translation transform
%   to IMM. It has the same size and pixel size as IMF.
%
%
% Y Keller, "Pseduo-polar based estimation of large translations rotations
% and scalings in images", Application of Computer Vision, 2005.
% WACV/MOTIONS 2005 Volume 1. 
%
% See also conv_fft2, imregcorr.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2015 University of Oxford
% Version: 0.3.0
% 
% University of Oxford means the Chancellor, Masters and Scholars of
% the University of Oxford, having an administrative office at
% Wellington Square, Oxford OX1 2JD, UK. 
%
% This file is part of Gerardus.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. The offer of this
% program under the terms of the License is subject to the License
% being interpreted in accordance with English Law and subject to any
% action against the University of Oxford being under the jurisdiction
% of the English Courts.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
narginchk(2, 3);
nargoutchk(0, 3);

% convert images to scimat structs, so that we can use pixel size and
% offset
if (~isstruct(imf))
    imf = scimat_im2scimat(imf);
end
if (~isstruct(imm))
    imm = scimat_im2scimat(imm);
end

if (size(imm.data, 5) == 3)
    imm.data = rgb2gray(squeeze(imm.data));
elseif (size(imm.data, 5) ~= 1)
    error('IMM must be a grayscale or RGB image')
end
if (size(imf.data, 5) == 3)
    imf.data = rgb2gray(squeeze(imf.data));
elseif (size(imf.data, 5) ~= 1)
    error('IMF must be a grayscale or RGB image')
end

% defaults
if (nargin < 3)
    alpha = 0;
end

% coordinates of the center of mass of the image. We are going to use this
% as the centre of rotation. Units are (row, col) pixels, not physical
% coordinates
cm = imcmass(single(imm.data));

% if we convolve imf with itself, then the convolution peak would be in the
% center of the image. An (N,M)-image convolved with itself produces a
% (2*N-1,2*M-1) convolution. The center of that convolution is at (N,M)
conv_centre = size(imf.data);

% template for the elastix transform for the rotation only
tElx = elastix_transform(imf.data, imm.data, cm([2 1]));

% loop angle values
cmax = -Inf; % convolution maximum for any angle so far
for I = 1:length(alpha)
    
    % update the angle
    tElx.TransformParameters(1) = alpha(I);
    
    % apply rotation to the moving image
    imm2 = transformix(tElx, imm.data);
    
    % convolution of both images
    aux = conv_fft2(single(imm2), single(imf.data(end:-1:1, end:-1:1)), 'full');
    
    % maximum of current convolution
    [auxmax, idx] = max(aux(:));
    
    % if we have improved on the convolution maximum, we keep the current
    % angle as best so far
    if (auxmax > cmax)
        
        % update convolution peak and best angles so far
        cmax = auxmax;
        alphaopt = alpha(I);
        
        % coordinates of convolution maximum (in row, col format)
        [ropt, copt] = ind2sub(size(aux), idx);
        
    end

end

% set rotation to the best match
tElx.TransformParameters(1) = alphaopt;

% elastix transform for the translation only
tElx2 = tElx;
tElx2.TransformParameters(1) = 0;
tElx2.TransformParameters(2:3) = [copt, ropt] - conv_centre([2 1]);

% combine the rotation and translation
tElx = elastix_compose_afftransf(tElx2, tElx);

% convert pixel to physical units
tElx.TransformParameters(2:3) = tElx.TransformParameters(2:3) ...
    .* [imf.axis([2 1]).spacing];
tElx.Spacing = [imf.axis([2 1]).spacing];
tElx.CenterOfRotationPoint = tElx.CenterOfRotationPoint ...
    .* [imf.axis([2 1]).spacing];

% if requested by the user, transform the moving image according to the
% registration solution
if (nargout > 2)
    imm = transformix(tElx, imm);
end

end

%% nested functions

%% elastix_transform:
%  imf, imm: fixed and moving images
%  rotc: centre of translation in (x, y) format
function tElx = elastix_transform(imf, imm, rotc)

% set up the translation transform to align the histology to the
% blockface
tElx.Transform = 'EulerTransform';
tElx.NumberOfParameters = 3;
tElx.TransformParameters = [0 0 0];
tElx.InitialTransformParametersFileName = 'NoInitialTransform';
tElx.HowToCombineTransforms = 'Compose';
tElx.FixedImageDimension = 2;
tElx.MovingImageDimension = 2;
tElx.FixedInternalImagePixelType = 'float';
tElx.MovingInternalImagePixelType = 'float';
tElx.Size([2 1]) = size(imf);
tElx.Index = [0 0];
tElx.Spacing = [1 1];
tElx.Origin = [0 0];
tElx.Direction = [1 0 0 1];
tElx.UseDirectionCosines = 'true';
tElx.CenterOfRotationPoint = rotc;
tElx.ResampleInterpolator = 'FinalLinearInterpolator';
tElx.Resampler = 'DefaultResampler';
tElx.DefaultPixelValue = median(imm(:));
tElx.ResultImageFormat = 'mha';
tElx.ResultImagePixelType = 'unsigned char';
tElx.CompressResultImage = 'true';

end
