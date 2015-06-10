function [tElx, cmax, imm] = regmatchedfilt(imf, imm)
% REGMATCHEDFILT  Matched filter registration.
%
% REGMATCHEDFILT uses a matched filter approach to find a global optimum
% for the translation registration of two images.
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
% TELX = REGMATCHEDFILT(IMF, IMM)
%
%   IMF, IMM are two grayscale images, not necessarily of the same size.
%   IMF is the fixed image, and IMM is the moving image.
%
%   TELX is the transform in Elastix format that register IMM onto IMF.
%
% [TELX, CMAX, IMMREG] = REGMATCHEDFILT(...)
%
%   CMAX is a scalar with the value of the convolution peak.
%
%   IMMREG is the image resulting from applying the translation transform
%   to IMM.
%
% See also conv_fft2, imregcorr.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2015 University of Oxford
% Version: 0.1.0
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
narginchk(2, 2);
nargoutchk(0, 3);

if (size(imm, 3) > 1)
    error('IMM must be a grayscale image')
end
if (size(imf, 3) > 1)
    error('IMF must be a grayscale image')
end

% convolution of both images
imc = conv_fft2(double(imm), double(imf(end:-1:1, end:-1:1)), 'full');

% coordinates of convolution maximum (in row, col format)
[rmax, cmax] = find(imc == max(imc(:)));

% if we convolve imf with itself, then the convolution peak would be in the
% center of the image. An (N,M)-image convolved with itself produces a
% (2*N-1,2*M-1) convolution. The center of that convolution is at (N,M)
centre = size(imf);

% set up the translation transform to align the histology to the
% blockface
tElx.Transform = 'TranslationTransform';
tElx.NumberOfParameters = 2;
tElx.TransformParameters([2 1]) = [rmax, cmax] - centre;
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
tElx.CenterOfRotationPoint([2 1]) = centre;
tElx.ResampleInterpolator = 'FinalLinearInterpolator';
tElx.Resampler = 'DefaultResampler';
tElx.DefaultPixelValue = median(imm(:));
tElx.ResultImageFormat = 'png';
tElx.ResultImagePixelType = 'unsigned char';
tElx.CompressResultImage = 'true';

% if requested by the user, transform the moving image according to the
% registration solution
if (nargout > 2)
    imm = transformix(tElx, imm);
end
