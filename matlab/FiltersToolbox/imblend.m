function [im, a0] = imblend(im1, im2)
% IMBLEND  Blend two image stacks to increase the apparent dynamic range
% and improve the signal to noise ratio
%
% A common problem in photography or microscopy is trying to image a target
% with strong contrast. If the settings are optimised for the highlights,
% information is lost as the shadow areas get underexposed. On the other
% hand, if the settings are optimised for the shadows, information is lost
% when the highlights get overblown.
%
% The most widely used method to deal with this problem in digital cameras
% is auto-bracketing. The camera takes several images at different light
% settings, and then blends them in an intelligent way to increase the
% apparent dynamic range of the sensor.
%
% This function implements a blending function for two images to increase
% the apparent dynamic range by removing under- and over-exposed voxels,
% and also improves the signal to noise ratio of the final image.
%
% IM = IMBLEND(IM1, IM2)
%
%   IM1, IM2 are 3D arrays with two images of the same scene. The images
%   may have been taken at different light settings, or with the same ones.
%
%   IM is the output 3D array that results from blending the input images.
%
% [IM, A0] = IMBLEND(IM1, IM2)
%
%   A0 is the scaling factor such that IM*A0 is closest in intensity values
%   to IM2, once saturated and underexposed voxels are removed.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.1.0
% $Rev$
% $Date$
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
error(nargchk(2, 2, nargin, 'struct'));
error(nargoutchk(0, 2, nargout, 'struct'));

% check image sizes
if any(size(im1) ~= size(im2))
    error('Input images must be the same size')
end

% saturation value for the input image types
IMAX1 = double(zeros(1, class(im1)) + Inf);
IMAX2 = double(zeros(1, class(im2)) + Inf);

% find voxels that are either saturated or underexposed
sat1 = im1 == IMAX1;
sat2 = im2 == IMAX2;
und1 = im1 == 0;
und2 = im2 == 0;
bad = sat1 | sat2 | und1 | und2;

% create vector with the voxels at two different levels of gain
v1 = double(im1(~bad));
v2 = double(im2(~bad));

% compute the median value in the 2nd image for each value in the first
% image
gmedian = zeros(1, IMAX1);
for I = 1:length(gmedian)
    gmedian(I) = median(v2(v1 == I));
end

% compute linearity measures as we increase the intensity
r = nan(1, IMAX1);
a = nan(1, IMAX1);
for I = 10:length(r)
    % correlation coefficient
    aux = corrcoef(1:I, gmedian(1:I));
    r(I) = aux(1, 2);
    
    % slope
    p = polyfit(1:I, gmedian(1:I), 1);
    a(I) = p(1);
end

% find the intensity with the maximum value for the correlation coefficient
[~, idx] = max(r);

% use that as the best estimate of the gain between the two images
a0 = a(idx);

% find maximum intensity we are going to have in the output image
IMAX = max(IMAX1 * a0, IMAX2);

% data type we need for the output image so that we don't clip any intesity
% values
if IMAX <= intmax('uint8')
    dataType = 'uint8';
elseif IMAX <= intmax('uint16')
    dataType = 'uint16';
elseif IMAX <= intmax('uint32')
    dataType = 'uint32';
else
    dataType = 'double';
end

% allocate memory for output image
im = zeros(size(im1), dataType);

% find voxels that are fine in each image
ok1 = ~und1 & ~sat1;
ok2 = ~und2 & ~sat2;

% voxels that are wrong in 1 image only will get their final value from the
% other image
idx = ok1 & ~ok2;
im(idx) = double(im1(idx)) * a0;

idx = ~ok1 & ok2;
im(idx) = im2(idx);

% voxels that are fine in both images, or bad in both images, will get a
% final value that is the average
idx = (ok1 & ok2) | (~ok1 & ~ok2);
im(idx) = (double(im1(idx)) * a0 + double(im2(idx))) / 2;
