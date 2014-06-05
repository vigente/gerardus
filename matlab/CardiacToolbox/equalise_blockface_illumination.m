function [imeq, illu] = equalise_blockface_illumination(im, polymask, ellipmask, ratio, thr, radheart, radpoly)
% EQUALISE_BLOCKFACE_ILLUMINATION  Correct illumination inhomogeneities in
% blockface photography.
%
% This function dilates a polygonal mask of the wax block so that it
% overlaps a bit with the background. Then, it obtains a rough segmentation
% of the heart using a threshold, and removes everything outside an
% elliptical mask. The wax block pixels minus the heart segmentation pixels
% are used to estimate the slow varying illumination of the block, and
% extrapolate to all pixels on the block. Finally, the image is corrected
% dividing by the estimated illumination.
%
% [IMEQ, ILLU] = equalise_blockface_illumination(IM, POLYMASK, ELLIPMASK)
%
%   IM is a grayscale or RGB blockface photograph of a wax embedded heart.
%
%   POLYMASK is a binary mask with the same size as IM. In POLYMASK, pixels
%   with value 1 belong to the wax block, and pixels with value 0 belong to
%   the background. This mask can be created drawing a polygonal line in
%   Matlab, and then converting it to a binary mask. This mask should be
%   drawn close enough to the edges so that when the algorithm dilates it,
%   it overlaps a bit with the background outside the wax block.
%
%     % display blockface image on the screen
%     him = imagesc(im);
%     % start polygonal line tool
%     h = impoly;
%     % now you can draw the mask on the screen with the mouse
%     polymask = createMask(h, him);
%     save('polymask.mat', 'polymask')
%
%   ELLIPMASK is a binary mask too. In ELLIPMASK, pixels equal to 1
%   correspond to the inside of an ellipse that contains the heart. This
%   mask can be generated similarly to POLYMASK, replacing impoly by
%
%     % start ellipse drawing tool
%     h = imellipse;
%
%   IMEQ is a grayscale image with the same size as IM, with the
%   illumination corrected.
%
%   ILLU is an image with the same size as IM, containing the estimate of
%   the illumination within POLYMASK.
%
% ... = equalise_blockface_illumination(..., RATIO, THR, RADHEART, RADPOLY)
%
%   RATIO is a scalar with the scale factor applied to reduce the image to
%   speed up the estimation of the illumination field. By default,
%   RATIO=1/16.
%
%   THR is a scalar with the intensity threshold for the rough heart
%   segmentation. By default, THR=45.
%
%   RADHEART is a scalar with the radius of the dilation applied to the
%   rought heart segmentation. By default, RADHEART=10.
%
%   RADPOLY is a scalar with the radius of the dilation applied to
%   POLYMASK.
%
%   Example:
%
%     % load blockface image
%     im = imread('test/R_55_0001.bmp');
%
%     % load polygonal and elliptical masks
%     load('test/polymask.mat')
%     load('test/ellipmask.mat')
%
%     % set parameters
%     ratio = 1/16;
%     thr = 45;
%     radheart = 10;
%     radpoly = 50;
%
%     % correct illumination
%     [imeq, illu] = equalise_blockface_illumination(im, polymask, ...
%     ellipmask, ratio, thr, radheart, radpoly);
%
%     % plot image and results
%     subplot(3, 1, 1)
%     imagesc(im)
%     subplot(3, 1, 2)
%     imagesc(imeq)
%     subplot(3, 1, 3)
%     imagesc(illu)

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

DEBUG = false;

% check arguments
narginchk(3, 7);
nargoutchk(0, 2);

% defaults
if ((nargin < 4) || isempty(ratio))
    ratio = 1/16;
end
if ((nargin < 5) || isempty(thr))
    thr = 45;
end
if ((nargin < 6) || isempty(radheart))
    radheart = 10;
end
if ((nargin < 7) || isempty(radpoly))
    radpoly = 50;
end

% convert to grayscale
if (size(im, 3) == 3)
    im = rgb2gray(im);
end

% plot image
if (DEBUG)
    subplot(2, 1, 1)
    hold off
    him = imagesc(im);
end

% threshold mask of the heart
heartmask = (im < thr) & ellipmask;
se = strel('disk', radheart);
heartmask = imdilate(heartmask, se);

% plot heart mask
if (DEBUG)
    subplot(2, 1, 2)
    hold off
    imagesc(heartmask)
end

% dilate the polynomial mask so that it includes a bit of the background
se = strel('disk', radpoly);
polymask = imdilate(polymask, se);

% combine both masks
mask = xor(polymask, heartmask);
 
% plot mask
if (DEBUG)
    subplot(2, 1, 2)
    hold off
    imagesc(mask)
end

% reduce image size, because the illumination is a rather slow varying
% field
imsmall = imresize(im, ratio, 'bicubic');
masksmall = imresize(mask, ratio, 'nearest');
polymasksmall = imresize(polymask, ratio, 'nearest');

% plot downsized image
if (DEBUG)
    subplot(2, 1, 2)
    hold off
    imagesc(imsmall)
end

% grid to have coordinates for the pixels
[ysmall, xsmall] = ndgrid(1:size(imsmall, 1), 1:size(imsmall, 2));

% fit smoothing surface to illumination
f = fit([xsmall(masksmall) ysmall(masksmall)], ...
    double(imsmall(masksmall)), 'poly55', 'Normalize', 'on');
illusmall = feval(f, xsmall(:), ysmall(:));
illusmall = reshape(illusmall, size(imsmall));
illusmall(~polymasksmall) = Inf;

if (DEBUG)
    % plot estimated illumination
    subplot(2, 1, 2)
    hold off
    imagesc(illusmall)
end


% resize illumination to the full resolution
illu = imresize(illusmall, size(im));
illu(isnan(illu)) = Inf;

% equalise illumination
imeq = double(im) ./ illu;

% plot illumination profiles
if (DEBUG)
    subplot(2, 1, 2)
    hold off
    imagesc(imeq)
end
