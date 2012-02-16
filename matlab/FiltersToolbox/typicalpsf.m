function [psf, n] = typicalpsf(im, sidesz, thr, minsz)
% TYPICALPSF  Estimate point spread function from microbeads image
%
% PSF = TYPICALPSF(IM, SIDESZ, THR, MINSZ)
%
%   PSF is an array with dimensions 2*SIDESZ+1 that contains the estimation
%   of the typical point spread function (PSF) obtained from IM.
%
%   IM is an array with an image of microbeads or any other spheres small
%   enough to be under the resolution power of the optical system, so that
%   their image can be considered roughly the PSF of the sytem plus some
%   noise.
%
%   THR is a scalar with the noise level. Any voxel with intensity < THR
%   will be made 0. By default, THR = QUANTILE(IM, .99). This assumes an
%   image with some bright points, but mostly background.
%
%   MINSZ is a scalar with the minimum size of a bead connected component
%   in voxels. Smaller components will be considered noise. By default,
%   MINSZ=4.
%
%   This function finds each potential bead in the image and counts the
%   number of voxels it contains. Only beads with size between the 1st and
%   3rd quantiles (25% to 75%) are kept. Those are aligned at their maximum
%   intensity, and then the median is used to find the typical value of
%   each voxel in the PSF.
%
%   It is necessary to provide an estimate of the bead image size. SIDESZ
%   contains the number of voxels at each side of the centre. That is, we
%   assume that the PSF is negligible beyond a box of size 2*SIDESZ+1.
%
% [PSF, N] = TYPICALPSF(...)
%
%   N is a scalar with the number of beads that were used to estimate the
%   PSF.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.2.0
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
error(nargchk(2, 4, nargin, 'struct'));
error(nargoutchk(0, 2, nargout, 'struct'));

% maximum intensity value for the image's data type
if isinteger(im)
    maxval = intmax(class(im));
elseif isnumeric(im)
    maxval = realmax(class(im));
else
    error('Image type is not numeric')
end

% defaults
if (nargin < 3 || isempty(thr))
    thr = quantile(im(:), .99);
else
    thr = quantile(im(:), thr);
end
if (nargin < 4 || isempty(minsz))
    minsz = 4;
end

% threshold to remove background noise
mask = im;
mask(im < thr) = 0;

% compute connected components
cc = bwconncomp(mask);
clear mask

% convert image to double
im = double(im);

% compute size of each component
stats = regionprops(cc, im, 'Area', 'MaxIntensity');

% remove connected components with small size, as they probably correspond to noise
idx = [stats.Area] <= minsz;
cc.PixelIdxList(idx) = [];
cc.NumObjects = nnz(~idx);
stats(idx) = [];

% get the 1st and 3rd quantiles of the size of the beads image
lo = quantile([stats.Area], .25);
hi = quantile([stats.Area], .75);

% list of components that are between the 1st and 3rd quantiles in size
idx = ([stats.Area] >= lo & [stats.Area] <= hi)';

% coordinates of the maximum of each object
m = zeros(cc.NumObjects, 3);
for I = 1:cc.NumObjects
    % whole image index of the voxel with maximum intensity
    aux = cc.PixelIdxList{I}(im(cc.PixelIdxList{I}) == stats(I).MaxIntensity);
    
    % idx = > r, c, s
    [m(I, 1), m(I, 2), m(I, 3)] = ind2sub(cc.ImageSize, aux(1));
end

% components too close to the edges will be ignored
bad = m(:, 1) <= sidesz(1) | m(:, 1) >= size(im, 1)-sidesz(1)+1 ...
    | m(:, 2) <= sidesz(2) | m(:, 2) >= size(im, 2)-sidesz(2)+1 ...
    | m(:, 3) <= sidesz(3) | m(:, 3) >= size(im, 3)-sidesz(3)+1;

% exclude components that have saturated intensity values
bad = bad | ([stats.MaxIntensity] >= maxval)';

idx = idx & ~bad;

% init matrix; each row of the matrix will contain the voxels in one box
psf = nan(nnz(idx), prod(2*sidesz+1));

% list of components to analyze
idx = find(idx)';

% number of beads that will be used to estimate the PSF
n = length(idx);

% loop every valid object
for I = 1:n
    
    % coordinates of centroid in an easier nomenclature
    r = m(idx(I), 1);
    c = m(idx(I), 2);
    s = m(idx(I), 3);
    
    % extract from the image a box centered on the centroid and with size
    % 2*sidesz+1
    box = im(r-sidesz(1):r+sidesz(1), ...
        c-sidesz(2):c+sidesz(2), ...
        s-sidesz(3):s+sidesz(3));
    
    % linearize box
    psf(I, :) = box(:)';
    
end

% get typical intensity values for each voxel of the PSF
psf = median(psf, 1);

% recover box shape
psf = reshape(psf, 2*sidesz+1);
