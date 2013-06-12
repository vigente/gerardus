function bw = bwsplitwithsurf(bw)
% BWSPLITWITHSURF  Split a volume into connected components separated by a
% surface
%
% BW2 = bwsplitwithsurf(BW)
%
%   BW is a 3-d segmentation that represents a curved surface. Background
%   voxels are equal to 0, and surface voxels are equal to 1.
%
%   BW2 is a segmentation of the same size as BW, where now all voxels to
%   one side of the surface are 0 and all voxels to the other side are 1.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
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
narginchk(1, 1);
nargoutchk(0, 1);

% to simplify the code, we tailor it to 3-d segmentations only
if (ndims(bw) ~= 3)
    error('BW must be a 3D array')
end

% dilate surface by 5 voxels
im = itk_imfilter('bwdilate', bw, 5, 1);

% invert the segmentation so that the dilated surface will be
% considered background
im = 1 - im;

% get a tight box around the surface, and expand with a 1 voxel margin to
% avoid disconnecting 
idx = bwbox(bw, 1);

% extract the box. By working on a smaller image, we can speed up
% computations
im2 = im(...
    idx(1, 1):idx(1, 2), ...
    idx(2, 1):idx(2, 2), ...
    idx(3, 1):idx(3, 2));

% find connected components
cc2 = bwconncomp(im2);

% assign a different label to each connected component's voxels
for I = 1:cc2.NumObjects
    im2(cc2.PixelIdxList{I}) = I;
end

% label the empty dilated surface with equal contributions from each
% component (30 min for 650x500x1650)
TODO = zeros(1, class(im2)) + I + 1;
im2(im2 == 0) = TODO;
im2 = bwregiongrow(im2, TODO);

%% extend the labels from the smaller image to the whole image

% find connected components in whole image
cc = bwconncomp(im);

if (cc.NumObjects ~= cc2.NumObjects)
    error('cc and cc2 should have the same number of components')
end

% loop components in the images
for I = 1:cc.NumObjects
    % index of one of the voxels in smaller image component
    ind2 = cc2.PixelIdxList{I}(1);
    
    % translate to whole image coordinates
    [r, c, s] = ind2sub(size(im2), ind2);
    r = r - 1 + idx(1, 1);
    c = c - 1 + idx(2, 1);
    s = s - 1 + idx(3, 1);
    ind = sub2ind(size(im), r, c, s);
    
    % find to which whole image component the current testing voxel belongs
    % to
    ind = cellfun(@(x) any(x==ind), cc.PixelIdxList);
    
    % labell pixels in the whole image component accordingly
    bw(cc.PixelIdxList{ind}) = I;
end

% put the smaller image into the larger image
bw(...
    idx(1, 1):idx(1, 2), ...
    idx(2, 1):idx(2, 2), ...
    idx(3, 1):idx(3, 2)) = im2;

% substract 1 from voxel labels so that one of the components defined by
% the surface can be considered background
bw = bw - 1;
