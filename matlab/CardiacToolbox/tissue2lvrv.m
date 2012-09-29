function [im, LV, RV, BKG, TODO] = tissue2lvrv(im, surf, RAD, blksz, border, numworkers)
% TISSUE2LVRV  Extract LV and RV cavities from tissue and atrioventricular
% surface segmentations
%
% This function takes at the input a segmentation of the cardiac muscle and
% a segmentation of the atrioventricular surface (or any other upper limit
% for the ventricles), and returns 
%
% [LVRV, LV, RV, BKG] = tissue2lvrv(TIS, SURF)
%
%   TIS and SURF are segmentation volumes of the same size. TIS contains a
%   segmentation of the cardiac tissue, and SURF a segmentation of the
%   atrioventricular surface (AVS). This surface is a natural way to set a
%   boundary to the ventricular cavities.
%
%   LVRV is a segmentation volume of the same size as TIS and SURF. It
%   contains the segmentations of the left ventricle (LV) and right
%   ventricle (RV) cavities, and the background.
%
%   LV, RV, BKG are scalars with the labels of each output segmentation:
%   LV = left ventricle cavity, RV = right ventricle cavity, BKG =
%   background.
%
% [..., TODO] = tissue2lvrv(..., RAD, BLKSZ, BORDER, NUMWORKERS)
%
%   RAD is a scalar. Quite often, the AVS is too thin to separate the
%   ventricles from the background, and segmentation errors may connect the
%   ventricles between them or to the background. In order to segment the
%   cavities, all the background is eroded by radius RAD, and the AVS is
%   dilated by radius RAD too. By default, RAD = 2.
%
%   BLKSZ, BORDER, NUMWORKERS are parameters to run the last part of the
%   algorithm in parallel by blocks. BLKSZ is a 3-vector with the size in
%   voxels of each block the segmentation will be split into. BORDER is a
%   3-vector with the overlap in voxels between blocks. NUMWORKERS is a
%   scalar with the maximum number of cores the algorithm is allowed to
%   use. By default, BLKSZ=SIZE(TIS) and no parallel processing is used. By
%   default too, BORDER=[0 0 0], NUMWORKERS=1.
%
%   TODO is the label of voxels that were not labelled by the algorithm.
%   This can happen with parallel processing, as each block can split the
%   segmentation into non-connected components.
%
%   Note that parallel processing is possible by selecting BLKSZ smaller
%   than the volume size, even if NUMWORKERS=1.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012 University of Oxford
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
narginchk(2, 6);
nargoutchk(0, 4);

if any(size(im) ~= size(surf))
    error('IM and SURF must have the same size')
end

% defaults
if (nargin < 3) || isempty(RAD)
    RAD = 2;
end
if (nargin < 4) || isempty(blksz)
    blksz = size(im);
end
if (nargin < 5) || isempty(border)
    border = [0 0 0];
end
if (nargin < 6) || isempty(numworkers)
    numworkers = 1;
end

% invert mask to obtain background (which includes the cavities)
im = uint8(~im);

% dilate the surface so that it's water-tight
surf = itk_imfilter('bwdilate', surf, RAD, 1);
    
% remove atrio-ventricular surface from segmentation
im(surf ~= 0) = 0;
    
% duplicate background segmentation
im2 = im;
    
% erode cavity segmentation
im2 = itk_imfilter('bwerode', im2, RAD, 1);
    
% find connected components
cc = bwconncomp(im2);
    
% size of each connected component
numPixels = cellfun(@numel, cc.PixelIdxList);

% sort from largest to smallest components
[~, idx] = sort(numPixels, 'descend');
cc.PixelIdxList = cc.PixelIdxList(idx);
    
% labels for the region grow algorithm
BKG = 1;
LV = 2;
RV = 3;
TODO = 4;

% initialise all segmented voxels as TODO
im(im ~= 0) = TODO;

% top part of the background
idx = find(cellfun(@(x) any(x==1), cc.PixelIdxList));
im(cc.PixelIdxList{idx}) = BKG;
cc.PixelIdxList(idx) = [];
cc.NumObjects = cc.NumObjects - 1;

% bottom part of the background
idx = sub2ind(size(im), 1, 1, size(im, 3));
idx = find(cellfun(@(x) any(x==idx), cc.PixelIdxList));
im(cc.PixelIdxList{idx}) = BKG;
cc.PixelIdxList(idx) = [];
cc.NumObjects = cc.NumObjects - 1;
    
% the ventricles will be the next largest components
cc.PixelIdxList(3:end) = [];
cc.NumObjects = 2;
    
% compute "solidity", that is, how many voxels the components has divided
% by the number of voxels in the bounding box
stats = regionprops(cc, 'Area', 'BoundingBox');
solidity(1) = stats(1).Area / prod(stats(1).BoundingBox(4:6));
solidity(2) = stats(2).Area / prod(stats(2).BoundingBox(4:6));

% the LV is more solid than the RV
islv = solidity > solidity(2:-1:1);

% tag the LV and RV
im(cc.PixelIdxList{islv}) = LV;
im(cc.PixelIdxList{~islv}) = RV;
    
% grow each region unil they cover the whole segmentation mask
%
% we enable parallel computing if the user wants to process by blocks
% smaller than the whole image. Note that selecting only 1 worker itself
% doesn't mean we want to process the whole image in one go. Due to memory
% limitations, we may want to process the image in smaller chuncks, but at
% the same time have only 1 worker doing it
if (all(blksz == size(im)))
    im = bwregiongrow(im, uint8(TODO));
else
    fun = @(x) bwregiongrow(x, uint8(TODO));
    im = blockproc3(im, blksz, fun, border, numworkers);
end
