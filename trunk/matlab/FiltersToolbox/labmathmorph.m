function im = labmathmorph(type, im, param)
% LABMATHMORPH  Mathematical morphology operators on a labelled
% segmentation, one label at a time
%
% This function runs mathematical morphology operations (dilation, erosion,
% opening, closing) on labelled segmentations. The operator is applied to
% each label independently. Care has been taken speed up operations for
% large images.
%
% IM2 = LABMATHMORPH(TYPE, IM, PARAM)
%
%   TYPE is a string with the operator name:
%
%     * 'dilate': Binary dilation
%
%     * 'erode': Binary erosion
%
%     * 'open': Erosion followed by dilation
%
%     * 'close': Dilation followed by erosion (i.e. smoothing)
%
%   IM is a 2D or 3D volume with the labelled segmentation. All voxels with
%   the same value will be considered to belong to the same label.
%
%   PARAM is a value or vector with the parameter(s) values for the
%   operator. See below for details.
%
% IM2 = LABMATHMORPH('dilate', IM, NDIL)
%
% IM2 = LABMATHMORPH('erode', IM, NERO)
%
% IM2 = LABMATHMORPH('close', IM, [NDIL NERO])
%
% IM2 = LABMATHMORPH('open', IM, [NERO NDIL])
%
% where NDIL, NERO are the radii to dilate and erode, respectively, in
% voxel units.

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
error(nargchk(3, 3, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% number of voxels we need to leave around the cropped branch. Note that we
% need to leave at least 1 voxel of empty space, otherwise the ITK
% morphological operators do weird things at the boundary
gap = 1;

% image size
sz = size(im);
if (length(sz) == 2)
    sz = [sz 1];
end

% get parameter values
switch type
    case 'dilate'
        if (length(param) == 1)
            ndil = param;
            gap = ndil + 1;
        else
            error('Dilation operator expects 1 value in PARAM')
        end
    case 'erode'
        if (length(param) == 1)
            nero = param;
            gap = nero + 1;
        else
            error('Erosion operator expects 1 value in PARAM')
        end
    case 'close'
        if (length(param) == 2)
            ndil = param(1);
            nero = param(2);
            gap = ndil + 1;
        else
            error('Closing operator expects a 2-vector in PARAM')
        end
    case 'open'
        if (length(param) == 2)
            nero = param(1);
            ndil = param(2);
            gap = ndil + 1;
        else
            error('Opening operator expects a 2-vector in PARAM')
        end
    otherwise
        error('Operator type not implemented')
end

% get indices and labels of the segmented voxels
idxlab = find(im);
lab = nonzeros(im);

% sort the label values
[lab, idx] = sort(lab);
idxlab = idxlab(idx);

% find where each label begins. The last index is "fake", i.e. it doesn't
% correspond to any label, but it is used to know where the last label ends
idxlab0 = [0 ; find(diff(lab)) ; length(lab)] + 1;

% list of labels with at least 1 voxel in the segmentation
LAB = unique(lab);

% free some memory
clear lab

% loop every branch
for I = 1:length(LAB)

    % list of voxels in current branch. The reason why we are not doing a
    % simple br = find(nrrd.data == LAB(I)); is because for a large volume,
    % that's a comparatively very slow operation
    br = idxlab(idxlab0(I):idxlab0(I+1)-1);
    
    % crop the part of the segmentation that contains the branch
    [r, c, s] = ind2sub(sz, br);
    
    from = min([r c s], [], 1);
    to = max([r c s], [], 1);
    
    from = max([from - gap ; 1 1 1], [], 1);
    to = min([to + gap ; sz], [], 1);
    
    imcrop0 = im(from(1):to(1), from(2):to(2), from(3):to(3));
    
    % voxels of the current branch
    imcrop = (imcrop0 == LAB(I));
    
    % remove current branch from the rest of the voxels
    imcrop0(imcrop) = 0;
    
    % run operator on the branch
    switch type
        case 'dilate'
            imcrop = itk_imfilter('bwdilate', imcrop, ndil, true);

        case 'erode'
            imcrop = itk_imfilter('bwerode', imcrop, nero, true);
            
        case 'close'
            imcrop = itk_imfilter('bwdilate', imcrop, ndil, true);
            imcrop = itk_imfilter('bwerode', imcrop, nero, true);
            
        case 'open'
            imcrop = itk_imfilter('bwerode', imcrop, nero, true);
            imcrop = itk_imfilter('bwdilate', imcrop, ndil, true);
    end
    
    % remove voxels from other branches from the result
    imcrop(imcrop0 ~= 0) = 0;
    
    % get indices of the resulting voxels
    [r, c, s] = ind2sub(size(imcrop), find(imcrop));
    
    % convert indices from crop to whole segmentation
    idx = sub2ind(size(im), r + from(1) - 1, c + from(2) - 1, ...
        s + from(3) - 1);
    
    % remove original branch from whole segmentation
    im(br) = 0;
    
    % add filtered branch to whole segmentation
    im(idx) = LAB(I);
    
end
