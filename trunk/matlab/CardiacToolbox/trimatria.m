function tis = trimatria(tis, surf)
% TRIMATRIA  Remove the atria from a whole cardiac tissue segmentation
%
% TIS2 = trimatria(TIS, SURF)
%
%   TIS is a 3D segmentation of a whole heart (voxels==0 correspond to
%   background, and voxels==1 correspond to tissue).
%
%   SURF is a 3D segmentation of the atrio-ventricular surface that
%   separates the ventricles from the atria.
%
%   TIS2 is TIS after having removed the atria from the segmentation, as
%   well as smaller non-connected components.

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
narginchk(2, 2);
nargoutchk(0, 1);

% dilate AVS
surf = itk_imfilter('bwdilate', surf, 2, 1);

% remove AVS from tissue segmentation
tis(surf ~= 0) = 0;

% compute connected components
cc = bwconncomp(tis);

% number of voxels in each component
n = cellfun(@(x) length(x), cc.PixelIdxList);

% component with the largets number of voxels
[~, idx] = max(n);

% keep only the largest component (ventricles)
tis = zeros(size(tis), 'uint8');
tis(cc.PixelIdxList{idx}) = 1;
