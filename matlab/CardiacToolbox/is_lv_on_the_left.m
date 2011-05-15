function flag = is_lv_on_the_left(nrrd, m, rot)
% IS_LV_ON_THE_LEFT  Check whether the left ventricle is on the left hand
% side of the image
%
% FLAG = IS_LV_ON_THE_LEFT(NRRD, M, ROT)
%
%   NRRD is a struct with the whole tissue segmentation.
%
%   M, ROT are a 3-vector and a (3, 3)-matrix that represent the centroid
%   and main axes of the heart, respectively.
%
%   FLAG is a boolean. It is true of the left ventricle is on the left hand
%   side of the image with respect to the coordinate system defined by ROT.
%
%
%   Note on SCI NRRD: Software applications developed at the University of
%   Utah Scientific Computing and Imaging (SCI) Institute, e.g. Seg3D,
%   internally use NRRD volumes to store medical data.
%
%   When label volumes (segmentation masks) are saved to a Matlab file
%   (.mat), they use a struct called "scirunnrrd" to store all the NRRD
%   information:
%
%   >>  scirunnrrd
%
%   scirunnrrd = 
%
%          data: [4-D uint8]
%          axis: [4x1 struct]
%      property: []
%

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010 University of Oxford
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
error( nargchk( 3, 3, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

% cut a whole tissue plane orthogonal to the vertical axis
im = scinrrd_intersect_plane(nrrd, m, rot(:, 3));

% the computed "X axis" of the orientation can be pointing more or less in
% the direction of of X or -X. We compute the angle
alpha = acosd(dot(rot(:, 1), [1 0 0]'));

% now we rotate the intersected plane to get an idea of whether the LV is
% going to be on the left or right hand side of the rotated volume
im = imrotate(im, alpha);

% label connected components
[lab, nlabs] = bwlabel(~im);

% get number of voxels in each label
nvox = zeros( nlabs, 1 );
for L = 1:nlabs
    nvox(L) = sum( lab(:) == L );
end

% sort labels by decreasing size
[nvox, idx] = sort(nvox, 1, 'descend');

% the largest component is the background; the 2nd and 3rd largest should
% be the RV and LV, in either order

% compute the eccentricity of each label (the opposite of the roundness)
stat2 = regionprops(double(lab == idx(2)), 'Eccentricity');
stat3 = regionprops(double(lab == idx(3)), 'Eccentricity');

% we vote whether the 2nd largest component is actually the RV: it is more
% eccentric
if (stat2.Eccentricity > stat3.Eccentricity)
    idxrv = idx(2);
    idxlv = idx(3);
else
    idxrv = idx(3);
    idxlv = idx(2);
end

% get the centroids of both components
crv = regionprops(double(lab == idxrv), 'Centroid');
clv = regionprops(double(lab == idxlv), 'Centroid');

% is the LV placed on the left hand side?
flag = clv.Centroid(1) < crv.Centroid(1);
