function scimat = scimat_closed_surf_to_bw(tri, x, scimat)
% SCIMAT_CLOSED_SURF_TO_BW  Segment the inside of the triangulation of a
% closed surface
%
% SCIMAT2 = scimat_closed_surf_to_bw(TRI, X, SCIMAT)
%
%   TRI is a 3-column matrix with a triangulation of a close surface. Each
%   element in TRI is an index to a row in X. Each row represents the three
%   vertices of a triangle on the surface.
%
%   X is a 3-column matrix with the Euclidean coordinates of the
%   triangulation vertices.
%
%   SCIMAT is a struct where the result will be saved to.
%
%   SCIMAT2 is the output, where the voxels within the closed surface have
%   been set to 1. Note that these voxels are _added_ to whatever previous
%   segmentation was present in SCIMAT.
%
% See also: delaunay_sphere.

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
narginchk(3, 3);
nargoutchk(0, 1);

if (size(tri, 2) ~= 3)
    error('TRI must be a 3-column matrix')
end
if (size(x, 2) ~= 3)
    error('X must be a 3-column matrix')
end

% find a tight box around the segmentation
% we only need to query the voxels within the box, whether they are inside
% or outside the surface. The rest, we know that they are outside
boxmin = min(x, [], 1);
boxmax = max(x, [], 1);

% real world coordinates to indices
idxmin = floor(scinrrd_world2index(boxmin, scimat.axis));
idxmax = ceil(scinrrd_world2index(boxmax, scimat.axis));

% make sure that we don't go outside the volume
idxmin = max([idxmin; 1 1 1], [], 1);
idxmax = max([idxmax; 1 1 1], [], 1);
idxmin = min([idxmin; size(scimat.data)], [], 1);
idxmax = min([idxmax; size(scimat.data)], [], 1);

% convert back to real world coordinates
boxmin = scinrrd_index2world(idxmin, scimat.axis);
boxmax = scinrrd_index2world(idxmax, scimat.axis);

% create the cell array with the vectors that the rectangular grid of the
% voxels that we are going to query
ci = {linspace(boxmin(1), boxmax(1), idxmax(2)-idxmin(2)+1), ...
    linspace(boxmin(2), boxmax(2), idxmax(1)-idxmin(1)+1), ...
    linspace(boxmin(3), boxmax(3), idxmax(3)-idxmin(3)+1)};

% query voxels to find which ones are inside and which ones outside the
% segmentation
scimat.data(idxmin(1):idxmax(1), ...
    idxmin(2):idxmax(2), ...
    idxmin(3):idxmax(3)) = cgal_insurftri(tri, x, ci, rand(3, 3));
