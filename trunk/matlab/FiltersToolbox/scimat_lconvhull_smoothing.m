function scimat = scimat_lconvhull_smoothing(scimat, alpha)
% SCIMAT_LCONVHULL_SMOOTHING  Smoothing of a binary image using a local
% convex hull
%
% SCIMAT2 = scimat_lconvhull_smoothing(SCIMAT, ALPHA)
%
%   SCIMAT is a SCI MAT volume with a binary image.
%
%   The perimeter voxels of the image are selected. The alpha-shape of
%   their coordinates produces a mesh triangulation with the local convex
%   hull.
%
%   ALPA is a scalar with the radius of the alpha shape, i.e. the size of
%   the convex hull neighbourhood. When the convex hull is computed, only
%   voxels within ALPHA distance can be connected. When ALPHA=Inf, the
%   convex hull is obtained.
%
%   The inside of the triangulation is converted to voxels using
%   cgal_insurftri().
%
%   ALPHA2 is the same SCI MAT volume with the local convex hull of SCIMAT.
%
% This function uses alphavol() by Jonas Lundgren.
%
% See also: cgal_insurftri.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012 University of Oxford
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
narginchk(2, 2);
nargoutchk(0, 1);

% coordinates of segmented voxels
[r, c, s] = ind2sub(size(scimat.data), find(scimat.data));
x = scinrrd_index2world([r c s], scimat.axis);
clear r c s

% compute alpha shape
[~, s] = alphavol(x, alpha);

% keep only the mesh surface
tri = s.bnd;
clear s

% % DEBUG: plot alpha shape's surface
% close all
% figure
% trisurf(tri, x(:, 1)*1e3, x(:, 2)*1e3, x(:, 3)*1e3, 'EdgeColor', 'none')
% axis xy equal
% view(-22, -4)
% camlight('headlight')
% lighting gouraud
% set(gca, 'FontSize', 18)
% xlabel('x (mm)')
% ylabel('y (mm)')
% zlabel('z (mm)')

% coordinates of first and last voxels in the image
box = scinrrd_index2world([1 1 1; size(scimat.data)], scimat.axis);

% vectors that define the sampling area of the image
xi = box(1, 1):scimat.axis(2).spacing:box(2, 1);
yi = box(1, 2):scimat.axis(1).spacing:box(2, 2);
zi = box(1, 3):scimat.axis(3).spacing:box(2, 3);

% convert the triangulation surface to a binary mask
scimat.data = cgal_insurftri(tri, x, {xi, yi, zi}, rand(5, 3));
