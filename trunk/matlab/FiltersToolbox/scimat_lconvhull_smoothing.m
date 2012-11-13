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
% Version: 0.3.3
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

% if the image has all voxels == 1, then the smoothed local convex hull is
% the image itself and we don't need to waste time doing other computations
%
% likewise, if the image has all voxels == 0, there's nothing to smooth
if ((nnz(scimat.data) == 0) || (numel(scimat.data) == nnz(scimat.data)))
    return
end

% coordinates of segmented voxels
[r, c, s] = ind2sub(size(scimat.data), find(scimat.data));
x = scinrrd_index2world([r c s], scimat.axis);
clear r c s

% if we don't have at least 4 unique points, the delaunay triangulation
% will give an error
if (size(x, 1) < 4)
    return
end

% compute alpha shape
try
    [~, s] = alphavol(x, alpha);
catch err
    % "Error computing the Delaunay triangulation. The points may be
    % coplanar or collinear"
    %
    % Just exit the function leaving the input segmentation untouched
    if (strcmp(err.identifier,'MATLAB:delaunay:EmptyDelaunay3DErrId'))
        return
   else
       % display any other errors as usual
      rethrow(err);
    end
end

% keep only the mesh surface
tri = s.bnd;
clear s

% very small segmentations can produce degenerated surfaces. In that case,
% there's nothing to smooth and we can just exit without modifying the
% input segmentation
if (size(tri, 1) < 4)
    return
end

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

% % initialise the output
% scimat.data = false(size(scimat.data));

% tight box around the surface
aux = x(tri(:), :);
box = [min(aux); max(aux)];

% world-coordinates to index values. We only need to search within this
% box, as voxels outside the box are definitely outside the surface too
boxidx = round(scinrrd_world2index(box, scimat.axis));

% all voxels in the original segmentation will be in the local convex hull.
% So we only need to check background voxels within the box

% linear indices of background voxels within the box
idx = find(~scimat.data(...
    boxidx(1, 1):boxidx(2, 1), ...
    boxidx(1, 2):boxidx(2, 2), ...
    boxidx(1, 3):boxidx(2, 3) ...
    ));

% box linear indices => box r, c, s
[r, c, s] = ind2sub(diff(boxidx) + 1, idx);

% convert box indices to whole image indices
r = r + boxidx(1, 1) - 1;
c = c + boxidx(1, 2) - 1;
s = s + boxidx(1, 3) - 1;

% r, c, s => world coordinates
xyz = scinrrd_index2world([r c s], scimat.axis);

% linear indices
idx = sub2ind(size(scimat.data), r, c, s);

% convert the triangulation surface to a binary mask
scimat.data(idx) = cgal_insurftri(tri, x, xyz, rand(3, 3));
