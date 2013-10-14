function [scimat, tri, x] = scimat_lconvhull_smoothing(scimat, rad, SAFE, dilrad)
% SCIMAT_LCONVHULL_SMOOTHING  Smoothing of a binary image using a local
% convex hull.
%
% SCIMAT2 = scimat_lconvhull_smoothing(SCIMAT, RAD)
%
%   SCIMAT is a SCI MAT volume with a binary image. This binary image is
%   assumed to not have holes. It can consist of several non-connected
%   components, but none of them can have any holes.
%
%   The perimeter voxels of the image are selected. The alpha-shape of
%   their coordinates produces a mesh triangulation with the local convex
%   hull.
%
%   Then, rays are traced from the centre of each voxel within the bounding
%   box of the mesh. The number of intersections with the mesh allows to
%   know whether the voxel is inside or outside of it. This part uses
%   cgal_insurftri(), which is very precise but can be a bit slow for very
%   large images.
%
%   RAD is a scalar with the radius of the alpha shape, i.e. the size of
%   the convex hull neighbourhood. When the convex hull is computed, only
%   voxels within RAD distance can be connected. When RAD=Inf, the
%   convex hull is obtained. Note that if you want to compare to other
%   alpha shape functions, e.g. cgal_alpha_shape3() or
%   cgal_fixed_alpha_shape3(), ALPHA=RAD^2.
%
%   SCIMAT2 is the local convex hull of SCIMAT.
%
% SCIMAT2 = scimat_lconvhull_smoothing(SCIMAT, RAD, false, DILRAD)
%
%   To allow for faster processing when the above method using
%   cgal_insurftri() is too slow, we allow an "unsafe" mode. This method
%   uses itk_tri_rasterization(), which is faster but less reliable.
%
%   The inside of the triangulation is converted to voxels using
%   itk_tri_rasterization(). As this function can misclassify quite a few
%   voxels that should be labelled as 1, we run it three times permuting
%   the coordinates of the mesh and the resulting image, and combining the
%   results with an OR. Finally, we run a hole filling algorithm in case
%   any voxels within the mesh were missed.
%
%   DILRAD is the radius of the dilation algorithm that is used internally
%   to fix artifacts by itk_tri_rasterization(). Larger values of DILRAD
%   are a safer bet that artifacts are removed, but it also makes the
%   function slower. By default DILRAD=10.
%
%
% This function uses alphavol() by Jonas Lundgren, which is faster than the
% alpha shape functions in the CGAL library because the latter need to use
% a rather slow Delaunay triangulation first.
%
% See also: cgal_insurftri, alphavol, itk_tri_rasterization,
% cgal_alpha_shape3, cgal_fixed_alpha_shape3, scimat_closed_surf_to_bw.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012-2013 University of Oxford
% Version: 0.6.0
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
narginchk(2, 4);
nargoutchk(0, 3);

% if the image has all voxels == 1, then the smoothed local convex hull is
% the image itself and we don't need to waste time doing other computations
%
% likewise, if the image has all voxels == 0, there's nothing to smooth
if ((nnz(scimat.data) == 0) || (numel(scimat.data) == nnz(scimat.data)))
    return
end

% if we don't have at least 4 unique points, the delaunay triangulation
% will give an error
if (nnz(scimat.data) < 4)
    return
end

% defaults
if (nargin < 3 || isempty(SAFE))
    SAFE = true;
end
if (nargin < 4 || isempty(dilrad))
    dilrad = 10;
end

% compute perimeter voxels
scimat.data = bwperim(scimat.data);

% coordinates of segmented voxels
[r, c, s] = ind2sub(size(scimat.data), find(scimat.data));
x = scinrrd_index2world([r c s], scimat.axis);
clear r c s

% compute alpha shape
try
    [~, s] = alphavol(x, rad);
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

% correct non-manifold orientation of triangles, non-connected
% vertices, etc
[tri, x] = tri_squeeze(tri, x);
[x, tri] = meshcheckrepair(x, tri, 'deep'); % correct orientation
tri = cgal_tri_fillholes(tri, x); % fill holes in the surface

% % compute alpha shape (alternative mode using CGAL functions): this method
% % is slower because of the CGAL Delaunay triangulation
% tic
% tri = cgal_fixed_alpha_shape3(x, rad.^2);
% tri = tri{1};
% toc

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

if (SAFE)
    
    % the safe way to rasterize the mesh interior is by using the CGAL
    % function. This is quite slow, though

    % init output
    scimat.data(:) = 0;
    
    % bounding box obtained from the mesh, rounded up outwards to the
    % closest voxel centres
    xmin = min(x);
    xmax = max(x);
    idxmin = floor(scinrrd_world2index(xmin, scimat.axis));
    idxmax = ceil(scinrrd_world2index(xmax, scimat.axis));
    xmin = scinrrd_index2world(idxmin, scimat.axis);
    xmax = scinrrd_index2world(idxmax, scimat.axis);
    
    % vectors of voxel centres to check whether they are inside or outside
    % the mesh
    cx = linspace(xmin(1), xmax(1), idxmax(2) - idxmin(2) + 1);
    cy = linspace(xmin(2), xmax(2), idxmax(1) - idxmin(1) + 1);
    cz = linspace(xmin(3), xmax(3), idxmax(3) - idxmin(3) + 1);
    
    % check whether voxels are inside or outside
    tic
    scimat.data(idxmin(1):idxmax(1), ...
        idxmin(2):idxmax(2), ...
        idxmin(3):idxmax(3)) ...
        = cgal_insurftri(tri, x, {cx, cy, cz});
    toc
    
else
    
    % rasterize mesh to binary segmentation. This is fast, but quite often
    % voxels that should be labelled as 1 get labelled as 0
    res = [scimat.axis.spacing]; % (r, c, s) format
    sz = size(scimat.data); % (r, c, s) format
    origin = scinrrd_index2world([1, 1, 1], scimat.axis);
    scimat.data = itk_tri_rasterization(tri, x, res, sz, origin);
    
    % that's why we re-compute the rasterization by permuting x and y
    % dimensions, and add the voxels found here to the voxels found before
    aux = itk_tri_rasterization(...
        tri, ...
        x(:, [2 1 3]), ...
        res([2 1 3]), ...
        sz([2 1 3]), ...
        origin([2 1 3]));
    aux = permute(aux, [2 1 3]);
    scimat.data = scimat.data | aux;
    
    % and again, permuting x <-> z
    aux = itk_tri_rasterization(...
        tri, ...
        x(:, [2 3 1]), ...
        res([3 1 2]), ...
        sz([3 1 2]), ...
        origin([2 3 1]));
    aux = permute(aux, [2 3 1]);
    scimat.data = uint8(scimat.data | aux);
    
    %% correct segmentation
    
    if (dilrad == 0)
        return
    end
    
    % the approach above is fast, but it often fails to set to 1 voxels within
    % the mesh. What we do is dilate the segmentation, fill holes, and then
    % re-check all the resulting extra voxels
    
    % dilate segmentation
    aux = itk_imfilter('bwdilate', scimat.data, dilrad, 1);
    [r, c, s] = ind2sub(size(aux), find(aux));
    if (any(r == 1) || any(r == size(aux, 1)) ...
            || any(c == 1) || any(c == size(aux, 2)) ...
            || any(s == 1) || any(s == size(aux, 3)) ...
            )
        warning('Gerardus:ImageOverflow', 'Dilation radius too large for image boundaries')
    end
    
    % bounding box obtained from the mesh
    xmin = min(x);
    xmax = max(x);
    idxmin = floor(scinrrd_world2index(xmin, scimat.axis));
    idxmax = ceil(scinrrd_world2index(xmax, scimat.axis));
    
    % remove voxels outside the bounding box, as we know for sure they are
    % part of the background
    aux(1:idxmin(1), :, :) = 0;
    aux(idxmax(1):end, :, :) = 0;
    aux(:, 1:idxmin(2), :) = 0;
    aux(:, idxmax(2):end, :) = 0;
    aux(:, :, 1:idxmin(3)) = 0;
    aux(:, :, idxmax(3):end) = 0;
    
    % fill holes in the segmentation
    aux = imfill(aux, 'holes');
    
    % get list of voxels that the dilation has added and their coordinates
    idx = find(xor(aux, scimat.data));
    [r, c, s] = ind2sub(size(aux), idx);
    xi = scinrrd_index2world([r, c, s], scimat.axis);
    
    % check those extra voxels with
    scimat.data(idx) = cgal_insurftri(tri, x, xi, rand(3));

end
