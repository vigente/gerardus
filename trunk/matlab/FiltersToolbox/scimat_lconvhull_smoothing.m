function [tri, x, scimat] = scimat_lconvhull_smoothing(scimat, rad, raster)
% SCIMAT_LCONVHULL_SMOOTHING  Mesh and rasterization of the smoothing of a
% binary image using a local convex hull.
%
% [TRI, X] = scimat_lconvhull_smoothing(SCIMAT, RAD)
%
%   SCIMAT is a SCI MAT volume with a binary image. This binary image is
%   assumed to have one closed surface only.
%
%   An alpha-shape of the coordinates of all perimeter voxels is computed.
%   The boundary of the alpha-shape is a triangular mesh. Non-manifold
%   vertices are identified and removed iteratively, until all vertices are
%   manifold.
%
%   All vertices that were identified as non-manifold are removed from the
%   list of perimeter voxels, and the alpha-shape is recomputed. This
%   process is repeated iteratively until the result is an alpha-shape with
%   only manifold vertices.
%
%   RAD is a scalar with the radius of the alpha shape, i.e. the size of
%   the convex hull neighbourhood. When the convex hull is computed, only
%   voxels within RAD distance can be connected. When RAD=Inf, the
%   convex hull is obtained. Note that if you want to compare to other
%   alpha shape functions, e.g. cgal_alpha_shape3() or
%   cgal_fixed_alpha_shape3(), ALPHA=RAD^2.
%
%   The output is a triangular mesh described by TRI, X.
%
%   TRI is a 3-column matrix. Each row represents the indices of the three
%   vertices that form a triangle. TRI as a whole represents the closed
%   surface.
%
%   X is a 3-column matrix. Each row represents the Cartesian coordinates
%   of a vertex on the surface, indexed by TRI values.
%
% [TRI, X, SCIMAT2] = scimat_lconvhull_smoothing(SCIMAT, RAD, RASTER)
%
%   This syntax computes a binary rasterization of the mesh too. Note that
%   these voxels are _added_ to whatever previous segmentation was present
%   in SCIMAT.
%
%   RASTER is a string with the rasterization method:
%
%     'cgal_insurftri' (default): More reliable, although much slower.
%
%     'itk_tri_rasterization': Fast, although sometimes lines or regions of
%     voxels that should be set to 1 are set to 0.
%
%
% This function uses alphavol() by Jonas Lundgren, which is faster than the
% alpha shape functions in the CGAL library because the latter need to use
% a rather slow Delaunay triangulation first.
%
% See also: alphavol, cgal_insurftri, itk_tri_rasterization,
% cgal_alpha_shape3, cgal_fixed_alpha_shape3, scimat_tri_to_raster.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012-2013 University of Oxford
% Version: 0.7.0
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
narginchk(2, 3);
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
if (nargin < 3 || isempty(raster))
    raster = 'cgal_insurftri';
end

% compute perimeter voxels
scimat.data = bwperim(scimat.data);

% coordinates of segmented voxels
[r, c, s] = ind2sub(size(scimat.data), find(scimat.data));
x = scinrrd_index2world([r c s], scimat.axis);
clear r c s

% matrix to keep track of the non-manifold vertices that we are going
% to remove
xremoved = [];

% loop to remove potential troublesome vertices
while (true)
    
    % coordinates of segmented voxels
    [r, c, s] = ind2sub(size(scimat.data), find(scimat.data));
    x = scinrrd_index2world([r c s], scimat.axis);
    clear r c s
    
    % remove vertices that have been found to be non-manifold in the
    % past
    if ~isempty(xremoved)
        x = setdiff(x, xremoved, 'rows');
    end
    
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
    tri = s.bnd;
    clear s
    
    % remove vertices that are not connected to any triangles
    [tri, x] = tri_squeeze(tri, x);
    
    % number of non-manifold vertices so far
    nnonmani = size(xremoved, 1);
    
    % remove vertices until there are no non-manifold vertices
    while (true)
        
        % find non-manifold vertices
        idx = tri_find_nonmanifold_vertex(tri, x, scimat.axis);
        
        % stop if all vertices are now manifold
        if isempty(idx)
            break
        end
        
        % keep track of removed vertices
        xremoved = [xremoved; x(idx, :)];
        
        % remove non-manifold vertices
        [tri, x] = remove_vertex_from_tri(tri, x, idx);
        
    end
    
    % if no more non-manifold vertices have been found with respect to
    % the previous iteration, we have finished cleaning up the mesh
    if (nnonmani == size(xremoved, 1))
        break
    end
    
    
end

% DEBUG
% disp(['Number of non-manifoled vertices removed: ' num2str(nnonmani)])

% label connected components
[ncomp, trilab, xlab] = triconncomp(tri, x);

% DEBUG
% plot the mesh colouring triangles according to their label
% hold off
% trisurf(tri, x(:, 1), x(:, 2), x(:, 3), trilab)
% axis equal

% largest connected component
nbin = hist(xlab, 1:ncomp);
[~, idx] = max(nbin);

% vertices to remove
idx = find(xlab ~= idx);

% delete smaller components
[tri, x] = remove_vertex_from_tri(tri, x, idx);

% double-check that there aren't any non-manifold vertices
idx = tri_find_nonmanifold_vertex(tri, x, scimat.axis);
if ~isempty(idx)
    warning(['Assertion fail. Non-manifold vertices: ' num2str(length(idx))])
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

% if the user hasn't asked for the rasterization, we don't waste time
% computing it
if (nargout < 3)
    return
end

% rasterize the mesh
scimat = scimat_tri_to_raster(tri, x, scimat, raster);
