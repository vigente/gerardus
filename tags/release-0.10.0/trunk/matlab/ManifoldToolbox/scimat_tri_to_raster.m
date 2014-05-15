function scimat = scimat_tri_to_raster(tri, x, scimat, raster)
% SCIMAT_TRI_TO_RASTER  Rasterize a closed triangular mesh
%
% This takes a triangular mesh and returns a raster image with voxels == 1
% inside the mesh and == 0 outside the mesh.
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
%   SCIMAT is the scimat-format struct with the details about the number of
%   voxels, voxel size and offset of the output raster image.
%
%   SCIMAT2 is the output, where the voxels within the closed surface have
%   been set to 1. Note that these voxels are _added_ to whatever previous
%   segmentation was present in SCIMAT.
%
% SCIMAT2 = scimat_closed_surf_to_bw(..., RASTER)
%
%   RASTER is a string with the rasterization method:
%
%     'cgal_insurftri' (default): More reliable, although much slower.
%
%     'itk_tri_rasterization': Fast, although sometimes lines or regions of
%     voxels that should be set to 1 are set to 0.
%
%
% See also: cgal_insurftri, itk_tri_rasterization,
% scimat_lconvhull_smoothing.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
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
narginchk(3, 4);
nargoutchk(0, 1);

if (size(tri, 2) ~= 3)
    error('TRI must be a 3-column matrix')
end
if (size(x, 2) ~= 3)
    error('X must be a 3-column matrix')
end

% defaults
if (nargin < 4 || isempty(raster))
    raster = 'cgal_insurftri';
end

% very small segmentations can produce degenerated surfaces. In that case,
% there's nothing to smooth and we can just exit without modifying the
% input segmentation
if (size(tri, 1) < 4)
    return
end

switch raster
    
    case 'cgal_insurftri'
        
        % the safe way to rasterize the mesh interior is by using the CGAL
        % function. This is quite slow, though
        
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
        scimat.data(idxmin(1):idxmax(1), ...
            idxmin(2):idxmax(2), ...
            idxmin(3):idxmax(3)) ...
            = cgal_insurftri(tri, x, {cx, cy, cz}, rand(3, 3));
        
    case 'itk_tri_rasterization'
        
        % rasterize mesh to binary segmentation. This is faster, although in
        % some situations we have observed some artifacts in the resulting
        % segmentation in the form of lines or regions of voxels == 0 that
        % should be == 1
        res = [scimat.axis.spacing]; % (r, c, s) format
        sz = size(scimat.data); % (r, c, s) format
        origin = scinrrd_index2world([1, 1, 1], scimat.axis);
        scimat.data = itk_tri_rasterization(tri, x, res, sz, origin);
        
    otherwise
        
        error('Gerardus:WrongInputValue', 'Invalid rasterizer')
        
end
