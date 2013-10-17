function idx = tri_find_nonmanifold_vertex(tri, x, scimat_axis)
% TRI_FIND_NONMANIFOLD_VERTEX  Find indices of non-manifold vertices in a
% triangular mesh
%
% IDX = tri_find_nonmanifold_vertex(TRI, RCS)
%
%   This function is a hack so that we can identify which vertices
%   meshcheckrepair() identifies as non-manifold.
%
%   TRI is a 3-column matrix. Each row represents the indices of the three
%   vertices that form a triangle. TRI as a whole represents the closed
%   surface.
%
%   RCS is a 3-column matrix. Each row represents the (row, column, slice)
%   coordinates of a vertex in the mesh.
%
%   IDX is a vector with the indices of the vertices in X that are
%   non-manifold.
%
% IDX = tri_find_nonmanifold_vertex(TRI, X, SCIMAT.AXIS)
%
%   Alternatively, vertex coordinates X can be provided as real world
%   (x,y,z) Cartesian coordinates. SCIMAT.AXIS is the axis field from a
%   SCIMAT structure with the offset and voxel size information.
%
% See also: scinrrd_world2index, scinrrd_index2world.

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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

% check arguments
narginchk(2, 3);
nargoutchk(0, 1);

% defaults
if (nargin < 3 || isempty(scimat_axis))
    scimat = scinrrd_im2nrrd(zeros(2), [1 1 1], [0 0 0]);
    scimat_axis = scimat.axis;
    clear scimat
end

% we convert the vertex coordinates to indices to avoid loss of precision
% when meshcheckrepair() writes to and reads from OFF file
rcs = round(scinrrd_world2index(x, scimat_axis));

% meshcheckrepair() is supposed to be used to repair a non-manifold mesh.
% What it actually does is that it duplicates non-manifold vertices, which
% I'm not sure it's a good repair. Anyway, here we are interested only in
% detecting the non-manifold vertices. We are going to do it in a bit of a
% convoluted way, because meshcheckrepair() just duplicated the vertices,
% but doesn't tell which ones they are

% run meshcheckrepair(). This will duplicate any non-manifold vertices We
% use evalc so that we can supress the output to the screen.
[~, rcs2] = evalc('meshcheckrepair(rcs, tri, ''deep'');');

% round values to avoid finite precision errors
rcs2 = round(rcs2);

% assign one index value to each unique vertex
[~, idx] = unique(rcs2, 'rows'); 

% delete all vertices except for the duplicates
rcs2(idx, :) = []; 

% sometimes a vertex may be repetead more than once. We want a list of
% unique vertices
rcs2 = unique(rcs2, 'rows'); 

% indices of non-manifold vertices
idx = find(ismember(rcs, rcs2, 'rows'));
