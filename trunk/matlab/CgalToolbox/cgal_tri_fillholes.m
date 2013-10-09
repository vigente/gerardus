function varargout = cgal_tri_fillholes(varargin)
% CGAL_TRI_FILLHOLES  Fill holes in a triangular mesh, with triangles.
%
% This function fills any holes in a triangular mesh. As holes are not
% necessarily triangular, any non-triangular facets are split into
% triangles. The resulting mesh is both closed and triangular.
%
% [TRI2, N] = cgal_tri_fillholes(TRI, X)
%
%   TRI is a 3-column matrix. Each row represents the indices of the three
%   vertices that form a triangle. TRI as a whole represents the closed
%   surface.
%
%   Note that this function requires that the triangles in TRI have a
%   manifold orientation. This can be achieved pre-processing with
%
%     [x, tri] = meshcheckrepair(x, tri, 'deep'); % correct orientation
%
%   X is a 3-column matrix. Each row represents the Cartesian coordinates
%   of a vertex on the surface, indexed by TRI values.
%
%   TRI2 is TRI plus the extra facets that fill the holes.
%
%   N is a scalar with the number of holes that have been filled. N <=
%   number of triangles added to TRI2.
%
% See also: meshcheckrepair.

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

error('MEX file not found')
