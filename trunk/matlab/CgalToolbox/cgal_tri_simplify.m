function varargout = cgal_tri_simplify(varargin)
% CGAL_TRI_SIMPLIFY  Reduce number of faces in triangular mesh using edge
% collapse.
%
% This function provides a Matlab interface to Fernando Cacciola's edge
% collapse method to simplify a triangulated surface mesh available in
% CGAL, using the default Lindstrom-Turk cost-strategy:
%
% http://doc.cgal.org/latest/Surface_mesh_simplification/index.html
%
% From the documentation: "Roughly speaking, the method consists of
% iteratively replacing an edge with a single vertex, removing 2 triangles
% per collapse."
%
% "Edges are collapsed according to a priority given by a user-supplied
% cost function, and the coordinates of the replacing vertex are determined
% by another user-supplied placement function. The algorithm terminates
% when a user-supplied stop predicate is met, such as reaching the desired
% number of edges."
%
% [TRI2, X2] = cgal_tri_simplify(TRI, X, RATIO)
%
%   TRI is a 3-column matrix. Each row contains the 3 nodes that form one
%   triangular facet in the mesh.
%
%   X is a 3-column matrix. X(i, :) contains the xyz-coordinates of the
%   i-th node in the mesh.
%
%   RATIO is a scalar with the stop criterion. The algorithm will stop when
%   the number of current undirected edges is RATIO * number of original
%   undirected edges.
%
%   TRI2, X2 is the description of the simplified output mesh.
%
% See also: cgal_surfsubdivision.

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

error('MEX function not found')
