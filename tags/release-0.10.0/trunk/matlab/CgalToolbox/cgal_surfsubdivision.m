function varargout = cgal_surfsubdivision(varargin)
% CGAL_SURFSUBDIVISION  3D Surface Approximating Subdivision Methods
%
% This function is a Matlab wrapper of the CGAL 3D Surface Subdivision
% Methods by Le-Jeng Andy Shiue.
%
% http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Subdivision_method_3/Chapter_main.html
%
% [TRI2, X2] = cgal_surfsubdivision(TRI, X, METHOD, ITER)
%
%   TRI is a 3-column matrix. Each row represents the indices of the three
%   vertices that form a triangle. TRI as a whole represents the closed
%   surface.
%
%   X is a 3-column matrix. Each row represents the Cartesian coordinates
%   of a vertex on the surface, indexed by TRI values.
%
%   METHOD is a string to select one of the following subdivision schemes
%   (all of them are approximating, not interpolating):
%
%     'CatmullClark': Catmull-Clark.
%     'Loop':         Loop
%     'DooSabin':     Doo-Sabin
%     'Sqrt3':        √3
%
%   ITER is an integer with the number of iterations to run. With ITER=0,
%   the output mesh is the same as the input mesh. With larger values of
%   ITER, the output mesh will be smoother, at the cost of having more
%   vertices. Usually, ITER=8 will produce already a very smooth mesh.
%
%   TRI2, X2 contain the output triangular mesh. Note that most subdivision
%   methods will produce some non-triangular facets, even when the starting
%   mesh is a triangulation. This function runs internally a triangulation
%   process to convert the generalized output polyhedron into a triangular
%   mesh.
%
% See also: cgal_tri_simplify.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
% Version: 0.1.1
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
