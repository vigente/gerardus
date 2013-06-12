function [idx, d, p] = cgal_closest_trifacet(tri, x, xi)
% CGAL_CLOSEST_TRIFACET  Closest triangular facet of a mesh to a point in 3D
%
%  Note that this C++ MEX implementation is 54 times faster than the
%  Matlab implementation closest_trifacet() in the Right Ventricle example
%  of matlab/test/test_cgal_closest_trifacet.m.
%
% [IDX, D, P] = cgal_closest_trifacet(TRI, X, XI)
%
%   TRI is a 3-column matrix. Each row contains the 3 nodes that form one
%   triangular facet in the mesh.
%
%   X is a 3-column matrix. X(i, :) contains the xyz-coordinates of the
%   i-th node in the mesh.
%
%   XI is a 3-column matrix. XI(i, :) are the xyz-coordinates of a test
%   point. This function finds the facet TRI(j,:) that is closest to
%   XI(i,:).
%
%   IDX is a vector with one element per point in XI. IDX(i) is the index
%   of the closest facet to XI(:,i). For instance, to obtain the nodes of
%   the facet closest to XI(i,:), run TRI(IDX(i), :).
%
%   D is a vector with the same length as IDX. D(i) is the distance of
%   point XI(i,:) to the mesh, or equivalently, to the closest facet.
%
%   P is a 3-column matrix. P(i,:) are the xyz-coordinates of the closest
%   point in the mesh to XI(i,:).
%
% See also: closest_trifacet (an inefficient Matlab implementation that
% mirrors this function)

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

error('MEX file not found')
