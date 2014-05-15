function a = trifacet_signed_area(tri, x)
% TRIFACET_SIGNED_AREA  Signed area of the facets in a 2D triangulation.
%
% A = trifacet_signed_area(TRI, X)
%
%   TRI is a 3-column matrix. Each row contains the 3 nodes that form one
%   triangular facet in the mesh.
%
%   X is a 2-column matrix. X(i, :) contains the xy-coordinates of the
%   i-th node in the mesh.
%
%   A is a vector where A(i) is the signed area of the i-th triangle.
%   Positive areas correspond to counter-clockwise triangles, and negative
%   areas, to clockwise triangles.

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
narginchk(2, 2);
nargoutchk(0, 1);

if (size(tri, 2) ~= 3)
    error('TRI must have 3 columns')
end
if (size(x, 2) ~= 2)
    error('X must have 2 columns')
end

% signed area of a triangle with vertices  V1, V2, V3 in counter-clockwise
% orientation
%
% a = 1/2*((x2-x1)(y3-y1) - (x3-x1)(y2-y1))

x2x1 = x(tri(:, 2), 1) - x(tri(:, 1), 1);
y3y1 = x(tri(:, 3), 2) - x(tri(:, 1), 2);
x3x1 = x(tri(:, 3), 1) - x(tri(:, 1), 1);
y2y1 = x(tri(:, 2), 2) - x(tri(:, 1), 2);

a = 0.5 * (x2x1 .* y3y1 - x3x1 .* y2y1);
