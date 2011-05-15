function x = intersect_line_plane(mx, nx, p0, p1, p2)
% INTERSECT_LINE_PLANE  Find intersection of a line with a plane
%
% X = INTERSECT_LINE_PLANE(MX, NX, P0, P1, P2)
%
%   MX, NX are 3-vectors that represent, respectively, a point and the
%   direction vector of a line in 3D.
%
%   P0, P1, P2 are 3-vectors that represent 3 non-collinear points on a
%   plane.
%
%   X is a 3-vector with the coordinates of the intersection between the
%   line and the plane.
%
%   All vectors are column vectors.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010 University of Oxford
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
error( nargchk( 5, 5, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

% make sure that inputs are column vectors
mx = mx(:);
nx = nx(:);
p0 = p0(:);
p1 = p1(:);
p2 = p2(:);

% equation to find the intersection point obtained from Wikipedia:
% http://en.wikipedia.org/wiki/Line-plane_intersection
k = [ -nx p1-p0 p2-p0 ] \ (mx-p0);
x = mx + nx * k(1);
