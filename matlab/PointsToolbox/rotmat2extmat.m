function b = rotmat2extmat(m, a, center)
% ROTMAT2EXTMAT  Convert centroid and rotation matrix pair to extended
% matrix form
%
% B = ROTMAT2EXTMAT(M, A)
%
%   A is a (3,3)-matrix that describes a rotation around a centroid M,
%   where M is a column 3-vector.
%
%     Y = A*(X - M) + M
%
%   B2 is a (4,4)-matrix that describes the same transformation in extended
%   form
%
%     [Y] = B2 [X]
%     [1]      [1]
%
%   B is B2 minus the bottom row.
%
% B = ROTMAT2EXTMAT(M, A, CENTER)
%
%   CENTER is a boolean (default CENTER=false). If true, the transformation
%   will move the centroid to 0:
%
%     Y = A*(X - M)
%
%
% See also: extmat2rotmat.

% Author: Ramon Casero
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
error( nargchk( 0, 3, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

% defaults
if (nargin < 1 || isempty(m))
    m = [0 0 0]';
end
if (nargin < 2 || isempty(a))
    a = eye(3);
end
if (nargin < 3 || isempty(center))
    center = false;
end

% make sure centroid is a column vector
m = m(:);

% extended format matrix
if (center)
    b = [a, -a*m];
else
    b = [a, (eye(3)-a)*m];
end
