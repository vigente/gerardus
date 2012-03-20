function a = vec2rotmat(v, a0)
% VEC2ROTMAT  Compute rotation matrix from Cartesian vector
%
% A = VEC2ROTMAT(V)
%
%   V is a 3-vector. We assume that V is a rotation of the Z-axis vector
%   [0,0,1]'.
%
%   A is one of the infinite rotation matrices that rotate the Z-axis
%   vector onto V
%
%     V = A * [0,0,1]'
%
%   In particular, A is the solution that maps the Cartesian coordinate
%   system onto a coordinate system that contains V and where the other two
%   vectors are closest to A in a Gram-Schmidt sense.
%
% A = VEC2ROTMAT(V, A0)
%
%   A0 is a (3,3)-matrix that can be seen as a rotation matrix in 3D, as
%   an orthonormal basis, or a first guess of A. Only the first two columns
%   are used. The two first columns are used together with V as the initial
%   non-orthogonal basis. When A0 is not provided, then the identity matrix
%   is used, i.e.g A0 = eye(3).

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
error( nargchk( 1, 2, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

% defaults
if ( nargin < 2 || isempty( a0 ) )
    a0 = eye( 3 );
end

% check that V is a vector
if ( size( v, 2 ) ~= 1 || size( v, 1 ) ~= 3 )
    error( 'V must be a 3-vector' )
end

% create a non-orthogonal Cartesian basis where the Z-axis vector is
% replaced by v
v = [ a0( :, 1:2 ), v ];

% now we want to apply Gram-Schmidt to this basis to create an orthonormal
% basis. We want to keep v, so we start the Gram-Schmidt process fixing it,
% and then remove the projections of the other two vectors
[a, foo] = qr( v( :, end:-1:1 ) ); % foo is necessary, because the result 
                                   % from a = qr() is different
a = a( :, end:-1:1 );

% check whether we want to keep each vector or the negative, so that they
% are closest to the original basis
for I = 1:3
    
    if ( dot( v(:,I), a(:,I) ) < 0 )
        a(:,I) = -a(:,I);
    end
    
end
