function rotmat = quaternion2matrix(q)
% QUATERNION2MATRIX  Convert quaternion to rotation matrix
%
% ROTMAT = QUATERNION2MATRIX(Q)
%
%   Q is a 4-vector encoding a rotation in 3D space.
%
%   ROTAM is the rotation matrix that implements the same rotation as the
%   quaternion.
%
%   See e.g. http://www.sjbrown.co.uk/2002/05/01/quaternions/

% Author: Vicente Grau
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
error( nargchk( 1, 1, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

% check length of input vector
if ( length( q ) ~= 4 )
    error( 'Q must be a 4-vector' )
end

% compute rotation matrix from quaternion
rotmat = [ 1 - 2*q(3).^2 - 2*q(4).^2,  ...
    2*q(2)*q(3) - 2*q(1)*q(4), ...
    2*q(4)*q(2) + 2*q(1)*q(3); ...
    ...
    2*q(2)*q(3) + 2*q(1)*q(4), ...
    1 - 2*q(2).^2 - 2*q(4).^2, ...
    2*q(3)*q(4) - 2*q(1)*q(2); ...
    ...
    2*q(4)*q(2) - 2*q(1)*q(3), ...
    2*q(3)*q(4) + 2*q(1)*q(2), ...
    1 - 2*q(2).^2 - 2*q(3).^2 ];
