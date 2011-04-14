function h = pts_rotmat(x1, x2)
% PTS_ROTMAT  Compute rotation matrix to minimize least squares error
%
% H = PTS_ROTMAT(X1, X2)
%
%    X1, X2 are two sets of points with dimensions (P, K) where P is the
%    number of points of dimension K.
% 
%    H is a rotation matrix with dimensions (K, K) of X2 onto X1, such that
%    the approximation between X1 and X2*H is optimal in the least squares
%    sense. Note: In 2D, the counterclockwise definition of H is
%
%      H = [cos(phi) sin(phi) ; -sin(phi) cos(phi)]
%
%    and the clockwise
%
%      H = [cos(phi) -sin(phi) ; sin(phi) cos(phi)]

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 1.0.0
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
error(nargchk(2, 2, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% compute rotation matrix
[u, sigma, v] = svd(x1' * x2);
h = v * sign(sigma) * u';
