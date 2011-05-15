function [m, a] = extmat2rotmat(b, center)
% EXTMAT2ROTMAT  Convert rotation in extended matrix form into centroid and
% matrix pair
%
% [M, A] = EXTMAT2ROTMAT(B)
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
%   Warning: Computing M from B2 requires the inversion of a matrix that
%   can be ill-conditioned or singular. In that case, an approximation
%   proposed by [1] is used, but some errors observed are as large as
%   O(1e-3).
%
%  [1] H.J. Kim, et al., “A new algorithm for solving ill-conditioned
%  linear systems”, IEEE Trans. Magnetics, 33: 1373-1376, 1996
%
% [M, A] = EXTMAT2ROTMAT(B, CENTER)
%
%   CENTER is a boolean (default CENTER=false). If true, the transformation
%   that is assumed is:
%
%     Y = A*(X - M)
%
%
% See also: rotmat2extmat.

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
error( nargchk( 0, 2, nargin, 'struct' ) );
error( nargoutchk( 0, 2, nargout, 'struct' ) );

% defaults
if (nargin < 1 || isempty(b))
    b = eye(4);
end
if (nargin < 2 || isempty(center))
    center = false;
end

% shift constant
k = 1e-8;

% extract rotation matrix from extended matrix
a = b(1:3, 1:3);

% linear system that needs to be solved to compute the centroid
if (center)
    aux = -a;
else
    aux = eye(3)-a;
end

% compute translation vector from extended matrix
% if the linear system is ill-conditioned, we need to use a trick to
% solve it
if (cond(aux) > 1e8)
    
    warning('Ill-conditioned matrix, using Kim et al. (1996) approximation')
    
    % compute eigevalues of the matrix
    d = eig(aux);
    
    % compute another auxiliary matrix
    aux2 = inv((aux + k * eye(3)));
    
    % obtain final approximation
    m = aux2 * (eye(3) + k * aux2) * b(1:3,4);
else
    m = aux \ b(1:3,4);
end
