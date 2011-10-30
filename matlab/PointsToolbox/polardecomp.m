function [p, u] = polardecomp(a)
% POLARDECOMP  Polar decomposition of a square complex matrix
%
% [P, U] = POLARDECOMP(A)
%
% (From Wikipedia [1]) The polar decomposition of a square complex matrix A
% is a matrix decomposition of the form
%
%    A = U*P = P'*U
%
% where U is a unitary matrix and P is a positive-semidefinite Hermitian
% matrix. Intuitively, the polar decomposition separates A into a component
% that stretches the space along a set of orthogonal axes, represented by
% P, and a rotation represented by U.
%
% This decomposition always exists; and so long as A is invertible, it is
% unique, with P positive-definite.
%
% In terms of the singular value decomposition of A, A = W*Σ*V', one has
%
%   P = V*Σ*V'
%   U = W*V'
%
% confirming that P is positive-definite and U is unitary.
%
% [1] http://en.wikipedia.org/wiki/Polar_decomposition

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
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
error(nargchk(1, 1, nargin, 'struct'));
error(nargoutchk(0, 2, nargout, 'struct'));

% polar decomposition
[w, s, v] = svd(a);
p = v * s * v';
u = w * v';
