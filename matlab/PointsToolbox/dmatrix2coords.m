function v = dmatrix2coords(d, m)
% DMATRIX2COORDS  Point coordinates from distance matrix
%
%   This function computes the coordinates of a set of points if all we
%   know is the distances between them.
%
% X = DMATRIX2COORDS(D, M)
%
%   D is an (N,N)-square matrix, where D(i,j) is the Euclidean distance
%   between the i-th and j-th points.
%
%   M is a scalar with the dimension of the points (this information cannot
%   be extracted from the distance matrix). For example, if the points are
%   3D, then M=3. By default, M=N.
%
%   X is an (N, M)-matrix where each row has the coordinates of a point.
%   If you compute the distance matrix of V, you obtain D again.
%
%     >> D = dmatrix(X');
%
%   Note that any rotation or translation on X produces the same D. Thus, X
%   is not necessarily the same configuration that produced D originally.
%
%   Note also that due to numerical errors, it's possible that some
%   small complex values creep into the solution.
%
% This function implements the method in Gower (1966), that is better
% explained in Gower (1968).
%
% JC Gower, "Some distance properties of latent root and vector methods
% used in multivariate analysis", Biometrika, 53(3 and 4): 325-338, 1966.
%
% JC Gower, "Adding a point to vector diagrams in multivariate analysis",
% Biometrika, 55(3):582-585, 1968.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.1.0
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
error(nargchk(1, 2, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% number of points
N = size(d, 1);
if (N ~= size(d, 2))
    error('Distance matrix must be square')
end

% defaults
if (nargin < 2 || isempty(m))
    m = N;
end

% compute matrix A
A = -.5 * d.^2;

% compute matrix B
am = mean(A, 1);
amm = mean(am);
am = repmat(am, N, 1);
B = A - am - am' + amm;

% compute the eigenvalues and eigenvectors
[v, d] = eig(B);

% scale vectors
v = v * sqrt(d);

% remove extra dimensions
v = v(:, 1:m);
