function z = coords_from_dist_gower(x, d, tol)
% COORDS_FROM_DIST_GOWER  Compute point coordinates from distances to a set
% of landmarks, as described in Gower (1986)
%
% Z = COORDS_FROM_DIST_GOWER(X, D)
%
%   For 2D points, X is a (2, P)-matrix with the coordinates of P
%   landmarks. Points can have any dimension, they are not limited to being
%   2D. Note that if the points in X are collinear, the result Z will be
%   wrong in general, as the solution is not unique.
%
%   D is a P-vector with the distance of point Z to each of the landmarks.
%   D can also be a matrix, if each column is a vector of distances.
%
%   Gower(1986) showed that the coordinates of Z can be computed from X and
%   D.
%
% Z = COORDS_FROM_DIST_GOWER(X, D, TOL)
%
%   TOL is the tolerance value for the collinearity test. By default,
%   TOL=1e-16. If TOL=[], then the collinearity test is skipped.
%
% J.C. Gower. Adding a point to vector diagrams in multivariate analysis.
% Biometrika, 55(3):582–585, 1968.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010-2011 University of Oxford
% Version: 0.2.0
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
error(nargchk(2, 3, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

if (size(d, 1) ~= size(x, 2))
    error('Vertical distance vector d must have a value for each landmark')
end

% defaults
if (nargin < 3)
    tol = 1e-16;
end

% test for collinearity (skip test if tol=[])
if (~isempty(tol) && iscollinear(x, tol))
    error('X points are collinear, output coordinates cannot be determined')
end

% compute centroid of landmarks
xmean = mean(x, 2);

% center landmarks
x = x - xmean(:, ones(1, size(d, 1)));

% compute landmarks squared norm
xnorm2 = sum(x.^2, 1)';

% replicate squared norm vector
xnorm2 = xnorm2(:, ones(1, size(d, 2)));

% compute point coordinates
z = xmean(:, ones(1, size(d, 2))) + .5 * pinv(x)' * (xnorm2 - d.^2);
