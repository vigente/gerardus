function [x, s] = pts_cn(x)
% PTS_CN  Center and normalize sets of points
%
% [Y, S] = PTS_CN(X)
%
%    X is a volume of points, i.e. a (P, K, N)-matrix where P is the number
%    of points of K dimensions, and N is the number of sets of points.
%
%    Y is a matrix like X, but with all sets centered around the origin and
%    with the square root of the sum of the squared distances of the points
%    to the centroid equal to 1.
%
%    S is a column vector with N elements. Each one is the scaling factor
%    that has been applied to one set of points.
%
%       Y = A * X * S
%
%    where A is EYE(P) - ONES(P) / P, as described in [1].
%
%  [1] Rohlf, F. & Slice, D. Extensions of the Procrustes method for the
%  optimal superimposition of landmarks. Systematic Zoology, 1990, 39,
%  40-59.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 1.0.0
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

% get sizes
[~, K, N] = size(x);

% init output
s = zeros(N, 1);

% centering transformation
for I = 1:N
    for J = 1:K
        x(:, J, I) = (x(:, J, I) - mean(x(:, J, I)));
    end
    s(I) = 1 / sqrt(sum(sum(x(:, :, I) .^ 2)));
    x(:, :, I) = x(:, :, I) * s(I);
end

% This is [1]'s elegant centering algorithm. However, it's awfully slow for
% even not very large volumes of data. And it needs loads of memory.
% 
% ip = eye(P) - ones(P) ./ P;
% 
% for I = 1:N
% 
%     % scale factors
%     s(I) = 1 / sqrt(trace(ip * x(:, :, I) * x(:, :, I)' * ip));
% 
%     % scaled and centered configurations
%     x(:, :, I) = ip * x(:, :, I) .* s(I);
%     
% end
