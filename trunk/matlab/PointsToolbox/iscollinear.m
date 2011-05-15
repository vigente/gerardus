function t = iscollinear(x, tol)
% ISCOLLINEAR  Test to determine whether a set of points are collinear
%
% T = ISCOLLINEAR(X)
%
%   X is a (D, N)-matrix with N points of D dimension. This test works with
%   points of any dimension.
%
%   T is a boolean that says whether X is a set of collinear points.
%
% T = ISCOLLINEAR(X, TOL)
%
%   TOL is the tolerance value for collinearity. By default, TOL=1e-16.

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
error(nargchk(1, 2, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% defaults
if (nargin < 2)
    tol = 1e-16;
end

% we need to have at least three points to test for collinearity, the other
% cases are trivial
if size(x, 2) == 0
    t = [];
    return
elseif size(x, 2) < 3
    t = true;
    return
end

% compute principal component analysis of the points
[~, d] = pts_pca(x);

% if there's only 1 eigenvalue, we know that the set is collinear
if (length(d) == 1)
    t = true;
    return
end

% normalize second eigenvalue by first one
t = sqrt(d(2) / d(1)) < tol;
