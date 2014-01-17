function [tri, x, map] = tri_squeeze(tri, x)
% TRI_SQUEEZE  Remove disconnected vertices from triangular mesh
%
% [TRI2, X2, I] = tri_squeeze(TRI, X)
%
%   TRI is a 3-column matrix. Each row contains the 3 nodes that form one
%   triangular facet in the mesh.
%
%   X is a 3-column matrix. X(i, :) contains the xyz-coordinates of the
%   i-th node in the mesh.
%
%   Assuming that TRI, X represents a triangular mesh where some of the
%   vertices in X are not part of any triangle in TRI, this function
%   returns TRI2, X2, that represent the same mesh, but with the
%   disconnected vertices removed.
%
%   I is the mapping of preserved vertices such that X2 = X(I), given as a
%   list of indices.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012-2014 University of Oxford
% Version: 0.2.0
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

% check arguments
narginchk(2, 2);
nargoutchk(0, 3);

% number of vertices
N = size(x, 1);

% list of all vertices that are part of at least a triangle
idx = unique(tri(:));

% vector that will do the mapping
map = zeros(N, 1);
map(idx) = (1:length(idx))';

% rename vertices
tri = map(tri);

% delete disconnected vertices
x(~map, :) = [];

if (nargout > 2)
    % format the mapping of preserved vertices as indices instead of a
    % boolean vector
    map = find(map);
end
