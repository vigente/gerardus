function d = dmatrix_complete(d)
% DMATRIX_COMPLETE  Fill in empty elements of non-Euclidean distance matrix
%
% D2 = DMATRIX_COMPLETE(D1)
%
%   D1 is a sparse distance matrix that represents the connectivity of a
%   set of points in a non-Euclidean space.
%
%   That is, not all points are connected to each other. (Note that lack of
%   connection in the sparse matrix is represented by a 0 element).
%
%   D2 is a full matrix with the distances between every pair of points.
%   For points that are not directly connected, their distance is computed
%   using Dijkstra's algorithm over the connection graph.
%
% This function uses a very fast implementation of Dijkstra's algorithm
% written by Mark Steyvers.

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
error(nargoutchk(0, 1, nargout, 'struct'));

% number of rows in the distance matrix
N = size(d, 1);
if (N ~= size(d, 2))
    error('Distance matrix must be square')
end

% compute shortest distances from each node in the graph
d = dijkstra(d, 1:N);
