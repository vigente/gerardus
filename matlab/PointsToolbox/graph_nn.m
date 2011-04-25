function nn = graph_nn(d, idxfrom, idxto)
% GRAPH_NN  Nearest neighbours from a subset of nodes in a graph to another
% subset
%
% This function establishes a corresponde between two subsets of nodes in a
% graph, so that each node in the first subset corresponds to the closest
% one in the second subset.
%
% NN = GRAPH_NN(D, IDXFROM, IDXTO)
%
%   D is a graph expressed as a sparse or full distance matrix, where
%   D(i,j) is the distance between the i-th and j-th nodes.
%
%   IDXFROM and IDTO are vectors that represent two subsets of nodes in the
%   graph. By default, the vectors contain every node in the graph. So if
%   you want to compute the nearest neighbours from all the nodes to a
%   subset, you can run
%
%     >> nn = graph_nn(d, [], idxto);
%
%   NN is a vector with the same length as IDXFROM. NN(i)==j means that j
%   is the closest node to node i.
%
% This function uses a very fast C++ implementation of the Dijkstra
% algorithm by Mark Steyvers.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.1.1
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
error(nargchk(1, 3, nargin, 'struct'));
error(nargoutchk(0, 2, nargout, 'struct'));

% default
if (nargin < 2 || isempty(idxfrom))
    idxfrom = 1:size(d, 1);
end
if (nargin < 3 || isempty(idxto))
    idxto = 1:size(d, 1);
end

% compute minimum distances from the skeleton points to every other node in
% the graph
mind = dijkstra(d, idxto);

% we want each point in the segmentation linked only to one point in the
% skeleton. This is the closest point in the skeleton (nnidx has the index
% of the nearest neighbour in the skeleton to each segmentation point)
[~, nn] = min(mind, [], 1);

% keep only nodes we are interested in
nn = nn(idxfrom);
