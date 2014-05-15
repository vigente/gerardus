function [N, cc] = graphcc(g)
% GRAPHCC  Split a graph into connected components.
%
% [N, CC] = GRAPHCC(G)
%
%   G is a square adjacency matrix of a graph. G(i,j)~=0 means that
%   vertices i, j are directly connected.
%
%   N is a scalar with the number of connected components in the graph.
%
%   CC is a cell array with N elements. CC{i} is a list of all the vertices
%   that belong to the i-th connected components. Vertices that belong to
%   different connected components have no paths between them in the graph.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

% check arguments
narginchk(1, 1);
nargoutchk(0, 2);

% list of indices
idx = find(sum(g, 2));

cc = [];
while ~isempty(idx)
    
    % pop-out a vertex
    v = idx(1);
    idx(1) = [];
    
    % find vertices connected to this one
    dv = dijkstra(g, v);
    
    % put all the connected vertices in a connected component
    cc{end+1} = find(~isinf(dv));
    
    % remove those vertices from the list of vertices
    idx = setdiff(idx, cc{end});
    
end

N = length(cc);
