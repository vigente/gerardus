% test_dijkstra.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012 University of Oxford
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

d = [...
    %  1   2   3   4   5   6
       0   7   9   0   0  14; ...
       7   0  10  15   0   0; ...
       9  10   0  11   0   2; ...
       0  15  11   0   6   0; ...
       0   0   0   6   0   9; ...
      14   0   2   0   9   0];
  
% run disjkstra's algorithm on first node
[dd, p] = dijkstra(sparse(d), 1)

% dd =
% 
%     0.0000    7.0000    9.0000   20.0000   20.0000   11.0000
% 
% p =
% 
%      0     1     1     3     6     3

% query several source nodes at the same time
[dd, p] = dijkstra(sparse(d), [1 3 6])

% dd =
% 
%     0.0000    7.0000    9.0000   20.0000   20.0000   11.0000
%     9.0000   10.0000    0.0000   11.0000   11.0000    2.0000
%    11.0000   12.0000    2.0000   13.0000    9.0000    0.0000
% 
% p =
% 
%      0     1     1     3     6     3
%      3     3     0     3     6     3
%      3     3     6     3     6     0

% query a source node, but we are only interested in the path to a single
% target node
[dd, p] = dijkstra(sparse(d), 1, 3)

% dd =
% 
%     0.0000    7.0000    9.0000   20.0000       Inf   11.0000
% 
% p =
% 
%      0     1     1     3   NaN     3
