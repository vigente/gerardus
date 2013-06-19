% test_dijkstra.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012 University of Oxford
% Version: 0.3.1
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic Dijkstra

% distance matrix to describe the graph in 
% http://en.wikipedia.org/wiki/Dijkstra's_algorithm
% (saved to matlab/test/Dijkstra_Animation_wikipedia.gif)

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
%      0     7     9    20    20    11
% 
% p =
% 
%      0     1     1     3     6     3

% query several source nodes at the same time
[dd, p] = dijkstra(sparse(d), [1 3 6])

% dd =
% 
%     0     7     9    20    20    11
%     9    10     0    11    11     2
%    11    12     2    13     9     0
% 
% p =
% 
%      0     1     1     3     6     3
%      3     3     0     3     6     3
%      3     3     6     3     6     0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dijkstra with limited number of targets

% query a source node, but we are only interested in the path to a single
% target node
[dd, p] = dijkstra(sparse(d), 1, 3)

% dd =
% 
%     0     7     9    20   Inf    11
% 
% p =
% 
%      0     1     1     3   NaN     3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dijkstra with piggyback

% secondary distance matrix. The algorithm will ignore this distances in
% terms of searching for shortest-paths, but will create another distance
% output adding the values in D2 while following the path given by D

d2 = [...
    %  1   2   3   4   5   6
       0  10   7   0   0   6; ...
      10   0   2  11   0   0; ...
       7   2   0   9   0   9; ...
       0  11   9   0  14   0; ...
       0   0   0  14   0  15; ...
       6   0   9   0  15   0];

[dd, p, dd2] = dijkstra(sparse(d), 1:6, 1:6, sparse(d2))

% dd2 =
% 
%     0    10     7    16    31    16
%    10     0     2    11    26    11
%     7     2     0     9    24     9
%    16    11     9     0    14    18
%    31    25    24    14     0    15
%    16    11     9    18    15     0
   