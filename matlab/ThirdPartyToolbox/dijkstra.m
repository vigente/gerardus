function [d, p] = dijkstra(g, s)
% DIJKSTRA  Shortests path tree from sparse matrix graph
%
% [D, P] = DIJKSTRA(G, S)
%
%   This function computes the shortest distance tree(s) from one (or more)
%   source node(s) in a graph.
%
%   G is a sparse matrix where element (i,j) contains the cost of going
%   from node i to node j in the graph. That is, G contains the edge
%   weights or distances between nodes.
%
%   Note: Matlab's convention is that the absence of an element in a sparse
%   matrix means edge weight = 0. However, in this function, the lack of an
%   entry in the sparse matrix is understood as a lack of edge (or edge
%   weight = Inf).
%
%   S is a vector with a list of source node indices. A shortest path tree
%   will be computed for each one of the nodes in S.
%
%   D, P are matrices where each row corresponds to a source node in S. D
%   is the shortest distance from each node to the source node.
%
%   P is the predecessor of each node in the shortest path tree.
%
%   Because this function is implemented as a MEX function in C++,
%   takes a sparse matrix at the input, and uses a Fibonacci heap
%   implementation [1], it is well suited for huge sparse graphs.
%
%   [1]
%   http://www.ahhf45.com/info/Data_Structures_and_Algorithms/resources/technical_artile/fibonacci_heap/fibonacci.htm

% Author: Mark Steyvers, Stanford University, 19 Dec 2000.
% Version: 0.1.1
% $Rev$
% $Date$
%
% Modified by Ramon Casero <rcasero@gmail.com>, University of Oxford,  23
% Mar 2010 to also provide the predecessor list at the output.
%
% This file is distributed as a derivative work of a third-party function
% with project Gerardus.
%
% http://code.google.com/p/gerardus/

error('Compiled MEX function has not been found')
