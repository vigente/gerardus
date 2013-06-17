function [d, p, d2] = dijkstra(g, s, t, g2)
% DIJKSTRA  Shortests path tree from sparse matrix graph
%
% [D, P] = dijkstra(G, S)
%
%   This function computes the shortest distance tree(s) from one (or more)
%   source node(s) in a graph.
%
%   G is a sparse matrix where element (i,j) contains the cost of going
%   from node i to node j in the graph. That is, G contains the edge
%   weights or distances between nodes.
%
%   G doesn't need to be symmetric. The connections from node i are given
%   by column G(:, i). Thus, it's possible that e.g. G(2,5)=1.5, but
%   G(5,2)=0.0.
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
%   is the shortest distance from each target node to the source node.
%
%   P is the predecessor of each node in the shortest path tree.
%
% ... = dijkstra(..., T)
%
%   T is a row matrix with a list of targets. For each source node S(i),
%   the algorithm will stop when it finds the shortest paths to all targets
%   in T. This is useful to save time when G is a large matrix, but we are
%   interested only in e.g. the distance between to nodes that are nearby.
%   By default, T is a list of all nodes.
%
%   Note that when running in this mode, D, P will still return values for
%   all target nodes, but you can only trust D, P values of nodes that are
%   on paths between the source and target nodes.
%
% [..., D2] = dijkstra(..., T, G2)
%
%   This syntax implements our Dijkstra with piggyback extension. The
%   motivation is being able to find shortest-paths according to the
%   intensities in a greyscale image, but also compute the actual length of
%   those paths without having to backtrack the list of predecessors with
%   graphpred2path(). In that case, G contains intensity-weighted
%   distances, and is the graph used for finding shortest-paths. G2
%   contains actual Euclidean distances.
%
%   G2 must be a sparse matrix of type double and same size as G.
%
%   D2 is an output matrix that "gets a piggyback" from the algorithm, in
%   the sense that it contains shortest-distances computed from G2 over the
%   paths computed from G.
%   
%
%
%   Note:
%
%   Because this function is implemented as a MEX function in C++,
%   takes a sparse matrix at the input, and uses a Fibonacci heap
%   implementation [1], it is well suited for huge sparse graphs. According
%   to wikipedia [2], "The implementation based on a min-priority queue
%   implemented by a Fibonacci heap and running in O(|E|+|V|log|V|) is due
%   to (Fredman & Tarjan 1984). This is asymptotically the fastest known
%   single-source shortest-path algorithm for arbitrary directed graphs
%   with unbounded non-negative weights."
%
% [1] http://www.ahhf45.com/info/Data_Structures_and_Algorithms/resources/technical_artile/fibonacci_heap/fibonacci.htm
%
% [2] http://en.wikipedia.org/wiki/Dijkstra's_algorithm
%
% Fredman, Michael Lawrence; Tarjan, Robert E. (1984). "Fibonacci heaps and
% their uses in improved network optimization algorithms". 25th Annual
% Symposium on Foundations of Computer Science (IEEE): 338–346.
% doi:10.1109/SFCS.1984.715934

% Original Author: Mark Steyvers, Stanford University, 19 Dec 2000.
%
% Modified by Ramon Casero <rcasero@gmail.com>, University of Oxford,  
% 23 Mar 2010 Also provide the predecessor list at the output
%  6 Mar 2013 Accept a list of targets, to limit unnecessary computations
% 14 Mar 2013 Dijkstra with piggyback extension
%
% This file is distributed as a derivative work of a third-party function
% with project Gerardus.
%
% http://code.google.com/p/gerardus/
%
% Version: 0.3.1
% $Rev$
% $Date$

error('Compiled MEX function has not been found')
