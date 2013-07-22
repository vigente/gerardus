function [u, v, stopCondition, err, dout] = lmdscale(d, u, v, opt)
% LMDSCALE  Local neighbourhood Multidimensional scaling on a plane
%
% This function solves a Multidimensional Scaling problem finding an
% embedding in R^2 of points according to a distance matrix D. Unlike
% classical MDS, D can be a sparse matrix, that is, each point's position
% only depends on a local neighbourhood of points connected to it.
%
% This function is based on the MDS optimisation idea proposed by Agarwal
% et al. (2010) for the case when the error metric is raw stress
%
%    sigma_r = sqrt(sum_{i,j} (d_ij - dout_ij)^2)
%
% We have improved it by using an exact one-step solution to the location
% of a point at a given distance from a set of neighbours by Gower (1968),
% replacing Agarwal's iterative PLACE() algorithm.
%
% [U, V, STOP, ERR, DOUT] = lmdscale(D)
%
%   D is a square matrix where D(i,j) is some distance (e.g. in meters)
%   between points i and j. Note that the points themselves are unknown, we
%   only know the distances between them. These distances could be geodesic
%   distances, or approximations thereof, e.g. chord-length distances.
%
%   D can be a full or sparse matrix. In both cases, values of D(i,j)=0 or
%   D(i,j)=Inf mean that points i and j are not neighbours, i.e. they are
%   not connected by an edge. Whether D is full or sparse can affect
%   performance. If D has lots of zeros, it's better to input D as a sparse
%   matrix. If D has few zeros, it's better to input D as a full matrix. 
%
%   The optimisation algorithm works by relocating each point i so that its
%   new distance on R^2 to each neighbour j is as close to D(i,j) as
%   possible. Note that non-neighbours (D=0 or D=Inf) are ignored in the
%   optimisation. This allows sparse matrices and local neighbourhood
%   optimisations.
%
%   U, V give an output parametrization of the unknown points, in the same
%   units as the distances. Note that U, V is one of inifite solutions, as
%   D is invariant to translations, rotations and reflections.
%
%   STOP is a cell-array with the condition/s that made the algorithm stop,
%   in string form.
%
%   ERR is a struct with several error measures at each algorithm
%   iteration:
%
%     ERR.rawstress: Raw stress (sigma_r), as defined above.
%
%     ERR.stress1:   Stress-1 (sigma_1).
%
%         sigma_1 = sqrt(sigma_r / sum_{i,j} (d_ij)^2)
%
%     ERR.maxalpha:  Maximum angular displacement of any point in each
%                    iteration.
%
%     ERR.medalpha:  Median angular displacement of all points in each
%                    iteration.
%
%   DOUT is the matrix of Euclidean distances between the points given by
%   U, V.
%
% [...] = smdscale(D, U0, V0, OPT)
%
%   U0, V0 are the user's initial guess for the point parametrization. By
%   default, a configuration of points randomly distributed on a
%   [0,0]x[1,1] square is used.
%
%   OPT is a struct with the stop conditions for the algorithm. The
%   algorithm stops when any of the stop conditions are met. The following
%   conditions are available:
%
%     opt.MaxIter: (default = 20) maximum number of iterations we allow the
%                  optimisation algorithm. Each iteration consists of an
%                  entire sweep of all points on the sphere.
%
%     opt.MaxInc:  inc is the Euclidean distance that each output point is
%                  moved in an iteration. The algorithm will stop if no
%                  point has been moved more than MaxInc.
%
%
% A. Agarwal et al. (2010), "Universal Multi-Dimensional Scaling",
% Proceedings of the 16th ACM SIGKDD International Conference on Knowledge
% Discovery and Data Mining, 53:1149-1158.
%
% J. C. Gower (1968), "Adding a point to vector diagrams in multivariate
% analysis", Biometrika, vol. 55, no. 3, pp. 582–585.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
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
narginchk(1, 4);
nargoutchk(0, 5);

% number of points
N = size(d, 1);
if (size(d, 1) ~= N || size(d, 2) ~= N)
    error('D must be a square matrix')
end

% defaults
if (nargin < 2 || isempty(u))
    % random distribution of points on a square [0, 0] x [1, 1]
    u = rand(1, N);
end
if (nargin < 3 || isempty(v))
    % random distribution of points on a square [0, 0] x [1, 1]
    v = rand(1, N);
end
if (nargin >= 4 && ~isstruct(opt))
    error('OPT must be a struct')
end
if (nargin < 4 || isempty(opt))
    opt.MaxIter = 20;
end

% make sure that u and v are column vectors
u = u(:);
v = v(:);

% check that the initial guess is a valid matrix
if (length(u) ~= N)
    error('U must be a vector with length(D) elements')
end
if (length(v) ~= N)
    error('V must be a vector with length(D) elements')
end

% we consider zero and infinite values in the distance matrix the same, as
% a lack of connection
d(isinf(d)) = 0;

% we computate stress measures only if the user requests them
if (nargout >= 4)
    
    if issparse(d) % local neighbourhoods

        % for local neighbourhoods, most points are not connected, so it's
        % inefficient to compute distance between them. Instead, we keep a
        % list of connected points
        [Icon, Jcon] = find(d ~= 0);
        
        % angular distance between points
        dout = distanceang(lat(Icon), lon(Icon), lat(Jcon), lon(Jcon));
        
        % compute different types of stress
        err.rawstress = norm(dout - d(d ~= 0));
        normdcon = norm(d(d ~= 0)); % no need to repeat this computation
        err.stress1 = sqrt(norm(dout - d(d ~= 0))/normdcon);
    
    else % fully connected graphs
        
        % distance between points
        dout = dmatrix([u v]', [u v]');
        
        % keep only distances between connected points
        dout(d == 0) = 0;
        
        % compute different types of stress
        err.rawstress = norm(dout - d);
        normdcon = norm(d); % no need to repeat this computation
        err.stress1 = sqrt(norm(dout - d)/normdcon);
    
    end
    
    % initialize vector for relocation distance
    err.maxinc = nan;
    err.medinc = nan;
end


% keep track of how much each point gets moved in an iteration
inc = zeros(N, 1);

% iterate points. Each point is rearranged to try to get as close as
% possible to the target distance of the points it's connected to
stopCondition = []; niter = 0;
while isempty(stopCondition)
    
    % an iteration consists in N steps. In each step, we re-adjust one
    % point with respect to only its neighbourhood
    niter = niter + 1;
    for I = 1:N
        
        % save the current point coordinates before we start relocating it
        u0 = u(I);
        v0 = v(I);
        
        % neighbours of current point
        idx = d(:, I) ~= 0;
        
        % relocation of current point
        aux = coords_from_dist_gower([u(idx) v(idx)]', d(idx, I));
        u(I) = aux(1);
        v(I) = aux(2);
        
        % save the distance it has moved
        inc(I) = norm([u0, v0] - [u(I) v(I)]);
        
    end % end loop of N points

    % update stress measures
    if (nargout >= 4)
        
        if issparse(d) % local neighbourhoods
            
            % angular distance between points
            dout = distanceang(lat(Icon), lon(Icon), lat(Jcon), lon(Jcon));
            
            % compute different types of stress
            err.rawstress(end+1) = norm(dout - d(d ~= 0));
            err.stress1(end+1) = sqrt(norm(dout - d(d ~= 0))/normdcon);
            
        else % fully connected graphs
            
            % distance between points
            dout = dmatrix([u v]', [u v]');
            
            % keep only distances between connected points
            dout(d == 0) = 0;
            
            % compute different types of stress
            err.rawstress(end+1) = norm(dout - d);
            err.stress1(end+1) = sqrt(norm(dout - d)/normdcon);
        end
        
        err.maxinc(end+1) = max(abs(inc));
        err.medinc(end+1) = median(abs(inc));
    end
        
    % check stop conditions (more than one may be true)
    if (isfield(opt, 'MaxIter') && (niter >= opt.MaxIter))
        stopCondition{end+1} = 'MaxIter';
    end
    
    if (isfield(opt, 'MaxInc') && (max(abs(inc)) < opt.MaxInc))
        stopCondition{end+1} = 'MaxInc';
    end
        
end

if (nargout >= 4)
    % reformat vector dout as sparse matrix if the input is also a sparse
    % matrix
    if issparse(d)
        dout = sparse(dout);
    end
end

end
