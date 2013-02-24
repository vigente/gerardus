function [lat, lon, sphrad, err, stopCondition, dsph, x] = smdscale(d, x, opt)
% SMDSCALE  Multidimensional scaling on a sphere
%
% This function solves a Multidimensional Scaling problem by finding a way
% to distribute a set of points on a sphere so that the matrix formed from
% the geodesic distance between all pairs of points is as close as possible
% (in the Frobenius norm sense) to an input distance matrix.
%
% [LAT, LON, SPHRAD, ERR, STOP, DSPH, X] = smdscale(D)
%
%   D is a square matrix where D(i,j) is the geodesic distance (e.g. in
%   meters) between the i-th and j-th points. If D is full, then D is taken
%   as provided by the user. If D is sparse, then this means that the user
%   has provided distances only in local neighbourhoods of the matrix (i.e.
%   we only know the distance of each point to a few of its neighbours). D
%   is then pre-processed running the Dijkstra shortest-path algorithm to
%   obtain a full matrix with a distance value for every D(i,j).
%
%   LAT, LON are vectors with the latitude and longitude coordinates of the
%   point configuration that is the solution to the problem in radians.
%
%   SPHRAD is the radius of the sphere that contains the LAT, LON
%   projection of the points. 
%
%   To compute the (x,y,z) Euclidean coordinates of the sphere points run
%
%     [x, y, z] = sph2cart(lon, lat, sphrad);
%
%   ERR is a struct with the error measure at each iteration of the
%   optimisation algorithm.
%
%     ERR.full is the Frobenius norm of the matrix (Dfull-DSPH)/2, where
%     Dfull is the full matrix after running Dijkstra on input matrix D
%     (If D was provided as a full matrix, then D==Dfull), and DSPH is
%     explained below.. This error must be monotonically decreasing in each
%     iteration of the algorithm.
%
%     ERR.sparse is the Frobenius norm of the matrix (D(idx)-DSPH(idx)),
%     where idx are the elements in the input sparse matrix. This result is
%     only returned if D was provided as a sparse matrix. This error is not
%     necessarily monotonically decreasing, but it is expected that the
%     general trend is that it decreases as the algorithm runs.
%
%   STOP is a string with the condition that made the algorithm stop.
%
%   DSPH is the matrix of arc length distances between the points given
%   by LAT, LON (units: meters).
%
%   X is a (3, P)-matrix where each column has the coordinates of one of
%   the points that were used to initialize the algorithm. These points can
%   be provided (see below), or else they will be randomly distributed on
%   the sphere.
%
% [LAT, LON, ERR, X] = smdscale(D, X, OPT)
%
%   X is the (3, P)-matrix mentioned above. If it is provided, then X at
%   the output is the same as X at the input.
%
%   OPT is a struct with the stop conditions for the algorithm. For the
%   moment, only two conditions are available:
%
%     opt.fronorm: target value of the Frobenius norm of the full distance
%                  matrix error, err.full (default 0.001).
%
%     opt.frorel:  relative decrease of the Frobenius norm, e.g. 0.1. This
%                  condition is only considered when the norm is
%                  decreasing, not if it increases.
%
%     opt.maxiter: maximum number of iterations we allow the optimisation
%                  algorithm. Each iteration consists of an entire sweep of
%                  all points on the sphere.
%
% This function implements the MDS optimisation idea proposed by Agarwal et
% al. (2010), but avoids using chords or approximating the mean on the
% sphere by the Euclidean mean followed by projection on the sphere (code
% not available, but explained in personal communication). Many thanks to
% A. Agarwal for his comments about his paper.
%
% A. Agarwal et al. (2010) "Universal Multi-Dimensional Scaling",
% Proceedings of the 16th ACM SIGKDD International Conference on Knowledge
% Discovery and Data Mining, 53:1149-1158.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012-2013 University of Oxford
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
narginchk(1, 3);
nargoutchk(0, 7);

% number of points
N = size(d, 1);
if (size(d, 1) ~= N || size(d, 2) ~= N)
    error('D must be a square matrix')
end

% keep a backup of the input distance matrix
dtarget = d;

% compute full matrix of distances from the sparse ajacency/distance matrix
SPARSE = issparse(d);
if SPARSE
    d = dijkstra(sparse(d), 1:N);
end

% defaults
if (nargin < 2 || isempty(x))
    % random distribution of points on the sphere
    lat = (rand(1, N)-.5) * pi;
    lon = rand(1, N) * pi * 2;
    r = ones(1, N);
    [x, y, z] = sph2cart(lon, lat, r);
    x = cat(1, x, y, z);
    clear y z
end
if (nargin < 3 || isempty(opt) || ~isfield(opt, 'fronorm'))
    opt.fronorm = .001;
end

% center points around zero
x = x * (eye(N) - ones(N)/N);

% project 3D Cartesian coordinates of the initial point configuration onto
% the surface of a sphere lon, lat given in radians
[lon, lat, r] = cart2sph(x(1, :), x(2, :), x(3, :));

% we use the mean radius of the points in the configuration as the radius
% of the sphere that best can contain their projections
sphrad = mean(r);

% matrix of arc length distances between the points in the sphere
dsph = arclen_dmatrix(lat, lon, sphrad);

% if the input distance matrix is sparse, it means that our actual target
% distances that we are trying to match are only those in the sparse
% matrix. However, note that the algorithm will monotonically minimise the
% Frobenius norm of the full matrix error. This doesn't imply that the
% minimisation of the Frobenius norm of the sparse matrix error will be
% monotomic too
if SPARSE
    % sparse matrix elements. Because the distance matrix is symmetrical,
    % we only need the upper triangular part of it
    idx = (dtarget ~= 0) & triu(ones(N), 1);
    
    % Frobenius norm of the sparse error matrix (divided by 2). This error
    % is not necessarily monotonically decreasing, although it is expected
    % to decrease in most iterations
    err.sparse = norm(dsph(idx) - dtarget(idx));
end

% Frobenius norm of the full error matrix (divided by 2). This error must
% be monotonically decreasing
err.full = norm(dsph - dtarget,'fro')/2;

% convert distances in the distance matrix from geodesic distances (e.g. in
% meters) to distances on the sphere given in radians. The reference sphere
% has the radius just computed above. Note that despite km2rad() expecting
% distances in km, as long as we give the distance and radius in the same
% units, the result is the same
drad = km2rad(d, sphrad);

% iterate points. Each point is rearranged to try to get as close as
% possible to the target distance of every other point
stopCondition = []; niter = 0;
while isempty(stopCondition)
    niter = niter + 1;
    for I = 1:N
        
        % compute azimuth of the great circles that conects our current
        % point to each other point. The azimuth is just a convenient way
        % of describing the geodesic from the current point to each other
        % point
        %
        % (great circle: intersection of the sphere and a plane which
        % passes through the center point of the sphere; great circles are
        % the geodesics of the sphere -
        % http://en.wikipedia.org/wiki/Great_circle)
        %
        % (azimuth: "angle at which a smooth curve crosses a meridian,
        % taken clockwise from north. The North Pole has an azimuth of 0º
        % from every other point on the globe." -
        % http://www.mathworks.com/help/map/ref/azimuth.html)
        az = azimuth(lat, lon, lat(I), lon(I), 'radians');
        
        % move each point towards current point, along the great circle
        % that connects them
        [lat2, lon2] = reckon(lat, lon, drad(I, :), az, 'radians');
        
        % remove the current point from the result (that won't have moved)
        lat2(I) = [];
        lon2(I) = [];
        
        % the current point's new position will be the mean of where each
        % point wants to go. The mean has to be a special mean computed on
        % the sphere
        [lat(I), lon(I)] = meanm(lat2, lon2, 'radians');
        
        % update the matrix of arc length distances
        dsph(I, :) = distance(lat(I), lon(I), lat, lon, [sphrad 0], 'radians');
        dsph(:, I) = dsph(I, :)';
        
        % compute Frobenius error at this iteration
        if SPARSE
            err.sparse(end+1) = norm(dsph(idx) - dtarget(idx));
        end
        err.full(end+1) = norm(dsph - dtarget,'fro')/2;
        
    end
    
    % check stop conditions (more than one may be true)
    if (isfield(opt, 'maxiter') && (niter >= opt.maxiter))
        stopCondition{end+1} = 'maxiter';
    end
    
    if (isfield(opt, 'frorel') && length(err.full) > 1 ...
            && err.full(end) < err.full(end-1) ...
            && (err.full(end-1)-err.full(end))/err.full(end-1) < opt.frorel)
        stopCondition{end+1} = 'frorel';
    end
    
    if (err.full(end) <= opt.fronorm)
        stopCondition{end+1} = 'fronorm';
    end
        
end

end

% matrix of distances between points on the sphere, given in arc length
% (i.e. units are meters)
function d = arclen_dmatrix(lat, lon, sphrad)

% number of points on the sphere
N = length(lat);

% compute distance matrix row by row
d = zeros(N);
for I = 1:N
    d(I, :) = distance(lat(I), lon(I), lat, lon, [sphrad 0], 'radians');
end

end
