function [lat, lon, err, stopCondition, dsph, sphrad] = smdscale(d, sphrad, lat, lon, opt)
% SMDSCALE  Multidimensional scaling on a sphere
%
% This function solves a Multidimensional Scaling problem by finding a way
% to distribute (embed) a set of points on a sphere so that the matrix
% formed from the geodesic distance between connected pairs of points is as
% close as possible (in the Frobenius norm sense) to an input distance
% matrix.
%
% [LAT, LON, ERR, STOP, DSPH, SPHRAD] = smdscale(D)
%
%   D is a square matrix where D(i,j) is some distance (e.g. in meters)
%   between points i and j. Note that the points themselves are unknown, we
%   only know the distances between them. These distances could be geodesic
%   distances, or approximations thereof, e.g. chord-length distances.
%   Values of D(i,j)=0 or D(i,j)=Inf mean that points i and j are not
%   neighbours, i.e. connected by an edge.
%
%   The optimisation algorithm works by relocating each point i so that its
%   new distance on the sphere to each neighbour j is as close to D(i,j) as
%   possible. Note that non-neighbours (D=0 or D=Inf) are ignored in the
%   optimisation. This allows sparse matrices and local neighbourhood
%   optimisations.
%
%   Global optimisations (i.e. the graph is complete, so every pair of
%   points is connected) can be achieved too by converting a local
%   neighbourhood graph to a complete graph:
%
%     D = dijkstra(sparse(D), 1:length(d));
%
%   LAT, LON are vectors with the latitude and longitude coordinates (in
%   radians) of a set of points on the sphere. The great circle distances
%   between the points approximate D as well as possible in the Frobenius
%   norm sense (ignoring Inf values).
%
%   To compute the (x,y,z) Euclidean coordinates of the sphere points run
%
%     [x, y, z] = sph2cart(LON, LAT, SPHRAD);
%
%   ERR is a vector with the Frobenius norm error measure at each step of
%   the optimisation algorithm (moving each point counts as a step). Note
%   that for local neighbourhood matrices, the Frobenius norm will be
%   smaller because there are more elements that don't contribute towards
%   the error (because they are 0 or Inf).
%
%   STOP is a cell-array with the condition/s that made the algorithm stop,
%   in string form.
%
%   DSPH (units: meters) is the matrix of great circle or geodesic
%   distances between the points given by LAT, LON, and SPHRAD. The great
%   circle distance between points i and j can be computed in Matlab as
%
%     dsph = distance(lat(i), lon(i), lat(j), lon(j), [sphrad 0], 'radians');
%
%   SPHRAD is a scalar with the radius of the sphere. Note that the radius
%   of the sphere is a critical parameter. If the sphere is too small or
%   too large, the algorithm will produce very poor results. By default, it
%   is estimated so that it can accommodate the distance from each point to
%   the most distant point. If D is a local neighbourhood matrix, this
%   requires running a full Dijkstra on D, which may be slow for large
%   matrices. The sphere's radius is estimated from the full D matrix as:
%
%      sphrad = median(max(D)) / pi;
%
% [...] = smdscale(D, SPHRAD, LAT0, LON0, OPT)
%
%   SPHRAD can be provided by the user, in which case it is not estimated
%   internally.
%
%   LAT0, LON0 are the user's initial guess for the latitude and longitude
%   of the point configuration. By default, a configuration of points
%   randomly distributed over the sphere is used.
%
%   OPT is a struct with the stop conditions for the algorithm. The next
%   conditions are available:
%
%     opt.maxiter: (default = 20) maximum number of iterations we allow the
%                  optimisation algorithm. Each iteration consists of an
%                  entire sweep of all points on the sphere.
%
%     opt.fronorm: target value of the Frobenius norm of the distance
%                  matrix error.
%
%     opt.frorel:  relative decrease of the Frobenius norm, e.g. 0.1. This
%                  condition is only considered when the norm is
%                  decreasing, not if it increases.
%
% This function implements the MDS optimisation idea proposed by Agarwal et
% al. (2010), but avoids using chords or approximating the mean on the
% sphere by the Euclidean mean followed by projection on the sphere (code
% not available, but explained in personal communication). It also improves
% on Agarwal et al in that it provides an estimate of the sphere's radius,
% and it can work with a full D matrix or a local neighbourhood D matrix.
% Many thanks to A. Agarwal for his comments about his paper.
%
% A. Agarwal et al. (2010) "Universal Multi-Dimensional Scaling",
% Proceedings of the 16th ACM SIGKDD International Conference on Knowledge
% Discovery and Data Mining, 53:1149-1158.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012-2013 University of Oxford
% Version: 0.4.0
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
narginchk(1, 5);
nargoutchk(0, 6);

% number of points
N = size(d, 1);
if (size(d, 1) ~= N || size(d, 2) ~= N)
    error('D must be a square matrix')
end

% defaults
if (nargin < 2 || isempty(sphrad))
    % is the distance matrix a local neighbourhood matrix?
    if nnz(d~=0 & ~isinf(d)) < numel(d) - length(d)
        % we need the full matrix to estimate the sphere's radius
        dfull = dijkstra(sparse(d), 1:N);
        
        % estimate the size of the sphere so that it can accommodate the
        % distance matrix
        %
        % first, we find the furthest point from each point. The median of
        % the corresponding distances give an estimate of the half
        % circumference of the sphere. The radius = half_circumference/pi
        % because circumference = 2*pi*radius
        sphrad = median(max(dfull)) / pi;
    else
        % we estimate the radius directly from the input distance matrix
        sphrad = median(max(d)) / pi;
    end
    
end
if (nargin < 3 || isempty(lat))
    % random distribution of points on the sphere
    lat = (rand(1, N)-.5) * pi;
end
if (nargin < 4 || isempty(lon))
    % random distribution of points on the sphere
    lon = rand(1, N) * pi * 2;
end
if (nargin < 5 || isempty(opt) || ~isfield(opt, 'maxiter'))
    opt.maxiter = 20;
end

% check that the initial guess is a valid matrix
if (length(lat) ~= N)
    error('LAT must be a vector with length(D) elements')
end
if (length(lon) ~= N)
    error('LON must be a vector with length(D) elements')
end

% make sure that lat and lon vectors are row vectors
lat = lat(:)';
lon = lon(:)';

% matrix of arc length distances between the points in the sphere
dsph = arclen_dmatrix(lat, lon, sphrad);

% when the optimisation algorithm re-positions a point, it only takes into
% accound its distance to the points it's directly connected to. This
% improves convergence to a better error. Here we compute the adjacency
% matrix (we assume that d=0 or d->Inf means that two points are not
% connected)
idx = (d ~= 0) & ~isinf(d);
    
% Frobenius norm of the distance error matrix (Inf values are ignored)
err = norm(dsph(idx) - d(idx));

% convert distances in the distance matrix from geodesic distances (e.g. in
% meters) to distances on the sphere given in radians. The reference sphere
% has the radius just computed above. Note that despite km2rad() expecting
% distances in km, as long as we give the distance and radius in the same
% units (e.g. both in meters), the result is the same
drad = km2rad(d, sphrad);

% iterate points. Each point is rearranged to try to get as close as
% possible to the target distance of the points it's connected to
stopCondition = []; niter = 0;
while isempty(stopCondition)

    % an iteration consists in N steps. In each step, we re-adjust one
    % point with respect to only its neighbourhood
    niter = niter + 1;
    for I = 1:N
        
        % neighbours of current point
        neigh = idx(I, :);
        
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
        az = azimuth(lat(neigh), lon(neigh), lat(I), lon(I), 'radians');
        
        % move each point towards current point, along the great circle
        % that connects them
        [lat2, lon2] = reckon(lat(neigh), lon(neigh), ...
            drad(I, neigh), az, 'radians');
        
        % the current point's new position will be the mean of where each
        % neighbour tells it to go. The mean has to be a special mean
        % computed on the sphere
        [lat(I), lon(I)] = meanm(lat2, lon2, 'radians');
        
        % update the matrix of arc length distances
        dsph(I, neigh) = distance(lat(I), lon(I), ...
            lat(neigh), lon(neigh), [sphrad 0], 'radians');
        dsph(neigh, I) = dsph(I, neigh)';
        
        % compute Frobenius norm of the error at this iteration (of the
        % local neighbourhood distance matrix)
        err(end+1) = norm(dsph(idx) - d(idx));
        
    end
    
    % check stop conditions (more than one may be true)
    if (isfield(opt, 'maxiter') && (niter >= opt.maxiter))
        stopCondition{end+1} = 'maxiter';
    end
    
    if (isfield(opt, 'frorel') && length(err) > 1 ...
            && err(end) < err(end-1) ...
            && (err(end-1)-err(end))/err(end-1) < opt.frorel)
        stopCondition{end+1} = 'frorel';
    end
    
    if (isfield(opt, 'fronorm') && err(end) <= opt.fronorm)
        stopCondition{end+1} = 'fronorm';
    end
        
end

% turn lat and lon into column vectors
lat = lat(:);
lon = lon(:);

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
