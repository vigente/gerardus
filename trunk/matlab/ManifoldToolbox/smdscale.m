function [lat, lon, stopCondition, err, dout, sphrad] = smdscale(d, sphrad, lat, lon, opt)
% SMDSCALE  Local neighbourhood Multidimensional scaling on a sphere
%
% This function solves a Multidimensional Scaling problem by finding a way
% to distribute (embed) a set of points on a sphere so that the matrix
% formed from the geodesic distance between connected pairs of points is as
% close as possible (in the Frobenius norm sense) to the input distance
% matrix.
%
% This function is based on the MDS optimisation idea proposed by Agarwal
% et al. (2010) for the sphere in case the error metric is raw stress
% (sigma_r) of the geodesic distances
%
%    sigma_r = sqrt(sum_{i,j} (d_ij - dout_ij)^2)
%
% No code was available with the paper, which does not provide the
% implementation details for the sphere either. Our implementation uses
% spherical trigonometry formulas to compute geodesic distances, means and
% point displacements over the sphere.
%
% The original paper assumes that the input distance matrix is full (i.e.
% all points are connected to all points). We allow for sparse input
% matrices, that represent local neighbourhoods of connectivity.
%
% [LAT, LON, STOP, ERR, DOUT, SPHRAD] = smdscale(D)
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
%   new distance on the sphere to each neighbour j is as close to D(i,j) as
%   possible. Note that non-neighbours (D=0 or D=Inf) are ignored in the
%   optimisation. This allows sparse matrices and local neighbourhood
%   optimisations.
%
%   LAT, LON give the output parametrization of the unknown points. LAT,
%   LON are vectors with the latitude and longitude coordinates (in
%   radians). Note that LAT, LON is one of infinite solutions, as D is
%   invariant to rotations.
%
%   To compute the (x,y,z) Euclidean coordinates of the sphere points run
%
%     [x, y, z] = sph2cart(LON, LAT, 1);
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
%   DOUT is the matrix of great circle or geodesic distances between the
%   points given by LAT, LON, and SPHRAD. The great circle distance between
%   points i and j can be computed in Matlab as
%
%     dsph = distance(lat(i), lon(i), lat(j), lon(j), sphrad, 'radians');
%
%   SPHRAD is a scalar with the radius of the sphere used in the MDS
%   algorithm. Note that the radius of the sphere is a critical parameter.
%   If the sphere is too small or too large, the algorithm will produce
%   very poor results. If not provided, SPHRAD is estimated so that it can
%   accommodate the distance from each point to the most distant point.
%   Note that this requires running a full Dijkstra, which may be slow for
%   very sparse large matrices. The sphere's radius is estimated from the
%   Dijkstra distance matrix DD as:
%
%      sphrad = median(max(DD)) / pi;
%
%   where DD=dijkstra(sparse(D), 1:N);
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
%   OPT is a struct with the stop conditions for the algorithm. The
%   algorithm stops when any of the stop conditions are met. The following
%   conditions are available:
%
%     opt.MaxIter: (default = 20) maximum number of iterations we allow the
%                  optimisation algorithm. Each iteration consists of an
%                  entire sweep of all points on the sphere.
%
%     opt.MaxAlpha: alpha is the angular distance (in radians) that each
%                   output point is moved in an iteration. The algorithm
%                   will stop if no point has been moved more than MaxAlpha
%                   radians.
%
%
% A. Agarwal et al. (2010) "Universal Multi-Dimensional Scaling",
% Proceedings of the 16th ACM SIGKDD International Conference on Knowledge
% Discovery and Data Mining, 53:1149-1158.
%
% J. P. Snyder (1987), "Map Projections: A Working Manual",  US Geological
% Survey Professional Paper 1395, US Government Printing Office,
% Washington, DC, 1987.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012-2013 University of Oxford
% Version: 0.5.1
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
    % we need the full matrix to estimate the sphere's radius. The user
    % could have provided a full distance matrix, but with lots of zeros
    % because not all elements are connected
    dfull = dijkstra(sparse(d), 1:N);
    
    % estimate the size of the sphere so that it can accommodate the
    % distance matrix
    %
    % first, we find the furthest point from each point. The median of
    % the corresponding distances give an estimate of the half
    % circumference of the sphere. The radius = half_circumference/pi
    % because circumference = 2*pi*radius
    sphrad = median(max(dfull)) / pi;
    clear dfull
end
if (nargin < 3 || isempty(lat))
    % random distribution of points on the sphere
    lat = (rand(1, N)-.5) * pi;
end
if (nargin < 4 || isempty(lon))
    % random distribution of points on the sphere
    lon = rand(1, N) * pi * 2;
end
if (nargin >= 5 && ~isstruct(opt))
    error('OPT must be a struct')
end
if (nargin < 5 || isempty(opt))
    opt.MaxIter = 20;
end

% check that the initial guess is a valid matrix
if (length(lat) ~= N)
    error('LAT must be a vector with length(D) elements')
end
if (length(lon) ~= N)
    error('LON must be a vector with length(D) elements')
end

% make sure that lat and lon vectors are column vectors
lat = lat(:);
lon = lon(:);

% we consider zero and infinite values in the distance matrix the same, as
% a lack of connection
d(isinf(d)) = 0;

% convert distances to arc distances. This way, all computations are
% independent of the radius
d = d / sphrad;

% we computate stress measures only if the user requests them
%
% note that distances will be in radians
if (nargout >= 4)

    if issparse(d) % local neighbourhoods

        % for local neighbourhoods, most points are not connected, so it's
        % inefficient to compute distance between them. Instead, we keep a
        % list of connected points
        [Icon, Jcon] = find(d ~= 0);
        
        % angular distance between points
        dout = distanceang(lat(Icon), lon(Icon), lat(Jcon), lon(Jcon));
        
        % compute different types of stress
        err.rawstress = sphrad * norm(dout - d(d ~= 0));
        normdcon = norm(d(d ~= 0)); % no need to repeat this computation
        err.stress1 = sqrt(norm(dout - d(d ~= 0))/normdcon);
    
    else % fully connected graphs
        
        % angular distance between points
        dout = dmatrix_ang(lat, lon);
        
        % keep only distances between connected points
        dout(d == 0) = 0;
        
        % compute different types of stress
        err.rawstress = sphrad * norm(dout - d);
        normdcon = norm(d); % no need to repeat this computation
        err.stress1 = sqrt(norm(dout - d)/normdcon);
    
    end
    
    
    % initialize vector for relocation angular distance
    err.maxalpha = nan;
    err.medalpha = nan;
end


% keep track of how much each point gets moved in an iteration
alpha = zeros(N, 1);

% iterate points. Each point is rearranged to try to get as close as
% possible to the target distance of the points it's connected to
stopCondition = []; niter = 0;
while isempty(stopCondition)
    
    % an iteration consists in N steps. In each step, we re-adjust one
    % point with respect to only its neighbourhood
    niter = niter + 1;
    for I = 1:N
        
        % save the current point coordinates before we start relocating it
        lat0 = lat(I);
        lon0 = lon(I);
        
        while (1) % relocation of current point
            
            % position of current point (target of the displacements)
            lat2 = lat(I);
            lon2 = lon(I);
            
            % neighbours of current point. These indices are obviously
            % necessary in the case of a sparse matrix. In the case of a
            % full matrix, it is more difficult to justify. First, this
            % prevents each point from being its own neighbour. Second, the
            % user may have provided a full matrix, but missing some
            % connections
            neigh = d(:, I) ~= 0;
            
            % position of neighbours (origin of the displacements)
            lat1 = lat(neigh);
            lon1 = lon(neigh);
            
            % compute azimuth of the great circles that conects neighbours
            % to current point. The azimuth is just a compact way of
            % describing a great circle with just a number. Formula from
            % Snyder (1987), p. 31
            %
            %
            % (great circle: intersection of the sphere and a plane which
            % passes through the center point of the sphere; great circles
            % are the geodesics of the sphere -
            % http://en.wikipedia.org/wiki/Great_circle)
            %
            % (azimuth: "angle at which a smooth curve crosses a meridian,
            % taken clockwise from north. The North Pole has an azimuth of
            % 0º from every other point on the globe." -
            % http://www.mathworks.com/help/map/ref/azimuth.html)
            az = atan2(...
                cos(lat2) .* sin(lon2-lon1),...
                cos(lat1) .* sin(lat2) - sin(lat1) .* cos(lat2) .* cos(lon2-lon1));
            
            % azimuths are undefined at the poles. Following Matlab's
            % convention, we make them zero at the north pole and pi at the
            % south pole
            az(lat1 <= -pi/2) = 0;
            az(lat2 >=  pi/2) = 0;
            az(lat2 <= -pi/2) = pi;
            az(lat1 >=  pi/2) = pi;

%             % DEBUG: using instead Matlab's azimuth function (slower)
%             az = azimuth(lat1, lon1, lat2, lon2, 'radians');
            
            % move each neighbour towards current point, along the great
            % circle that connects them, according to the target distance
            [lat1, lon1] = reckongc(lat1, lon1, az, full(d(neigh, I)));
            
%             % DEBUG: using instead Matlab's reckon function (slower)
%             [lat1, lon1] = reckon(lat1, lon1, d(neigh, I), az, 'radians');
            
            % each neighbour has been moved to a new location. We now move
            % the current point to the mean of the neighbours new locations
            [x, y, z] = sph2cart(lon1, lat1, 1);
            [lon(I), lat(I)] = cart2sph(sum(x), sum(y), sum(z));
            
%             % DEBUG: using instead Matlab's meanm (slower)
%             [lat(I), lon(I)] = meanm(lat1, lon1, 'radians');
            
            % if the current point relocation has stabilized, then we move
            % on to the next point (0.0873 rad = 5 degrees)
            if (distanceang(lat2, lon2, lat(I), lon(I)) < 0.0873)
                break
            end
            
        end % while loop for relocation of current point
        
        % once we have stopped moving the current point around, we save the
        % total angular distance it has moved
        alpha(I) = distanceang(lat0, lon0, lat(I), lon(I));
        
    end % end loop of N points
    
    % note that distances will be in radians
    if (nargout >= 4)
        
        if issparse(d) % local neighbourhoods
            
            % angular distance between points
            dout = distanceang(lat(Icon), lon(Icon), lat(Jcon), lon(Jcon));
            
            % compute different types of stress
            err.rawstress(end+1) = sphrad * norm(dout - d(d ~= 0));
            err.stress1(end+1) = sqrt(norm(dout - d(d ~= 0))/normdcon);
            
        else % fully connected graphs
            
            % angular distance between points
            dout = dmatrix_ang(lat, lon);
            
            % keep only distances between connected points
            dout(d == 0) = 0;
            
            % compute different types of stress
            err.rawstress(end+1) = sphrad * norm(dout - d);
            err.stress1(end+1) = sqrt(norm(dout - d)/normdcon);
        end
        
        err.maxalpha(end+1) = max(abs(alpha));
        err.medalpha(end+1) = median(abs(alpha));
    end
    
    % check 
        
    % check stop conditions (more than one may be true)
    if (isfield(opt, 'MaxIter') && (niter >= opt.MaxIter))
        stopCondition{end+1} = 'MaxIter';
    end
    
    if (isfield(opt, 'MaxAlpha') && (max(abs(alpha)) < opt.MaxAlpha))
        stopCondition{end+1} = 'MaxAlpha';
    end
        
end

% turn lat and lon into column vectors
lat = lat(:);
lon = lon(:);

if (nargout >= 4)
    % convert output distances from angular to geodesic distances
    dout = dout * sphrad;
    
    % reformat vector dout as sparse matrix if the input is also a sparse
    % matrix
    if issparse(d)
        dout = sparse(dout);
    end
end

end

% angular distances between points on the sphere
function alpha = distanceang(lat1, lon1, lat2, lon2)

% haversine formula with the modified formulation discussed by Bob
% Chamberlain for better numerical accuracy over large distances
a = sin((lat2-lat1)/2).^2 ...
    + cos(lat2).*cos(lat1).*sin((lon2-lon1)/2).^2;

alpha = 2*atan2(sqrt(a), sqrt(1-a));

end

% angular distances matrix between points on the sphere
function dalpha = dmatrix_ang(lat, lon)

% convert vectors to matrix form
lat = repmat(lat, 1, size(lat, 1));
lon = repmat(lon, 1, size(lon, 1));

dalpha = distanceang(lat, lon, lat', lon');

end

% Matlab has a function reckon() to compute the destination if we move from
% a point on the sphere along a certain great circle a certain distance.
% This function is a wrapper of an internal private function called
% greatcirclefwd(). We can accelerate computations a bit if we simply avoid
% using reckon() and reimplement the functionality of greatcirclefwd() from
% Snyder (1987), p. 31.
%
% Here lat, long are latitude, longitude; az is the great circle's azimuth;
% and drad is the angular distance we want to travel. All distances and
% coordinates in radians
function [lat2, lon2] = reckongc(lat, lon, az, drad)
% On a sphere of radius A, compute points on a great circle at specified
% azimuths and ranges.  PHI, LAMBDA, PHI0, LAMBDA0, and AZ are angles in
% radians, and RNG is a distance having the same units as R.

% points near the pole may have bad estimations of azimuth. Make sure they
% have valid values
epsilon = 10*epsm('radians');    % set tolerance
az(lat >= pi/2-epsilon) = pi;    % starting at north pole
az(lat <= epsilon-pi/2) = 0;     % starting at south pole

% target latitude
lat2 = asin(sin(lat) .* cos(drad) + cos(lat) .* sin(drad) .* cos(az));

% target longitude
lon2 = lon + atan2(sin(drad) .* sin(az), ...
    cos(lat) .* cos(drad) - sin(lat) .* sin(drad) .* cos(az));
end
