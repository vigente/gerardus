function [lat, lon, err, stopCondition, x] = smdscale(d0, x, opt)
% SMDSCALE  Multidimensional scaling on the surface of a sphere
%
% This function solves a Multidimensional Scaling problem by finding a way
% to distribute a set of points on the surface of a sphere so that their
% geodesic distance matrix is close in the Frobenius norm sense to a
% distance matrix provided at the input.
%
% [LAT, LON, ERR, STOP, X] = smdscale(D0)
%
%   D0 is a square matrix where D0(i, j) is the geodesic distance measured
%   in degrees of arc between the i-th and j-th points.
%
%   LAT, LON are vectors with the latitude, longitude coordinates of the
%   point configuration that is the solution to the problem. The matrix of
%   geodesic distances between the points can be computed as
%
%     d = gdmatrix(LAT, LON);
%
%   ERR is a vector of the Frobenius norm of the distance matrix error,
%   i.e. norm(d-D0, 'fro') at each iteration of the algorithm.
%
%   STOP is a string with the condition that made the algorithm stop.
%
%   X is a (3, P)-matrix where each column has the coordinates of one of
%   the points that were used to initialize the algorithm. These points can
%   be provided (see below), or else they are randomly distributed on the
%   sphere.
%
% [LAT, LON, ERR, X] = smdscale(D0, X, OPT)
%
%   X is the (3, P)-matrix mentioned above. If it is provided, then X at
%   the output is the same as X at the input.
%
%   OPT is a struct with the stop conditions for the algorithm. For the
%   moment, only two conditions are available:
%
%     opt.fronorm: value of the Frobenius norm of the distance matrix error
%                  (default 0.001).
%
%     opt.frorel:  relative decrease of the Frobenius norm, e.g. 0.1. This
%                  condition is only considered when the norm is
%                  decreasing, not if it increases.
%
% Previous work on computing MDS on the sphere can be found in e.g.
%
% T.F. Cox & M.A.A. Cox (1991) "Multidimensional scaling on a sphere",
% Communications in Statistics - Theory and Methods, 20(9):2943-2953.
%
% A. Agarwal et al. (2010) "Universal Multi-Dimensional Scaling",
% Proceedings of the 16th ACM SIGKDD International Conference on Knowledge
% Discovery and Data Mining, 53:1149-1158.

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

% check arguments
error(nargchk(1, 3, nargin, 'struct'));
error(nargoutchk(0, 5, nargout, 'struct'));

% number of points
N = size(d0, 1);
if (size(d0, 1) ~= N || size(d0, 2) ~= N)
    error('D0 must be a square matrix')
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

% project 3D Cartesian coordinates onto the surface of a sphere
% lon, lat given in radians
[lon, lat] = cart2sph(x(1, :), x(2, :), x(3, :));

% convert sphere coordinates units from radians to degrees
lon = lon / pi * 180;
lat = lat / pi * 180;

% distance matrix computed from the current configuration of points
d = gdmatrix(lat, lon);

% initial error
err = norm(d - d0, 'fro');

% iterate points. For each point, we rearrange the rest of the points to
% approach the target distance matrix
stopCondition = [];
while isempty(stopCondition)
    for I = 1:N
        
        % compute azimuth of the great circle that conects our current point to
        % each other point
        az = azimuth(lat(I), lon(I), lat, lon);
        
        % move each point towards/away from current point, so that they are at
        % the target distance, moving them along the geodesics
        [lat, lon] = reckon(lat(I), lon(I), d0(I, :), az);
        
        % distance matrix computed from the current configuration of points
        d = gdmatrix(lat, lon);
        
        % compute the Frobenius norm of the difference between the current and
        % target distance matrices
        err(end+1) = norm(d - d0, 'fro');
        
        % check stop conditions
        if (isfield(opt, 'frorel') && length(err) > 1 ...
                && err(end) < err(end-1) ...
                && (err(end-1)-err(end))/err(end-1) < opt.frorel)
            stopCondition = 'frorel';
            break
        elseif (err(end) <= opt.fronorm)
            stopCondition = 'fronorm';
            break
        end
        
    end
end
