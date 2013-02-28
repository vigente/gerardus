function [xi, uv, x, ui, vi] = surface_interpolation(x, param, interp)
% SURFACE_INTERPOLATION  Interpolate a surface from a scattered set of
% points
%
%   This function takes a scattered data set X in 3D and finds a 2D
%   parameterisation U, V (on a plane or sphere). It then computes three
%   interpolants to map (U, V) -> X: X = F(U, V).
%
%   Finally, it creates a regular grid UI, VI and computes the
%   corresponding interpolated XI, XI = F(UI, VI). This allow to visualise
%   the surface resulting from interpolating the points in X.
%
% [XI, UV, X2, UI, VI] = surface_interpolation(X)
%
%   X is a 3-column matrix. Each row has the coordinates of a point that
%   belongs to the surface we want to interpolate.
%
%   XI is a (:,:,3)-volume with the interpolation result on a regular grid
%   automatically created by this function. Each (i,j,:)-vector contains
%   the (x,y,z) Euclidean coordinates of the interpolated point that
%   corresponds to the parameterisation UI(i,j), VI(i,j). To plot the
%   interpolated surface run
%
%     surf(xi(:, :, 1), xi(:, :, 2), xi(:, :, 3))
%
%   UV is a 2-column matrix with the parameterisation of X, (U, V)->X. In
%   planar parameterisations, UV has units of meters. In spherical
%   parameterisations, UV=[LON, LAT], in units of radians.
%
%   X2: Some methods (e.g. 'mbae') add extra points to X before computing
%   the interpolants. X2 contains the points actually used for
%   interpolation.
%
%   UI, VI are the grid for the interpolated domain. UI, VI are rectangular
%   matrices. Their size depends on the sampling rate of the domain
%   (resolution) and the domain's size. In spherical parameterisations,
%   UI=LON and VI=LAT, in units of radians.
%
% ... = surface_interpolation(X, PARAM, INTERP)
%
%   PARAM is a struct that describes the method used to parametrise the
%   surface and the set of points X, and its parameters. PARAM needs to
%   have at least a field PARAM.type that describes the parameterisation
%   method:
%
%     'xy' (default): No change, the X coordinates are kept the same.
%
%        No extra parameters.
%
%     'pca': The points in X are rotated according to their eigenvectors to
%     make the dominant plane of the points X as horizontal as possible
%     before interpolating.
%
%        No extra parameters.
%
%     'isomap': Use the Isomap method by [1] to "unfold" the curved surface
%     defined by X before interpolating. (This option requires function
%     IsomapII).
%
%        PARAM.size: Neighbourhood size (maximum number of closest
%                    neighours allowed for each point). By default,
%                    PARAM.size=4. Small neighbourhoods allow to capture
%                    higher curvature in the surface, but too small
%                    neighbourhoods create several connected components and
%                    will result in an error.
%
%        PARAM.d:    Distance matrix between points. If this matrix is not
%                    provided, it is computed internally as the Euclidean
%                    matrix of distances between every pair of points.
%                    Providing the matrix allows the user to define
%                    specific neighbourhoods and topologies.
%
%     'sphisomap': Our extension to the Isomap method so that points are
%     parameterised on a sphere rather than on a plane.
%
%        PARAM.dtype: Type of distance matrix to optimise: 
%
%                     'none': (default) Points are simply projected onto
%                     the sphere along the radius to the centroid. No MDS
%                     optimisation.
%
%                     'full': The full distance matrix is optimised, after
%                     applying Dijkstra's algorithm to obtain shortest
%                     distances between every pair of points.
%
%                     'sparse': Optimisation is only run on the local
%                     neighbourhoods defined by the distance matrix.
%
%                     'full+sparse': First 'full', then 'sparse'.
%
%        PARAM.*: Any other fields are passed to function
%                 smdscale(...,OPT=PARAM). This allows to pass e.g.
%                 stopping conditions to the spherical MDS algorithm. See
%                 help smdscale for all options.
%
%   INTERP is a struct that describes the interpolation method, and its
%   parameters:
%
%     'tsi'  (default): Matlab's TriScatteredInterp() function. Compact
%     support, limited to the convex hull of the scattered points 2D
%     projection on the interpolation domain.
%
%        INTERP.res:  Resolution of the surface grid (default, [1 1]).
%
%        INTERP.klim: Scalar factor for the extension of the interpolation
%                     domain beyond the tight box that contains UV 
%                     (default 1).
%
%     'tps': Thin-plate spline. Global support.
%
%        See parameters above.
%
%     'gridfit': John D'Errico's gridfit() function [3] (note:
%     approximation, rather than interpolation). Compact support with
%     extrapolation outside the convex hull.
%
%        See parameters above.
%
%     'mba': Multilevel B-Spline Approximation Library by SINTEF ICT [4].
%     Compact support, limited to a rectangle that tighly contains the
%     scattered points 2D projection on the interpolation domain.
%
%        INTERP.res:  Resolution of the surface grid (default, [1 1]).
%
%        INTERP.klim: Scalar factor for the extension of the interpolation
%                     domain beyond the tight box that contains UV 
%                     (default 1).
%
%        INTERP.nlev: Number of levels in the hierarchical construction of
%                     the B-spline (default 7).
%
%     'mbae': Like the 'mba' method, but first a thin-plate spline is used 
%     to extrapolate values on the interpolation domain boundary. Then MBA
%     is used for the local support interpolation of X and the boundary
%     values set by the thin-plate spline.
%
%        See parameters above.
%
%     'sphspline': Sphere spline interpolation [5].
%
%        INTERP.res:  Resolution of the surface grid given as arc increment
%                     (in radians) when sampling the sphere 
%                     (default 2*pi/20).
%
%        INTERP.tension: Tension applied to the sphere spline. More tension
%                     makes the spline smoother (default 0). 
%
%
%
% [1] J.B. Tenenbaum, V. de Silva and J.C. Langford, "A Global Geometric
% Framework for Nonlinear Dimensionality Reduction", Science 290(5500):
% 2319-2323, 2000.
%
% [2] Isomap Homepage, http://isomap.stanford.edu/
%
% [3] Surface Fitting using gridfit by John D'Errico 11 Nov 2005 (Updated
% 29 Jul 2010). Code covered by the BSD License
% http://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
%
% [4] MBA - Multilevel B-Spline Approximation Library
% http://www.sintef.no/Projectweb/Geometry-Toolkits/MBA/
%
% [5] Wessel, P., and J. M. Becker, 2008, Gridding of spherical data using
% a Green's function for splines in tension, Geophys. J. Int., 174, 21-28,
% doi:10.1111/j.1365-246X.2008.03829.x
% http://www.soest.hawaii.edu/wessel/sphspline/

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010-2013 University of Oxford
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
narginchk(1, 3);
nargoutchk(0, 5);

if (size(x, 2) ~= 3)
    error('X must be a 3-column matrix')
end

% general defaults
if (nargin < 2 || isempty(param) || ~isfield(param, 'type') ...
        || isempty(param.type))
    param.type = 'xy';
end
if (nargin < 3 || isempty(interp) || ~isfield(interp, 'type') ...
        || isempty(interp.type))
    interp.type = 'tsi';
end

%% Parameterisation: map the 3D points (x,y,z) to a 2D domain (u,v)

% method-specific defaults
switch param.type
    
    case {'xy', 'pca'}
        
        % nothing to do here
        
    case 'isomap'
        
        % if distance matrix is not provided by the user, it is computed
        % internally as the full Euclidean distance matrix. Note that we
        % can still define local neighbourdhoods with the size parameter
        if (~isfield(param, 'd') || isempty(param.d))
            
            param.d = dmatrix(x', x', 'euclidean');
            
        end
        
        % size of the local neighbourhood (number of nearest neighbours
        % each point is connected to)
        if (~isfield(param, 'size') || isempty(param.size))
            
            param.size = 4;
            
        end
        
    case 'sphisomap'
        
        if (~isfield(param, 'maxiter') || isempty(param.maxiter))
            
            param.maxiter = 50;
            
        end
        if (~isfield(param, 'dtype') || isempty(param.dtype))
            
            param.dtype = 'full';
            
        end
        
        
    otherwise
        
        error('Parametrization method not implemented')
        
end

% obtain a 2D parameterisation for the input 3D points
switch param.type
    
    case 'xy'
        
        % (u,v) is simply (x,y)
        uv = x(:, 1:2);
        
    case 'pca'
        
        % rotate valve points to make the valve surface as horizontal as
        % possible
        m = mean(x, 1);
        uv = x - m(ones(1, size(x, 1)), :);
        eigv = pts_pca(uv');
        uv = uv * eigv;
        uv = uv(:, 1:2);
        
    case 'isomap'

        % compute 2-d projection of the 3-d data
        options.dims = 2;
        options.display = 0;
        options.overlay = 0;
        options.verbose = 0;
        options.dijkstra = 1;
        uv = IsomapII(param.d, 'k', param.size, options);
        uv = uv.coords{1}';
        if (size(uv, 1) ~= size(x, 1))
            error('The neighbourhood size (param.size) is too small to connect all points')
        end
        
    case 'sphisomap'
        
        % number of points to embed
        N = size(x, 1);
        
        % initial guess for the sphere embedding by simply projecting each
        % point along the radius to the centroid
        [lat, lon, sphrad] = proj_on_sphere(x);
        
        if (strcmp(param.dtype, 'full') || strcmp(param.dtype, 'full+sparse'))
            
            % embbed the point set on the sphere using the initial guess, using
            % a full distance matrix as the input to Multidimensional Scaling
            % (MDS)
            [lat, lon] = ...
                smdscale(dijkstra(sparse(param.d), 1:N), ...
                sphrad, lat, lon, param);
            
        end
        
        if (strcmp(param.dtype, 'sparse') || strcmp(param.dtype, 'full+sparse'))

            % use the solution to the full distance matrix to initialise MDS on
            % the sparse distance matrix
            [lat, lon] = ...
                smdscale(param.d, sphrad, lat, lon, param);
            
        end
        
        % create output with the parameterisation values for each point
        % lon = em(1, :)
        % lat = em(2, :)
        uv = [lon lat];
        
    otherwise
        
        error('Parametrization method not implemented')
end

%% Interpolation: compute domain

% method-specific defaults
switch interp.type
    
    case {'tps', 'tsi', 'gridfit'}
        
        if (~isfield(interp, 'res') || isempty(interp.res))
            interp.res = [1 1];
        end
        if (~isfield(interp, 'klim') || isempty(interp.klim))
            interp.klim = 1;
        end
        
    case {'mba', 'mbae'}
        
        if (~isfield(interp, 'res') || isempty(interp.res))
            interp.res = [1 1];
        end
        if (~isfield(interp, 'klim') || isempty(interp.klim))
            interp.klim = 1;
        end
        if (~isfield(interp, 'nlev') || isempty(interp.nlev))
            interp.nlev = 7;
        end
        
    case 'sphspline'
        
        if (~isfield(interp, 'res') || isempty(interp.res))
            interp.res = 2*pi/20; % arc increment (in radians) when sampling the sphere
        end
        if (~isfield(interp, 'tension') || isempty(interp.tension))
            interp.tension = 0;
        end
        
    otherwise
        
        error('Interpolation method not implemented')
end

switch interp.type
    
    case {'tps', 'tsi', 'gridfit', 'mba', 'mbae'}
        
        % find box that contains embedded coordinates
        uvmin = min(uv, [], 1);
        uvmax = max(uv, [], 1);
        
        % box size and centroid
        delta = uvmax - uvmin;
        boxm = mean([uvmax; uvmin], 1);
        
        % extend the box
        uvmin = boxm - delta/2*interp.klim;
        uvmax = boxm + delta/2*interp.klim;
        
        % generate grid for the embedding box
        [vi, ui] = ndgrid(...
            uvmin(2):interp.res(1):uvmax(2), ...
            uvmin(1):interp.res(2):uvmax(1));
        
    case 'sphspline'
        
        % number of sampling points in each direction on the sphere to
        % achieve radian increment requested by the user
        Nsph = 2*pi / interp.res;
        
        % regular sampling of a sphere
        [xisph, yisph, zisph] = sphere(Nsph);
        
        % convert sphere points to lat and lon coordinates (in radians)
        % lon = ui
        % lat = vi
        [ui, vi] = cart2sph(xisph, yisph, zisph);
        
    otherwise
        
        error('Interpolation method not implemented')
end

%% Interpolation: compute surface

% interpolate
switch interp.type
    
    case 'tps' % thin-plate spline
        
        xi = pts_tps_map(uv, x, [ui(:) vi(:)]);
        xi = reshape(xi, [size(ui) 3]);

    case 'tsi' % Matlab's TriScatteredInterp
        
        fx = TriScatteredInterp(uv(:, 1), uv(:, 2), x(:, 1), 'natural');
        fy = TriScatteredInterp(uv(:, 1), uv(:, 2), x(:, 2), 'natural');
        fz = TriScatteredInterp(uv(:, 1), uv(:, 2), x(:, 3), 'natural');
        fx = fx(ui, vi);
        fy = fy(ui, vi);
        fz = fz(ui, vi);
        xi = cat(3, fx, fy, fz);
        
    case 'gridfit'
        
        fx = gridfit(uv(:, 1), uv(:, 2), x(:, 1), ...
            uvmin(1):interp.res(2):uvmax(1), ...
            uvmin(2):interp.res(1):uvmax(2), ...
            'tilesize', 150);
        fy = gridfit(uv(:, 1), uv(:, 2), x(:, 2), ...
            uvmin(1):interp.res(2):uvmax(1), ...
            uvmin(2):interp.res(1):uvmax(2), ...
            'tilesize', 150);
        fz = gridfit(uv(:, 1), uv(:, 2), x(:, 3), ...
            uvmin(1):interp.res(2):uvmax(1), ...
            uvmin(2):interp.res(1):uvmax(2), ...
            'tilesize', 150);
        xi = cat(3, fx, fy, fz);
        
    case 'mba' % Multilevel B-Spline Approximation Library
        
        fx = mba_surface_interpolation(uv(:, 1), uv(:, 2), x(:, 1), ...
            ui(:), vi(:), interp.nlev);
        fy = mba_surface_interpolation(uv(:, 1), uv(:, 2), x(:, 2), ...
            ui(:), vi(:), interp.nlev);
        fz = mba_surface_interpolation(uv(:, 1), uv(:, 2), x(:, 3), ...
            ui(:), vi(:), interp.nlev);
        xi = cat(3, ...
            reshape(fx, size(ui)), ...
            reshape(fy, size(ui)), ...
            reshape(fz, size(ui)));
        
    case 'mbae' % MBA extrapolated using a TPS
        
        % use TPS to interpolate the points, but we are only interested in
        % the edges and corners of the interpolation, where the TPS is
        % extrapolating linearly; also, we decimate the points in the edges
        uv2 = [ui(1, 1:end)' vi(1, 1:end)' ; ... % top edge
            ui(2:end, end) vi(2:end, end) ; ... % right edge
            ui(end, 1:end-1)' vi(end, 1:end-1)' ; ... % bottom edge
            ui(2:end-1, 1) vi(2:end-1, 1) ; ... % left edge
            ];
        xi = pts_tps_map(uv, x, uv2);
        x = [x; xi(1:10:end, :)]; % don't keep all points on the boundary
        % local support interpolation with boundary conditions provided by
        % the TPS
        uv = [uv; uv2(1:10:end, :)];
        
        fx = mba_surface_interpolation(uv(:, 1), uv(:, 2), x(:, 1), ...
            ui(:), vi(:), interp.nlev);
        fy = mba_surface_interpolation(uv(:, 1), uv(:, 2), x(:, 2), ...
            ui(:), vi(:), interp.nlev);
        fz = mba_surface_interpolation(uv(:, 1), uv(:, 2), x(:, 3), ...
            ui(:), vi(:), interp.nlev);
        xi = cat(3, ...
            reshape(fx, size(ui)), ...
            reshape(fy, size(ui)), ...
            reshape(fz, size(ui)));
        
    case 'sphspline'
        
        % constant to convert radians to degrees, because sphsplinet()
        % expects degrees as the input unit
        aux = 180/pi;
        
        % compute spherical interpolation
        xi = cat(3, ...
            sphsplinet(lon*aux, lat*aux, x(:, 1), ui*aux, vi*aux, interp.tension), ...
            sphsplinet(lon*aux, lat*aux, x(:, 2), ui*aux, vi*aux, interp.tension), ...
            sphsplinet(lon*aux, lat*aux, x(:, 3), ui*aux, vi*aux, interp.tension) ...
            );
        
    otherwise
        
        error('Interpolation method not implemented')
end
