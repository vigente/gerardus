function [xi, tri] = surface_interp(uv, x, interp, ui, vi)
% SURFACE_INTERP  Interpolate an open or closed surface from a parametrized
% scattered set of points.
%
%   This function takes a scattered data set (X, Y, Z) in 3D and their 2D
%   parametrization U, V (on a plane or sphere) and computes three
%   interpolants to map (U, V) -> (X, Y, Z):
%
%     X = F1(U, V)
%     Y = F2(U, V)  => (X, Y, Z) = F(U, V)
%     Z = F3(U, V)
%
%   Then, it computes surface coordinates at (UI, VI), i.e.
%
%     XI = F(UI, VI)
%
%   Note: The (U, V) parametrization can be computed with surface_param().
%
% XI = surface_interp(UV, X, INTERP, UI, VI)
%
%   UV is a 2-column matrix with the parameterisation of X, (U, V)->X. In
%   planar parameterisations, UV has units of meters. In spherical
%   parameterisations, UV=[LATITUDE, LONGITUDE]=[ELEVATION, AZIMUTH], in
%   units of radians.
%
%   X is a 3-column matrix. Each row has the coordinates of a point that
%   belongs to the surface we want to interpolate.
%
%   UV and X must have the same number of rows.
%
%   INTERP is a string or struct to select the interpolation method. The
%   struct also allows to pass parameters to the interpolation method. 
%
%     Interpolation methods summary (INTERP.type):
%
%     Open surfaces:
%                 'tsi':      Matlab's TriScatteredInterp() function.
%                 'tps':      Thin-plate spline.
%                 'gridfit':  John D'Errico's [3] approximating gridfit()
%                             function.
%                 'mba':      Multilevel B-Spline Approximation Library by
%                             SINTEF ICT [4].
%                 'mbae':     'mba' with extrapolation over whole domain.
%
%     Closed surfaces:
%                 'sphspline': Spherical spline interpolation by Wessel et
%                             al. [5].
%
%   UI, VI are vectors or arrays of the same size with the surface
%   coordinates at which we want to interpolate. They have the same units
%   as UV.
%
%   XI is a 3-column matrix with the 3D coordinates of the interpolated
%   surface points.
%   
%=========================== OPEN SURFACES ================================
%
%     'tsi'  (default): Matlab's TriScatteredInterp() function. Compact
%     support, limited to the convex hull of the scattered points 2D
%     projection on the interpolation domain.
%
%     'tps': Thin-plate spline. Global support.
%
%     'gridfit': John D'Errico's gridfit() function [3] (note:
%     approximation, rather than interpolation). Compact support with
%     extrapolation outside the convex hull.
%
%     'mba': Multilevel B-Spline Approximation Library by SINTEF ICT [4].
%     Compact support, limited to a rectangle that tighly contains the
%     scattered points 2D projection on the interpolation domain.
%
%        INTERP.nlev: Number of levels in the hierarchical construction of
%                     the B-spline (default 7).
%
%     'mbae': Like the 'mba' method, but first a thin-plate spline is used 
%     to extrapolate values on the interpolation domain boundary. Then MBA
%     is used for the local support interpolation of X and the boundary
%     values set by the thin-plate spline.
%
%========================== CLOSED SURFACES ===============================
%
%     'sphspline': Spherical spline interpolation by Wessel et al. [5].
%
%        INTERP.tension: Tension applied to the sphere spline. More tension
%                     makes the spline smoother (default 0). 
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
% Copyright © 2010-2014 University of Oxford
% Version: 0.0.1
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
narginchk(5, 5);
nargoutchk(0, 2);

if (size(x, 2) ~= 3)
    error('X must be a 3-column matrix')
end
if (any(size(ui) ~= size(vi)))
    error('UI and VI must have the same size')
end


% general defaults
if (isempty(interp) || ~isfield(interp, 'type') || isempty(interp.type))
    error('No INTERP.type provided')
end

%% default parameters for each method

switch interp.type
    
    case {'tsi', 'tps', 'gridfit'}
        
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
        
        if (~isfield(interp, 'N') || isempty(interp.N))
            interp.N = 200;
        end
        if (~isfield(interp, 'tension') || isempty(interp.tension))
            interp.tension = 0;
        end
        
    otherwise
        
        error('Interpolation method not implemented')
end

%% testing warning

% this function is a cannibalization of old surface_interpolation(). Only
% interpolation type 'sphspline' has been tested. We give a warning for
% other interpolation methods, to make the user aware of this
if ~strcmp(interp.type, 'sphspline')
    warning('Interpolation method has not been tested with this function yet')
end

%% compute interpolation

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
        
        % convert radians to degrees, because sphsplinet()
        % expects degrees as the input unit. Also, use lat/lon variables to
        % make uv more readable
        lat = uv(:, 1) * (180/pi);
        lon = uv(:, 2) * (180/pi);
        ui = ui * (180/pi);
        vi = vi * (180/pi);
        
        % compute spherical interpolation
        xi = cat(ndims(ui), ...
            sphsplinet(lon, lat, x(:, 1), vi, ui, interp.tension), ...
            sphsplinet(lon, lat, x(:, 2), vi, ui, interp.tension), ...
            sphsplinet(lon, lat, x(:, 3), vi, ui, interp.tension) ...
            );
        
    otherwise
        
        error('Interpolation method not implemented')
end
