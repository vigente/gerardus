function [xi, em, gx, gy] = surface_interpolation(x, PARAM, INTERP, res, KLIM)
% SURFACE_INTERPOLATION  Interpolate a surface from a scattered set of points
%
% [XI, EM, GX, GY] = SURFACE_INTERPOLATION(X)
%
%   X is a 3-row matrix. Each column has the coordinates of a point that
%   belongs to the surface we want to interpolate.
%
%   XI is a 3-row matrix with the coordinates of the interpolated points.
%
%   EM is a 2-row matrix with the coordinates of the X points projected
%   onto the interpolation domain. Note that the interpolation domain may
%   change between methods, so don't expect results to be aligned if you do
%   e.g. plot3(gx(:), gy(:), y(:, 3), 'o').
%
%   GX, GY are the grid for the box that contains EM.
%
% ... = SURFACE_INTERPOLATION(X, PARAM, INTERP, RES, KLIM)
%
%   PARAM is a string with the method used to parametrize the surface and
%   X:
%
%     'xy' (default): No change, the X coordinates are kept the same.
%
%     'pca': X points are rotated according to their eigenvectors to make
%     the dominant plane of the points X as horizontal as possible before
%     interpolating.
%
%     'isomap': Use the Isomap method by [1] to "unfold" the curved surface
%     defined by X before interpolating. (This option requires function
%     IsomapII).
%
%   INTERP is a string with the interpolation method:
%
%      'tps' (default): Thin-plate spline. Global support.
%
%      'tsi': Matlab's TriScatteredInterp() function. Local support,
%      limited to the convex hull of the scattered points 2D projection on
%      the interpolation domain.
%
%      'gridfit': John D'Errico's gridfit() function [3] (note:
%      approximation, rather than interpolation). Local support with
%      extrapolation outside the convex hull.
%
%      'mba': Multilevel B-Spline Approximation Library by SINTEF ICT [4].
%      Local support, limited to a rectangle that tighly contains the
%      scattered points 2D projection on the interpolation domain.
%
%      'mbae': Like the 'mba' method, but first a thin-plate spline is used 
%      to extrapolate values on the interpolation domain boundary. Then MBA
%      is used for the local support interpolation of X and the boundary
%      values set by the thin-plate spline.
%
%   RES is a 2-vector with the grid spacing in the x- and y-directions. By
%   default, RES=[1 1].
%
%   KLIM is a scalar factor for the extension of the interpolation domain.
%   By default, KLIM=1 and the interpolation domain is a rectangle that
%   tightly contains X. Sections of the interpolated surface that protude
%   from the image volume are removed.
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

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010-2011 University of Oxford
% Version: 0.1.2
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
error(nargchk(2, 5, nargin, 'struct'));
error(nargoutchk(0, 4, nargout, 'struct'));

% defaults
if (nargin < 2 || isempty(PARAM))
    PARAM='xy';
end
if (nargin < 3 || isempty(INTERP))
    INTERP='tps';
end
if (nargin < 4 || isempty(res))
    res=[1 1];
end
if (nargin < 5 || isempty(KLIM))
    KLIM=1;
end

%% map the 3D points (x,y,z) to a 2D domain (u,v)

% this is analogous to computing the knot vector for a curve interpolation.
% The idea is that (x,y)->z is not necessarily a function, due to the valve
% folding over. However, we hope that (u,v)->(x,y,z) is a function
switch PARAM
    
    case 'xy'
        
        % (u,v) is simply (x,y)
        em = x(1:2, :);
        
    case 'pca'
        
        % rotate valve points to make the valve surface as horizontal as
        % possible
        m = mean(x, 2);
        em = x - m(:, ones(1, size(x, 2)));
        eigv = pts_pca(em);
        em = eigv' * em;
        em = em(1:2, :);
        
    case 'isomap'
        
        % compute distance matrix
        d = dmatrix(x, x, 'euclidean');
        
        % compute 2-d projection of the 3-d data
        options.dims = 2;
        options.display = 0;
        options.overlay = 0;
        options.verbose = 0;
        em = IsomapII(d, 'k', round(size(x, 2)/3), options);
        em = em.coords{1};
    
    otherwise
        error('Parametrization method not implemented')
end

%% compute interpolation domain

% find box that contains embedded coordinates
emmin = min(em, [], 2);
emmax = max(em, [], 2);

% box size and centroid
delta = emmax - emmin;
boxm = mean([emmax emmin], 2);

% extend the box
emmin = boxm - delta/2*KLIM;
emmax = boxm + delta/2*KLIM;

% generate grid for the embedding box
[gy, gx] = ndgrid(emmin(2):res(1):emmax(2), emmin(1):res(2):emmax(1));


%% compute interpolating surface

% source and target points that will define the warp
%s = em; % don't duplicate data in memory
%t = x; % don't duplicate data in memory

% interpolate
switch INTERP
    case 'tps' % thin-plate spline
        xi = pts_tps_map(em', x', [gx(:) gy(:)]);
    case 'tsi' % Matlab's TriScatteredInterp
        fx = TriScatteredInterp(em(1, :)', em(2, :)', x(1, :)', 'natural');
        fy = TriScatteredInterp(em(1, :)', em(2, :)', x(2, :)', 'natural');
        fz = TriScatteredInterp(em(1, :)', em(2, :)', x(3, :)', 'natural');
        fx = fx(gx, gy);
        fy = fy(gx, gy);
        fz = fz(gx, gy);
        xi = [fx(:) fy(:) fz(:)];
    case 'gridfit'
        fx = gridfit(em(1, :)', em(2, :)', x(1, :)', ...
            emmin(1):res(2):emmax(1), emmin(2):res(1):emmax(2), ...
            'tilesize', 150);
        fy = gridfit(em(1, :)', em(2, :)', x(2, :)', ...
            emmin(1):res(2):emmax(1), emmin(2):res(1):emmax(2), ...
            'tilesize', 150);
        fz = gridfit(em(1, :)', em(2, :)', x(3, :)', ...
            emmin(1):res(2):emmax(1), emmin(2):res(1):emmax(2), ...
            'tilesize', 150);
        xi = [fx(:) fy(:) fz(:)];
    case 'mba' % Multilevel B-Spline Approximation Library
        xi = [...
            mba_surface_interpolation(em(1, :)', em(2, :)', x(1, :)', gx(:), gy(:)) ...
            mba_surface_interpolation(em(1, :)', em(2, :)', x(2, :)', gx(:), gy(:)) ...
            mba_surface_interpolation(em(1, :)', em(2, :)', x(3, :)', gx(:), gy(:))];
    case 'mbae' % MBA extrapolated using a TPS
        % use TPS to interpolate the points, but we are only interested in
        % the edges and corners of the interpolation, where the TPS is
        % extrapolating linearly; also, we decimate the points in the edges
        em2 = [gx(1, 1:end)' gy(1, 1:end)' ; ... % top edge
            gx(2:end, end) gy(2:end, end) ; ... % right edge
            gx(end, 1:end-1)' gy(end, 1:end-1)' ; ... % bottom edge
            gx(2:end-1, 1) gy(2:end-1, 1) ; ... % left edge
            ]';
        xi = pts_tps_map(em', x', em2')';
        x = [x xi(:, 1:10:end)];
        % local support interpolation with boundary conditions provided by
        % the TPS
        em = [em em2(:, 1:10:end)];
        xi = [...
            mba_surface_interpolation(em(1, :)', em(2, :)', x(1, :)', gx(:), gy(:)) ...
            mba_surface_interpolation(em(1, :)', em(2, :)', x(2, :)', gx(:), gy(:)) ...
            mba_surface_interpolation(em(1, :)', em(2, :)', x(3, :)', gx(:), gy(:))];
    otherwise
        error('Interpolation method not implemented')
end
