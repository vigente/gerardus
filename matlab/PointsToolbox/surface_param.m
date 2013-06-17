function [uv, out] = surface_param(x, param)
% SURFACE_PARAM  Parametrization of a scattered set of points that belong
% to a surface
%
%   This function takes a scattered data set X in 3D that is supposed to
%   represent a sampling of a surface, and finds a 2D parameterisation U, V
%   for them. The surface can be open or closed, in which case the
%   parametrization domain will be the R^2 plane or the unit sphere,
%   respectively.
%
%   This function provides several parametrization methods, some
%   implemented by us, some by other authors.
%
% [UV, OUT] = surface_param(X, PARAM)
%
%   X is a 3-column matrix. Each row has the coordinates of a point that
%   belongs to the surface we want to interpolate.
%
%   UV is a 2-column matrix with the parameterisation of X, (U, V)->X. In
%   planar parameterisations, UV has units of meters. In spherical
%   parameterisations, UV=[LAT, LON], in units of radians.
%
%   OUT is a struct with the extra output from some of the parametrization
%   methods.
%
%   PARAM is a string or struct that selects the parametrization method
%   and possibly arguments for the method.
%
%     Parametrizations summary:
%
%     PARAM =
%                 'xy' [default]
%                 'pca'
%
%     PARAM.type =
%                 'isomap'
%                 'sphisomap'
%
%     Parametrizations details:
%
%     * 'xy' (default): Projection on the xy-plane, i.e. xy-coordinates of
%       the points in X.
%
%     * 'pca': Projection on the plane defined by the two eigenvectors with
%       the largest eigenvalues.
%
%     * 'isomap': Isometric mapping (Tenenbaum et al. [1][2]). This is a
%       non-linear method that can be applied to complex open surfaces that
%       bend on themselves or are rolled-up. This option uses the
%       third-party function IsomapII.
%
%       PARAM.d:     [Opt] Distance matrix between points. By default,
%                    PARAM.d is the Euclidean matrix of distances between
%                    every pair of points.
%
%                    Providing the matrix allows the user to define
%                    specific neighbourhoods and topologies, or even use
%                    non-Euclidean distances.
%
%       PARAM.neigh: [Opt] Type of neighbourhood ('epsilon' or 'k'). 
%       PARAM.size:  [Opt] Size of neighbourhood.
%
%                    Each point in X is directly connected only to
%                    neighbours within a distance PARAM.size (for
%                    'epsilon'), or the PARAM.size closest neighbours (for
%                    'k'), according to the distances in PARAM.d. This
%                    allows to further restrict the neighbourhoods in
%                    PARAM.d. Note that in the 'k' case, in order to obtain
%                    a symmetric local neighbourhood distance matrix, the
%                    number of nearest neighbours may be increased.
%
%                    Small neighbourhoods capture higher curvature in the
%                    surface, but too small neighbourhoods create several
%                    connected components and will produce an error.
%
%       PARAM.options: [Opt] Extra arguments passed to IsomapII. Most users
%                    won't need to change this. For details, see
%                    help IsomapII.
%
%     * 'sphisomap': Our extension to the Isomap method for closed
%       surfaces. The result is a parametrization on the unit sphere.
%
%       PARAM.d:     See 'isomap' above.
%       PARAM.neigh: See 'isomap' above.
%       PARAM.size:  See 'isomap' above.
%
%       PARAM.init:   Initial parametrization. It can be a 2-column matrix,
%                     in which case those are the coordinates directly
%                     used, or a string:
%
%                     'sphproj': (default) Surface points are projected
%                     onto the sphere along the radii to their centroid.
%
%                     'random': Uniformly random distribution of points on
%                     the sphere. Counter-intuitively, this can work better
%                     than sphproj, because it may prevent MDS from getting
%                     trapped in local minima.
%
%        PARAM.*:     Any other fields are passed to function
%                     smdscale(...,OPT=PARAM). This allows to pass e.g.
%                     stopping conditions to the spherical MDS algorithm.
%                     See "help smdscale" for all options.
%
%        OUT.err:     Frobenius norm of the distance matrix error
%                     D_xy - D_uv at each step of the optimisation
%                     algorithm (moving each point counts as a step).
%
%
%
% [1] J.B. Tenenbaum, V. de Silva and J.C. Langford, "A Global Geometric
% Framework for Nonlinear Dimensionality Reduction", Science 290(5500):
% 2319-2323, 2000.
%
% [2] Isomap Homepage, http://isomap.stanford.edu/

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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
narginchk(2, 2);
nargoutchk(0, 2);

if (size(x, 2) ~= 3)
    error('X must be a 3-column matrix')
end

if ischar(param)
    aux = param; param = [];
    param.type = aux;
end

out = [];

% defaults
if (isempty(param) || ~isfield(param, 'type') ...
        || isempty(param.type))
    param.type = 'xy';
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

        % if distance matrix is not provided by the user, it is computed
        % internally as the full Euclidean distance matrix. Note that we
        % can still define local neighbourdhoods with the size parameter
        if (~isfield(param, 'd') || isempty(param.d))
            param.d = dmatrix(x', x', 'euclidean');
        end
        
        % compute 2-d projection of the 3-d data
        options.dims = 2;
        options.display = 0;
        options.overlay = 0;
        options.verbose = 0;
        options.dijkstra = 1;
        uv = IsomapII(param.d, param.neigh, param.size, options);
        uv = uv.coords{1}';
        if (size(uv, 1) ~= size(x, 1))
            error('The neighbourhood size (param.size) is too small to connect all points')
        end
        
    case 'sphisomap'

        % if distance matrix is not provided by the user, it is computed
        % internally as the full Euclidean distance matrix. Note that we
        % can still define local neighbourdhoods with the size parameter
        if (~isfield(param, 'd') || isempty(param.d))
            param.d = dmatrix(x', x', 'euclidean');
        end
        
        % number of points to embed
        N = size(x, 1);
        
        % reduce the fully connected matrix of Euclidean distances to a
        % sparse local neighbourhood matrix
        if (param.size < Inf)
            switch param.neigh
                
                case 'epsilon'
                    
                    param.d(param.d > param.size) = 0;
                    
                case 'k'
                    
                    % the code here reproduces the behaviour of the
                    % equivalent in IsomapII, only here by columns instead
                    % of by rows. The symmetrization step makes the result
                    % the same
                    
                    % column by column, sort neighbours in oder of
                    % increasing distance
                    [~, idx] = sort(param.d, 1, 'ascend');
                    
                    % we want to keep only the k+1 nearest neighbours (k
                    % neighbours + the point itself)
                    idx = idx(1:param.size+1, :);
                    
                    % all the indices just computed refer to their
                    % corresponding column, i.e., they are all from 1 to N.
                    % Make the indines refer to the whole matrix
                    idx = idx + repmat((0:N-1)*N, param.size+1, 1);
                    
                    % keep only the k-nearest neighbour distances
                    aux = zeros(size(param.d));
                    aux(idx) = param.d(idx);
                    param.d = aux;
                    
                    % make sure the distance matrix is symmetric
                    param.d = min(param.d, param.d');
                    
                otherwise
                    error('Invalid neighbourhood type')
            end
        end
        param.d = sparse(param.d);
        
        % defaults
        if (~isfield(param, 'maxiter') || isempty(param.maxiter))
            param.maxiter = 50;
        end
        if (~isfield(param, 'init') || isempty(param.init))
            param.init = 'sphproj';
        end

        % sMDS initialisation
        if ischar(param.init) % initial parametrization as a string
            switch param.init
                case 'sphproj'
                    
                    % initial guess for the sphere embedding by simply
                    % projecting each point along the radius to the
                    % centroid
                    [lat, lon, sphrad] = proj_on_sphere(x);
                    
                case 'random'
                    
                    % kludge: we use this to estimate the sphere's radius
                    [~, ~, sphrad] = proj_on_sphere(x);
                    
                    % counter-intuitively, starting from a completely
                    % random distribution usually gives good results
                    lat = (rand(N, 1) - .5) * pi;
                    lon = 2 * (rand(N, 1) - .5) * pi;
                    
                otherwise
                    error('Not valid initialisation for spherical Isomap parameterisation')
            end
        else % initial parametrization as a point set
            if (size(param.init, 1) ~= N || size(param.init, 1) ~= 2)
                error('Initial parametrization provided, but it is not a 2-column matrix with the same number of points as X')
            end
            
            lat = param.init(:, 1);
            lon = param.init(:, 2);
        end
        
        % compute the full matrix of distances from the local neighbourhood
        % one
        param.d = dijkstra(param.d, 1:N);
        
        % compute the parametrization that optimises isometry
        [lat, lon, out.err] = smdscale(param.d, sphrad, lat, lon, param);
            
        
%         % DEBUG: plot parametrization
%         [xsph, ysph, zsph] = sph2cart(lon, lat, sphrad);
%         [~, xyzsph] = procrustes(x, [xsph ysph zsph], 'Scaling', false);
%         tri = DelaunayTri(xyzsph);
%         tri = freeBoundary(tri);
%         trisurf(tri, x(:, 1), x(:, 2), x(:, 3));
%         axis equal
%         
%         % DEBUG: plot cloud of points for parameterisation
%         plot3(xsph, ysph, zsph, '.')
%         axis equal
        
        % create output with the parameterisation values for each point
        % lat = em(1, :)
        % lon = em(2, :)
        uv = [lat lon];
        
    otherwise
        
        error('Parametrization method not implemented')
end
