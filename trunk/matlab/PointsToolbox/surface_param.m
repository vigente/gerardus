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
%     XY plane:
%                 'xy' [default]
%                 'pca'
%                 'isomap'
%
%     Unit sphere:
%                 'sphproj
%                 'cald'
%                 'sphisomap'
%
%     Parametrizations details:
%
%   * 'xy' (default): Projection on the xy-plane, i.e. xy-coordinates of
%     the points in X.
%
%   * 'pca': Projection on the plane defined by the two eigenvectors with
%     the largest eigenvalues.
%
%   * 'isomap': Isometric mapping (Tenenbaum et al. [1][2]). This is a
%     non-linear method that can be applied to complex open surfaces that
%     bend on themselves or are rolled-up. This option uses the
%     third-party function IsomapII. Because isomap uses Dijkstra's
%     shortest-path algorithm internally, it suffers strongly from
%     metrication errors in regular meshes. On the other hand, it can be
%     applied to any adjacency graph, not only meshes.
%
%     Small neighbourhoods capture higher curvature in the surface, but
%     too small neighbourhoods create several connected components and
%     will produce an error.
%
%     PARAM.d: [Opt] Distance matrix between points. The parameter can have
%       the following formats:
%       (1) a full N x N matrix with distances between all pairs of points.
%       In this case, you can define the local neighbourhoods using
%       PARAM.neigh and PARAM.size, so that some of the connections in the
%       full matrix will be deleted and replaced by Disjkstra's
%       shortest-path lengths.
%       (2) a sparse N x N matrix (missing entries are treated as Inf). In
%       this case, you usually want to make PARAM.size=Inf, because the
%       local neighbourhood is already given by the sparse matrix. The
%       missing connections will be replaced internally by Disjkstra's
%       shortest-path lengths.
%       (3) the name of a function (e.g. 'd_fun') that takes one argument,
%       i, and returns a row vector containing the distances from all N
%       points to point i.
%
%       If not provided, PARAM.d is the full Euclidean matrix of distances
%       between every pair of points.
%
%       As a particular case, note that in case of a triangular mesh, the
%       sparse local distance that corresponds to the mesh can be computed
%       using dmatrix_mesh().
%
%     PARAM.neigh: [Opt] Type of neighbourhood ('epsilon' or 'k'). 
%     PARAM.size:  [Opt] Size of neighbourhood.
%
%       Each point in X is directly connected only to neighbours within a
%       distance PARAM.size (for 'epsilon'), or the PARAM.size closest
%       neighbours (for 'k'), according to the distances in PARAM.d. This
%       allows to further restrict the neighbourhoods in PARAM.d.
%
%     PARAM.options: [Opt] Extra arguments passed to IsomapII. Most users
%       won't need to change this. For details, see help to IsomapII().
%
%   * 'cmdsmap': Open surface classical MDS mapping improved with a Fast
%     Marching method (Zigelman et al. [3]). The fast marching algorithm
%     doesn't suffer from metrication errors on regular meshes as badly as
%     Dijkstra's. On the other hand, the current implementation is limited
%     to triangular meshes, while Dijkstra's will work with any adjacency
%     graph. This option uses classical MDS on a full distance matrix, so
%     it's not valid for local neighbourhoods or closed surfaces. For local
%     neighbourhoods, see 'lmdscale' below. For closed surfaces, see 
%
%     PARAM.d or PARAM.tri must be provided. If PARAM.d is provided,
%     PARAM.tri and PARAM.options are ignored. Otherwise, PARAM.d is
%     computed from PARAM.tri and PARAM.options.
%
%     PARAM.d: [Req/Opt] This is expected to be the full distance matrix
%       between mesh vertices that can be computed with
%         [~, param.d] = dmatrix_mesh(param.tri, x, 'fastmarching');
%       Note that for comparison purposes, Dijkstra's algorithm can be used
%       instead of Fast Marching
%         [~, param.d] = dmatrix_mesh(param.tri, x, 'dijkstra');
%       The latter is equivalent to the 'isomap' method above.
%
%     PARAM.tri: [Req/Opt] 3-column matrix. Each row contains the 3 nodes
%       that form one triangular facet in the mesh.
%
%     PARAM.options: [Opt] Extra arguments passed to the fast marching
%       method:
%       - options.constraint_map: Exploration radius to reduce the set of
%         explored points. Only points with current distance smaller than L
%         will be expanded. Set some entries of L to -Inf to avoid any
%         exploration of these points.
%       - options.W: non-uniform speed.
%       - options.end_points: stop when these points are reached.
%       - options.nb_iter_max: stop when a given number of iterations is
%         reached.
%       - options.heuristic: heuristic that tries to guess the distance
%         that remains from a given node to a given target. This is an
%         array of same size as W.
%
%
%   * 'sphproj':   Direct projection on unit sphere centered on point set
%     centroid.
%
%     OUT.rmed:    Median value of the radii of all points in the point
%                    set.
%
%   * 'cald':      Li Shen's Control Area and Length Distortions (CALD)
%     spherical parametrization [4][5].
%
%     PARAM.MeshGridSize: The interpolation mesh used to smooth the
%                  spherical parametrization has length
%                  2*PARAM.MeshGridSize+1. By default,
%                  PARAM.MeshGridSize = 50.
%
%     PARAM.MaxSPHARMDegree: Degree of the spherical harmonics used for
%                 interpolation to smooth the spherical parametrization.
%                 By default, PARAM.MaxSPHARMDegree = 6.
%
%     PARAM.Tolerance: Scalar. In the interpolation smoothing, grid
%                 points with a height larger than PARAM.Tolerance*gmin,
%                 where gmin is the lowest point in the grid, will be
%                 truncated. By default, PARAM.Tolerance = 2.
%
%     PARAM.Smoothing: Scalar. Before parametrization smoothing,
%                parametrization triangle areas relative to triangular
%                mesh areas are smoothed by elevating to
%                (1/PARAM.Smoothing). By default, PARAM.Smoothing = 2.
%
%     PARAM.Iteration: Maximum number of smoothing iterations. Smoothing
%                will stop before reaching PARAM.Iteration if the
%                algorithm converges to a solution. By default,
%                PARAM.Iteration = 100.
%
%     PARAM.LocalIteration: Number of smoothing iteration in the local
%               smoothing algorithm. By default, PARAM.LocalIteration = 10.
%
%   * 'sphisomap': Our extension to the Isomap method for closed
%     surfaces. The result is a parametrization on the unit sphere.
%
%     PARAM.d:     See 'isomap' above.
%     PARAM.neigh: See 'isomap' above.
%     PARAM.size:  See 'isomap' above.
%
%     PARAM.init:   Initial parametrization. It can be a 2-column matrix,
%                   in which case those are the coordinates directly
%                   used, or a string:
%
%                   'sphproj': (default) Surface points are projected
%                   onto the sphere along the radii to their centroid.
%
%                   'random': Uniformly random distribution of points on
%                   the sphere. Counter-intuitively, this can work better
%                   than sphproj, because it may prevent MDS from getting
%                   trapped in local minima.
%
%      PARAM.*:     Any other fields are passed to function
%                   smdscale(...,OPT=PARAM). This allows to pass e.g.
%                   stopping conditions to the spherical MDS algorithm.
%                   See "help smdscale" for all options.
%
%      OUT.err:     Frobenius norm of the distance matrix error
%                   D_xy - D_uv at each step of the optimisation
%                   algorithm (moving each point counts as a step).
%
%
%
% [1] J.B. Tenenbaum, V. de Silva and J.C. Langford, "A Global Geometric
% Framework for Nonlinear Dimensionality Reduction", Science 290(5500):
% 2319-2323, 2000.
%
% [2] Isomap Homepage, http://isomap.stanford.edu/
%
% [3] G. Zigelman, R. Kimmel, and N. Kiryati, "Texture mapping using
% surface flattening via multidimensional scaling,” Visualization and
% Computer Graphics, IEEE Transactions on, vol. 8, no. 2, pp. 198–207,
% 2002.
%
% [4] Shen and Makedon, "Spherical mapping for processing of 3D closed
% surfaces", Image and Vision Computing, 24(7):743–761, 2006.
% http://www.sciencedirect.com/science/article/pii/S0262885606000485
%
% [5] CALD and SPHARM-MAT homepage,
% http://www.iupui.edu/~shenlab/software.html

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
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
        
        % IsomapII takes the following distance matrices:
        %  (1) a full N x N matrix (as in isomap.m)  
        %  (2) a sparse N x N matrix (missing entries are treated as INF)
        %  (3) the name of a function (e.g. 'd_fun') that takes
        %        one argument, i, and returns a row vector containing the 
        %        distances from all N points to point i.        
        
        % default parameters
        if (~isfield(param, 'neigh') || isempty(param.neigh))
            param.neigh = 'epsilon';
        end
        if (~isfield(param, 'size') || isempty(param.size))
            param.size = Inf;
        end

        % if distance matrix is not provided by the user, it is computed
        % internally as the full Euclidean distance matrix
        if (~isfield(param, 'd') || isempty(param.d))
            param.d = dmatrix(x', x', 'euclidean');
        end
        
        % compute 2-d projection of the 3-d data
        if (~isfield(param, 'options') || isempty(param.options))
            param.options.dims = 2;
            param.options.display = 0;
            param.options.overlay = 0;
            param.options.verbose = 0;
            param.options.dijkstra = 1;
        end
        uv = IsomapII(param.d, param.neigh, param.size, param.options);
        uv = uv.coords{1}';
        if (size(uv, 1) ~= size(x, 1))
            error('The neighbourhood size (param.size) is too small to connect all points')
        end
        
    case 'cmdsmap'
        
        % if distance matrix is not provided
        if (~isfield(param, 'd') || isempty(param.d))
            
            % the user hasn't provided the mesh triangulation either
            if (~isfield(param, 'tri') || isempty(param.tri))
                error('Either PARAM.d or PARAM.tri must be provided')
            end
            
            % default for fast marching method options
            if (~isfield(param, 'options') || isempty(param.options))
                param.options = [];
            end

            % remove fast marching options that don't apply to this method
            if (isfield(param.options, 'constraint_map'))
                warning('Gerardus:InvalidInputArg', 'PARAM.options.constraint_map ignored, this method only works with a full distance matrix')
                param.options = rmfield(param.options, 'constraint_map');
            end
            if (isfield(param.options, 'end_points'))
                warning('Gerardus:InvalidInputArg', 'PARAM.options.end_points ignored, this method only works with a full distance matrix')
                param.options = rmfield(param.options, 'end_points');
            end
            if (isfield(param.options, 'nb_iter_max'))
                warning('Gerardus:InvalidInputArg', 'PARAM.options.nb_iter_max ignored, this method only works with a full distance matrix')
                param.options = rmfield(param.options, 'nb_iter_max');
            end
            
            % compute full distance matrix from the mesh
            [~, param.d] = dmatrix_mesh(param.tri, x, 'fastmarching', param.options);
            
        end
        
        % if the user provided a sparse distance matrix, we assume that
        % it's a mistake
        if (issparse(param.d))
            error('PARAM.d must be a full distance matrix, not a sparse mesh adjacency-distance matrix');
        end
        
        % use classic MDS for the optimization
        uv = cmdscale(param.d);
        
        % the output parametrization is 2D
        uv = uv(:, 1:2);
            
    case 'lmdscale'
        
        % extract parameters for the fast marching algorithm
        options = [];
        if (~isfield(param, 'W') || isempty(param.W))
            % do nothing
        else
            options.W = param.W;
        end
        if (~isfield(param, 'end_points') || isempty(param.end_points))
            % do nothing
        else
            options.end_points = param.end_points;
        end
        if (~isfield(param, 'nb_iter_max') || isempty(param.nb_iter_max))
            % do nothing
        else
            options.nb_iter_max = param.nb_iter_max;
        end
        if (~isfield(param, 'heuristic') || isempty(param.heuristic))
            % do nothing
        else
            options.heuristic = param.heuristic;
        end
        if (~isfield(param, 'constraint_map') || isempty(param.constraint_map))
            % do nothing
        else
            options.constraint_map = param.constraint_map;
        end
        
        % local neighbourhoods can be used in two ways:
        %
        % * the user provides a distance matrix (typically, sparse) with 0
        % values for vertices outside the local neighbourhood.
        %
        % * the user does not provide a distance matrix, but specifies
        %   the neighbourhood size with param.size
        %
        % if distance matrix is not provided by the user, it is computed
        % extract parameters
        if (~isfield(param, 'd') || isempty(param.d))
            
            % number of vertices in the mesh
            N = size(x, 1);
            
            % compute full distance matrix from the mesh
            d = zeros(N);
            for I = 1:N
                d(:, I) = perform_fast_marching_mesh(x, param.tri, I, options);
            end
            
        end
        

    case 'sphproj'
        
        % init output
        uv = zeros(size(x, 1), 1);
        
        % project 3D Cartesian coordinates onto the surface of a sphere
        % lon, lat given in radians, centered around 0
        [uv(:, 2), uv(:, 1), r] = cart2sph(...
            x(:, 1) - mean(x(:, 1)), ...
            x(:, 2) - mean(x(:, 2)), ...
            x(:, 3) - mean(x(:, 3)));
        
        % we use the median radius of the points in the configuration as the radius
        % of the sphere that best can contains their projections
        out.rmed = median(r);
        
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
                    [lat, lon] = proj_on_sphere(x);
                    
                case 'random'
                    
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
        
        % estimate the size of the sphere so that it can accommodate the
        % distance matrix
        %
        % first, we find the furthest point from each point. The median of
        % the corresponding distances give an estimate of the half
        % circumference of the sphere. The radius = half_circumference/pi
        % because circumference = 2*pi*radius
        sphrad = median(max(param.d)) / pi;
        
        % compute the parametrization that optimises isometry
        [lat, lon, out.err, out.stopCondition, out.dsph, out.sphrad] = ...
            smdscale(param.d, sphrad, lat, lon, rmfield(param, 'd'));
        
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
        
    case 'cald'
        
        confs.vars = {'MeshGridSize', 'MaxSPHARMDegree', 'Tolerance', ...
            'Smoothing', 'Iteration', 'LocalIteration', 't_major', ...
            'SelectDiagonal', 'OutDirectory'};
        confs.args = [1 1 1 1 1 1 10 10 200];
        confs.inFilter = {
            '*_bim.mat;*_fix.mat;*_obj.mat' 'Binary and Mesh Objects (*_bim.mat;*_fix.mat;*_obj.mat)'
            '*_obj.mat'                     'Mesh Objects (*_obj.mat)'
            '*_bim.mat;*_fix.mat'           'Binary Objects (*_bim.mat, *_fix.mat)'
            '*.mat'                         'All Objects (*.mat)'
            };
        confs.default = {'50'  '6'  '2'  '2'  '100'  '10'  ''  ''  './'};
        if (~isfield(param, 'MeshGridSize') || isempty(param.MeshGridSize))
            confs.MeshGridSize = 50;
        else
            confs.MeshGridSize = param.MeshGridSize;
        end
        if (~isfield(param, 'MaxSPHARMDegree') || isempty(param.MaxSPHARMDegree))
            confs.MaxSPHARMDegree = 6;
        else
            confs.MaxSPHARMDegree = param.MaxSPHARMDegree;
        end
        if (~isfield(param, 'Tolerance') || isempty(param.Tolerance))
            confs.Tolerance = 2;
        else
            confs.Tolerance = param.Tolerance;
        end
        if (~isfield(param, 'Smoothing') || isempty(param.Smoothing))
            confs.Smoothing = 2;
        else
            confs.Smoothing = param.Smoothing;
        end
        if (~isfield(param, 'Iteration') || isempty(param.Iteration))
            confs.Iteration = 100;
        else
            confs.Iteration = param.Iteration;
        end
        if (~isfield(param, 'LocalIteration') || isempty(param.LocalIteration))
            confs.LocalIteration = 10;
        else
            confs.LocalIteration = param.LocalIteration;
        end
        confs.t_major = 'x';
        confs.SelectDiagonal = 'ShortDiag';
        confs.OutDirectory = '.';
        
        % initial parametrization
        [sph_verts, name3] = initParamCALD(x, param.tri, '', confs);
        
        % smooth the parameterization
        [~, ~, sph_verts, new_name] = smootheCALD(x, param.tri, sph_verts, name3, confs);
        
        % convert xyz coordinates to spherical coordinates
        [lon, lat] = cart2sph(sph_verts(:, 1), sph_verts(:, 2), sph_verts(:, 3));
        
        % delete files created by the SPHARM functions
        if exist(new_name, 'file')
            delete(new_name) % CALD_smo.mat
        end
        if exist('initParamCALD', 'file')
            rmdir('initParamCALD', 's')
        end
        
        % create output
        uv = [lat lon];
        
    otherwise
        
        error('Parametrization method not implemented')
end
