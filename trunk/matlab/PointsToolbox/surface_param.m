function [uv, out] = surface_param(x, param)
% SURFACE_PARAM  2D parametrization of a scattered set of 3D points or of
% the vertices of a triangular mesh, for open and closed surfaces.
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
%                 'xy' [default]: Simple projection on XY.
%                 'pca': Projection on main PCA plane.
%                 'isomap': Isomap (with Dijkstra's shortest path).
%                 'cmdsmap': MDSmap, global neighbourhood, classical MDS.
%                 'lmdscale': MDSmap, local neighbourhood, local MDS.
%
%     Unit sphere:
%                 'sphproj': Simple projection on sphere around centroid.
%                 'cald'
%                 'sphisomap'
%
%=========================== OPEN SURFACES ================================
%
%   * 'xy' (default): Projection on the xy-plane, i.e. xy-coordinates of
%     the points in X.
%
%--------------------------------------------------------------------------
%
%   * 'pca': Projection on the plane defined by the two eigenvectors with
%     the largest eigenvalues.
%
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
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
%     If PARAM.d is not provided:
%
%     PARAM.tri: [Req/Opt] 3-column matrix. Each row contains the 3 nodes
%       that form one triangular facet in the mesh.
%
%     PARAM.dmethod: [Opt] String with the method to compute distances
%       between non-neighbours: 'fastmarching' (default) or 'dijkstra'. The
%       latter is faster, but suffers much more from metrication errors.
%
%     PARAM.options: [Opt] Extra arguments passed to the fast marching
%       method. From the help to perform_fast_marching_mesh():
%       - options.W: non-uniform speed.
%       - options.heuristic: heuristic that tries to guess the distance
%         that remains from a given node to a given target. This is an
%         array of same size as W.
%       - parameters that create a local neighbourhood are ignored with a
%         warning.
%
%--------------------------------------------------------------------------
%
%   * 'lmdscale': Open surface local neighbourhood MDS mapping (Schwartz,
%     Wolfson and Shaw [4, 5]), improved with a Fast Marching method
%     (Zigelman et al. [3]). This is like 'cmdsmap' above, but with local
%     neighbourhoods. Local neighbourhoods were already proposed in the
%     first publications [4, 5] on what we are calling MDSmap. We use our
%     own implementation of this idea, lmdscale().
%
%     PARAM.d or PARAM.tri must be provided. If PARAM.d is provided,
%     PARAM.tri and PARAM.options are ignored. Otherwise, PARAM.d is
%     computed from PARAM.tri and PARAM.options.
%
%     PARAM.d: [Req/Opt] Can be a sparse or full distance matrix. Distances
%       greater than 0 reflect that two vertices are within a common local
%       neighbourhood.
%         [~, param.d] = dmatrix_mesh(param.tri, x, 'fastmarching');
%
%     If PARAM.d is not provided:
%
%     PARAM.tri: [Req/Opt] 3-column matrix. Each row contains the 3 nodes
%       that form one triangular facet in the mesh.
%
%     PARAM.dmethod: [Opt] String with the method to compute distances
%       between non-neighbours: 'fastmarching' (default) or 'dijkstra'. The
%       latter is faster, but suffers much more from metrication errors,
%       and produces only full matrices in the current implementation, so
%       it's not valid for local neighbourhoods.
%
%     PARAM.options: [Opt] Extra arguments passed to the fast marching
%       method. From the help to perform_fast_marching_mesh():
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
%     PARAM.options2: [Opt] Extra arguments passed to the local
%       neighbourhood MDS algorithm. See help to lmdscale().
%       - opt.MaxIter: (default = 20) maximum number of iterations we allow
%         the optimisation algorithm. Each iteration consists of an entire
%         sweep of all points on the sphere.
%       - opt.MaxInc:  inc is the Euclidean distance that each output point
%         is moved in an iteration. The algorithm will stop if no point has
%         been moved more than MaxInc.
%
%     OUT.stopCondition: cell-array with the condition/s that made the
%       algorithm stop, in string form.
%
%     OUT.err:  struct with several error measures at each algorithm
%     iteration. See help to lmdscale() for details.
%
%
%========================== CLOSED SURFACES ===============================
%
%   * 'sphproj':   Direct projection on unit sphere centered on point set
%     centroid.
%
%     OUT.medrad:  Median value of the radii of all points in the point
%       set.
%
%--------------------------------------------------------------------------
%
%   * 'cald': Li Shen's Control Area and Length Distortions (CALD)
%     spherical parametrization [6, 7].
%
%     PARAM.tri: [Req] 3-column matrix. Each row contains the 3 nodes that
%       form one triangular facet in the mesh.
%
%     PARAM.options: [Opt] Extra arguments passed to CALD:
%       - options.MeshGridSize: The interpolation mesh used to smooth the
%         spherical parametrization has length 2*PARAM.MeshGridSize+1. By
%         default, MeshGridSize = 50.
%       - options.MaxSPHARMDegree: Degree of the spherical harmonics used
%         for interpolation to smooth the spherical parametrization. By
%         default, MaxSPHARMDegree = 6.
%       - options.Tolerance: Scalar. In the interpolation smoothing, grid
%         points with a height larger than PARAM.Tolerance*gmin, where gmin
%         is the lowest point in the grid, will be truncated. By default,
%         Tolerance = 2.
%       - options.Smoothing: Scalar. Before parametrization smoothing,
%         parametrization triangle areas relative to triangular mesh areas
%         are smoothed by elevating to (1/PARAM.Smoothing). By default,
%         Smoothing = 2.
%       - options.Iteration: Maximum number of smoothing iterations.
%         Smoothing will stop before reaching Iteration if the algorithm
%         converges to a solution. By default, Iteration = 100.
%       - options.LocalIteration: Number of smoothing iteration in the
%         local smoothing algorithm. By default, LocalIteration = 10.
%
%--------------------------------------------------------------------------
%
%   * 'smdscale': Local neighbourhood Multidimensional scaling on a sphere.
%     Our extension to lmdscale() 
%
%     PARAM.d or PARAM.tri must be provided. If PARAM.d is provided,
%     PARAM.tri and PARAM.options are ignored. Otherwise, PARAM.d is
%     computed from PARAM.tri and PARAM.options.
%
%     PARAM.d: [Req/Opt] Can be a sparse or full distance matrix. Distances
%       greater than 0 reflect that two vertices are within a common local
%       neighbourhood.
%         [~, param.d] = dmatrix_mesh(param.tri, x, 'fastmarching');
%
%     If PARAM.d is not provided:
%
%     PARAM.tri: [Req/Opt] 3-column matrix. Each row contains the 3 nodes
%       that form one triangular facet in the mesh.
%
%     PARAM.dmethod: [Opt] String with the method to compute distances
%       between non-neighbours: 'fastmarching' (default) or 'dijkstra'. The
%       latter is faster, but suffers much more from metrication errors,
%       and produces only full matrices in the current implementation, so
%       it's not valid for local neighbourhoods.
%
%     PARAM.options: [Opt] Extra arguments passed to the fast marching
%       method. From the help to perform_fast_marching_mesh():
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
%     PARAM.init: [Opt] Parameterization initialisation:
%       - 'sphproj' (default): Direct projection on unit sphere centered on
%         point set centroid, as the option above.
%       - 'random': Random distribution of points on the sphere. This is
%         useful if 'sphproj' produces a fold-over of points that
%         corresponds to a local minimum.
%       - [lat lon]: A 2-column matrix with the initial latitude and
%         longitude coordinates of the parameterization.
%
%     PARAM.sphrad: [Opt] Parameterization sphere's radius. If not
%       provided, it is computed as the median of the distance of the
%       points in X to their centroid.
%
%     PARAM.options2: [Opt] Extra arguments passed to the spherical local
%       neighbourhood MDS algorithm. See help to smdscale().
%       - opt.MaxIter: (default = 20) maximum number of iterations we allow
%         the optimisation algorithm. Each iteration consists of an entire
%         sweep of all points on the sphere.
%       - opt.MaxAlpha:  inc is the angular distance in radians that each
%         output point is moved in an iteration. The algorithm will stop if
%         no point has been moved more than MaxAlpha.
%
%     OUT.stopCondition: cell-array with the condition/s that made the
%       algorithm stop, in string form.
%
%     OUT.err: Error measures vs. iterations:
%       - err.rawstress: Raw stress (sigma_r), as defined above.
%       - err.stress1:   Stress-1 (sigma_1).
%         sigma_1 = sqrt(sigma_r / sum_{i,j} (d_ij)^2)
%       - err.maxalpha:  Maximum angular displacement of any point in each
%         iteration.
%       - err.medalpha:  Median angular displacement of all points in each
%         iteration.
%
%--------------------------------------------------------------------------
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
% [4] E. Wolfson and E.L. Schwartz, "Computing minimal distances on
% polyhedral surfaces", IEEE Pattern Analysis and Machine Intelligence,
% 11(9):1001-1005, 1989.
%
% [5] E.L. Schwartz, A. Shaw, E. Wolfson, "A numerical solution to the
% generalized mapmaker’s problem: flattening nonconvex polyhedral", IEEE
% Pattern Analysis and Machine Intelligence, 11(9):1005-1008, 1989.
%
% [6] Shen and Makedon, "Spherical mapping for processing of 3D closed
% surfaces", Image and Vision Computing, 24(7):743–761, 2006.
% http://www.sciencedirect.com/science/article/pii/S0262885606000485
%
% [7] CALD and SPHARM-MAT homepage,
% http://www.iupui.edu/~shenlab/software.html

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
% Version: 0.6.0
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
    
%--------------------------------------------------------------------------

    case 'xy'
        
        % check arguments
        nargoutchk(0, 1);
        
        % (u,v) is simply (x,y)
        uv = x(:, 1:2);
        
%--------------------------------------------------------------------------

    case 'pca'
        
        % this function doesn't provide the "out" output
        nargoutchk(0, 1);
        
        % rotate valve points to make the valve surface as horizontal as
        % possible
        m = mean(x, 1);
        uv = x - m(ones(1, size(x, 1)), :);
        eigv = pts_pca(uv');
        uv = uv * eigv;
        uv = uv(:, 1:2);
        
%--------------------------------------------------------------------------

    case 'isomap'
        
        % check arguments
        nargoutchk(0, 1);

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
        
%--------------------------------------------------------------------------

    case 'cmdsmap'
        
        % check arguments
        nargoutchk(0, 1);

        % if distance matrix is not provided we need to compute it from the
        % mesh description
        if (~isfield(param, 'd') || isempty(param.d))
            
            % but if the user hasn't provided the mesh triangulation
            % either, we cannot, so we give an error
            if (~isfield(param, 'tri') || isempty(param.tri))
                error('Either PARAM.d or PARAM.tri must be provided')
            end
            
            % default for fast marching method options
            if (~isfield(param, 'dmethod') || isempty(param.dmethod))
                param.dmethod = 'fastmarching';
            end
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
            [~, param.d] ...
                = dmatrix_mesh(param.tri, x, param.dmethod, param.options);
            
            if (issparse(param.d))
                error('Assertion fail: param.d should be a ');
            end
            
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
            
%--------------------------------------------------------------------------

    case 'lmdscale'
        
        % check arguments
        nargoutchk(0, 2);

        % if distance matrix is not provided we need to compute it from the
        % mesh description
        if (~isfield(param, 'd') || isempty(param.d))
            
            % but if the user hasn't provided the mesh triangulation
            % either, we cannot, so we give an error
            if (~isfield(param, 'tri') || isempty(param.tri))
                error('Either PARAM.d or PARAM.tri must be provided')
            end
            
            % default for distance method
            if (~isfield(param, 'options') || isempty(param.options))
                param.options = [];
            end
            
            % default for fast marching method options
            if (~isfield(param, 'dmethod') || isempty(param.dmethod))
                param.dmethod = 'fastmarching';
            end
            if (~isfield(param, 'options') || isempty(param.options))
                param.options = [];
            end
            
            % compute neighbourhood distance matrix from the mesh
            [~, param.d] ...
                = dmatrix_mesh(param.tri, x, param.dmethod, param.options);
            
        end
        
        % at this point, the distance matrix can be full or sparse, both
        % are valid, unlike in the cases above, where it had to be full
        
        % compute parameterization using local neighbourhood method
        [u, v, out.stopCondition, out.err] ...
            = lmdscale(param.d, [], [], param.options2);
        uv = [u, v];

%--------------------------------------------------------------------------

    case 'sphproj'
        
        % check arguments
        nargoutchk(0, 2);

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
        out.medrad = median(r);
        
%--------------------------------------------------------------------------

    case 'cald'
        
        % check arguments
        nargoutchk(0, 1);

        % required parameters
        if (~isfield(param, 'tri') || isempty(param.tri))
            error('PARAM.tri must be provided')
        end

        % set up parameters for the CALD algorithm. We have to create a
        % structure with the right format. We can then replace default
        % values by values provided by the user
        param.options.vars = {'MeshGridSize', 'MaxSPHARMDegree', 'Tolerance', ...
            'Smoothing', 'Iteration', 'LocalIteration', 't_major', ...
            'SelectDiagonal', 'OutDirectory'};
        param.options.args = [1 1 1 1 1 1 10 10 200];
        param.options.inFilter = {
            '*_bim.mat;*_fix.mat;*_obj.mat' 'Binary and Mesh Objects (*_bim.mat;*_fix.mat;*_obj.mat)'
            '*_obj.mat'                     'Mesh Objects (*_obj.mat)'
            '*_bim.mat;*_fix.mat'           'Binary Objects (*_bim.mat, *_fix.mat)'
            '*.mat'                         'All Objects (*.mat)'
            };
        param.options.default = {'50'  '6'  '2'  '2'  '100'  '10'  ''  ''  './'};
        if (~isfield(param.options, 'MeshGridSize') || isempty(param.options.MeshGridSize))
            param.options.MeshGridSize = 50;
        end
        if (~isfield(param.options, 'MaxSPHARMDegree') || isempty(param.options.MaxSPHARMDegree))
            param.options.MaxSPHARMDegree = 6;
        end
        if (~isfield(param.options, 'Tolerance') || isempty(param.options.Tolerance))
            param.options.Tolerance = 2;
        end
        if (~isfield(param.options, 'Smoothing') || isempty(param.options.Smoothing))
            param.options.Smoothing = 2;
        end
        if (~isfield(param.options, 'Iteration') || isempty(param.options.Iteration))
            param.options.Iteration = 100;
        end
        if (~isfield(param.options, 'LocalIteration') || isempty(param.options.LocalIteration))
            param.options.LocalIteration = 10;
        end
        param.options.t_major = 'x';
        param.options.SelectDiagonal = 'ShortDiag';
        param.options.OutDirectory = '.';
        param.options = orderfields(param.options, ...
            {'vars', 'args', 'inFilter', 'default', ...
            'MeshGridSize', 'MaxSPHARMDegree', 'Tolerance', ...
            'Smoothing', 'Iteration', 'LocalIteration', 't_major', ...
            'SelectDiagonal', 'OutDirectory'});
        
        % initial parametrization
        [sph_verts, name3] ...
            = initParamCALD(x, param.tri, '', param.options);
        
        if (any(isnan(sph_verts)))
            error('Initial parameterization produced NaN values')
        end
        
        % smooth the parameterization
        [~, ~, sph_verts, new_name] ...
            = smootheCALD(x, param.tri, sph_verts, name3, param.options);
        
        % convert xyz coordinates to spherical coordinates
        [lon, lat] ...
            = cart2sph(sph_verts(:, 1), sph_verts(:, 2), sph_verts(:, 3));
        
        % delete files created by the SPHARM functions
        if exist(new_name, 'file')
            delete(new_name) % CALD_smo.mat
        end
        if exist('initParamCALD', 'file')
            rmdir('initParamCALD', 's')
        end
        
        % create output
        uv = [lat lon];
        
%--------------------------------------------------------------------------

    case 'smdscale'

        % check arguments
        nargoutchk(0, 2);

        % number of points to embed
        N = size(x, 1);
        
        % if distance matrix is not provided we need to compute it from the
        % mesh description
        if (~isfield(param, 'd') || isempty(param.d))
            
            % but if the user hasn't provided the mesh triangulation
            % either, we cannot, so we give an error
            if (~isfield(param, 'tri') || isempty(param.tri))
                error('Either PARAM.d or PARAM.tri must be provided')
            end
            
            % if no options are provided for fastmarching, we need to have
            % at least an empty field
            if (~isfield(param, 'options') || isempty(param.options))
                param.options = [];
            end
            
            % default for fast marching method options
            if (~isfield(param, 'dmethod') || isempty(param.dmethod))
                param.dmethod = 'fastmarching';
            end
            if (~isfield(param, 'options') || isempty(param.options))
                param.options = [];
            end
            
            % compute neighbourhood distance matrix from the mesh
            [~, param.d] ...
                = dmatrix_mesh(param.tri, x, param.dmethod, param.options);
            
        end
  
        % defaults
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
                    
                    % if the user didn't provide an estimate for the
                    % parameterization sphere's radius, we can use sphrad
                    if (~isfield(param, 'sphrad') || isempty(param.sphrad))
                        param.sphrad = sphrad;
                    end
                    
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
        
        % if the user has not provided an estimate for the parameterization
        % sphere's radius, we have to compute it internally
        if (~isfield(param, 'sphrad') || isempty(param.sphrad))
            [~, ~, param.sphrad] = proj_on_sphere(x);
        end
        
        % if no options are provided for the MDS algorithm, we need to have
        % at least an empty field
        if (~isfield(param, 'options2') || isempty(param.options2))
            param.options = [];
        end
        
        % compute the parametrization that optimises isometry
        [lat, lon, out.stopCondition, out.err] = ...
            smdscale(param.d, param.sphrad, lat, lon, param.options2);
        
        % create output with the parameterisation values for each point
        % lat = em(1, :)
        % lon = em(2, :)
        uv = [lat lon];
        
%--------------------------------------------------------------------------

    otherwise
        
        error('Parametrization method not implemented')
end
