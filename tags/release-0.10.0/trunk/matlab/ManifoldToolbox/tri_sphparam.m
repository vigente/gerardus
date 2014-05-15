function [y, yIsValid, stopCondition, sigma, sigma0, t] = tri_sphparam(tri, x, method, d, y0, sphparam_opts, smacof_opts, scip_opts)
% TRI_SPHPARAM  Spherical parametrization of closed triangular mesh.
%
% [Y, YISVALID, STOPCONDITION, SIGMA, SIGMA0, T] = tri_sphparam(TRI, X, METHOD)
%
%   TRI is a 3-column matrix with a surface mesh triangulation. Each row
%   gives the indices of one triangle. The mesh needs to be a 2D manifold,
%   that can be embedded in 2D or 3D space. The orientation of the
%   triangles does not matter.
%
%   X is a 3-column matrix with the coordinates of the mesh vertices. Each
%   row gives the (x,y,z)-coordinates of one vertex.
%
%   METHOD is a string that selects the parametrization method:
%
%     'cmdscale': Classical Multidimensional Scaling (MDS).
%
%     'smacof':   Unconstrained SMACOF, followed by projection of points on
%                 sphere.
%
%     'consmacof-local': Constrained SMACOF with local untangling of
%                 connected vertices.
%
%     'consmacof-global': Constrained SMACOF with untangling of all
%                 vertices simulataneously (too slow except for very small
%                 problems).
%
%   Y is a 3-column matrix with the coordinates of the spherical
%   parametrization of the mesh. Each row contains the (x,y,z)-coordinates
%   of a point on the sphere. Y is the best valid solution found by the
%   algorithm within all iterations (including the intial guess). If no
%   valid solution could be found, the initial guess is returned. To
%   compute the spherical coordinates of the points, run
%
%     [lon, lat, r] = cart2sph(y(:, 1), y(:, 2), y(:, 3));
%
%   YISVALID is a boolean flag that indicates whether a valid solution
%   could be found and returned.
%
%   STOPCONDITION is a cell array with a string for each stop condition
%   that made the algorithm stop at the last iteration.
%
%   SIGMA is a vector with the stress value at each iteration. A value of
%   NaN means that the algorithm could not find a valid solution in the
%   corresponding iteration.
%
%   SIGMA0 is a scalar with the stress value of the initial guess Y0. Note
%   that we have no guarantee that Y0 is a valid solution, so be careful
%   when comparing SIGMA0 to SIGMA. In particular, SIGMA0 may be smaller
%   than any SIGMA, but may not be a valid solution according to the
%   constraints.
%
%   T is a vector with the time between the beginning of the algorithm and
%   each iteration. Units in seconds.
%
%   In 'consmacof-local', STOPCONDITION, SIGMA and T are cell arrays, with
%   the output parameters for each connected component untangled by the
%   algorithm.
%
%
% ... = tri_sphparam(..., D, Y0, SPHPARAM_OPTS, SMACOF_OPTS, SCIP_CONS)
%
%   D is the square distance matrix. For 'cmdscale', D must be a full
%   matrix. For the other methods, it can be a sparse or full matrix.
%   D(i,j)=0 means that the distance between vertices i and j is not
%   considered for the stress measure.
%
%   Y0 is an initial guess for the output parametrization. For 'cmdscale',
%   Y0 must be empty. For SMACOF methods, the choice of Y0 is important
%   because the algorithm can be trapped into local minima.
%
%   SPHPARAM_OPTS is a struct with parameters to tweak the spherical
%   parametrization algorithm.
%
%     'Display': (default = 'off') Do not display any internal information.
%                'iter': display internal information at every iteration.
%
%     'TopologyCheck': (default false) Check that parametrization has no
%                self-intersections and that all triangles have a positive
%                orientation.
%
%     'volmin':  (default 0) Only used by constrained SMACOF methods.
%                Minimum volume allowed to the oriented spherical
%                tetrahedra at the output, formed by the triangles and the
%                centre of the sphere. Note that if volmin>0, then all
%                output triangles have outwards-pointing normals.
%
%     'volmax':  (default Inf) Only used by constrained SMACOF methods.
%                Maximum volume of the output tetrahedra (see 'volmin').
%
%     'LocalConvexHull': (default true) Only used by 'consmacof-local' 
%                method. When a local neighbourhood is tangled, untangle
%                the convex hull that contains it. The reason is that
%                untangling a convex domain is simpler, and it produces
%                better quality solutions.
%
%   SMACOF_OPTS is a struct with parameters to tweak the SMACOF algorithm.
%   See cons_smacof_pip for details.
%
%   SCIP_OPTS is a struct with parameters to tweak the SCIP algorithm. See
%   cons_smacof_pip for details.
%
% See also: cmdscale, cons_smacof_pip, qcqp_smacof.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.3.0
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

%% Process input to the function

% check arguments
narginchk(3, 8);
nargoutchk(0, 6);

% start clock
tic;

% common defaults
if (nargin < 4 || isempty(d))
    % by default, we compute the full distance matrix between vertices in
    % the mesh using Fast Marching, as a linear approximation to geodesic
    % distances on the manifold
    [~, d] = dmatrix_mesh(tri, x, 'fastmarching');
end
if (nargin < 6)
    sphparam_opts = [];
end    
if (nargin < 7)
    smacof_opts = [];
end
if (nargin < 8)
    scip_opts = [];
end

% number of vertices
N = size(x, 1);

% sphparam_opts defaults
if (~isfield(sphparam_opts, 'sphrad'))
    % estimate the output parametrization sphere's radius, such that the
    % sphere's surface is the same as the total surface of the mesh
    sphparam_opts.sphrad = estimate_sphere_radius(tri, x);
end
if (~isfield(sphparam_opts, 'volmin'))
    sphparam_opts.volmin = 0;
end
if (~isfield(sphparam_opts, 'volmax'))
    sphparam_opts.volmax = Inf;
end
if (~isfield(sphparam_opts, 'Display'))
    sphparam_opts.Display = 'none';
end
if (~isfield(sphparam_opts, 'LocalConvexHull'))
    sphparam_opts.LocalConvexHull = true;
end
if (~isfield(sphparam_opts, 'TopologyCheck'))
    sphparam_opts.TopologyCheck = false;
end

% smacof_opts defaults
if (~isfield(smacof_opts, 'MaxIter'))
    smacof_opts.MaxIter = 100;
end
if (~isfield(smacof_opts, 'Epsilon'))
    smacof_opts.Epsilon = 1e-4;
end
if (~isfield(smacof_opts, 'Display'))
    smacof_opts.Display = 'none';
end
if (~isfield(smacof_opts, 'TolFun'))
    smacof_opts.TolFun = 1e-6;
end

% scip_opts defaults
if (~isfield(scip_opts, 'limits_solutions'))
    % from SCIP we only need that it enforces the constraints, and we let
    % SCMACOF optimize the stress
    scip_opts.limits_solutions = 1;
end
if (~isfield(scip_opts, 'display_verblevel'))
    % by default, be silent
    scip_opts.display_verblevel = 0;
end

% if Y0 not provided, we make it empty. This simplifies a bit the code
% after
if (nargin < 5)
    y0 = [];
end

if (strcmp(method, 'cmdscale'))
    
    % enforce that the classical MDS method cannot have an initial guess
    if (~isempty(y0))
        error('With CMDSCALE method, initial guess cannot be provided')
    end
    
else
    
    if (isempty(y0))
        
        % in non-classical MDS methods, if no initial guess is provided,
        % compute a random sampling of the sphere
        y0 = rand(N, 3);
        y0 = y0 ./ repmat(sqrt(sum(y0.^2, 2)), 1, 3) * sphparam_opts.sphrad;
        
    else
        
        % if initial guess is provided, make sure that it corresponds to
        % points on the sphere surface
        if any(abs(sqrt(sum(y0.^2, 2)) - sphparam_opts.sphrad) > 1e-10)
            error(['Initial guess points are not on a sphere of radius ' num2str(sphparam_opts.sphrad)])
        end
        
    end
    
end

% check whether the initial guess is a valid solution
if isempty(y0)
    y0IsValid = [];
else
    y0IsValid = all(sphtri_signed_vol(tri, y0) > 0);
end

% check inputs dimensions
if (size(tri, 2) ~= 3)
    error('TRI must have 3 columns')
end
if (size(x, 2) ~= 3)
    error('X must have 3 columns')
end
if ((N ~= size(d, 1)) || (N ~= size(d, 2)))
    error('D must be a square matrix with the same number of rows as X')
end

%% Compute output parametrization with one of the implemented methods

if (strcmp(sphparam_opts.Display, 'iter'))
    fprintf('Parametrization method: %s\n', method)
end

% all methods implemented here operate with Euclidean chord distances,
% instead of the geodesic distances on the surface of the sphere
d = arclen2chord(d, sphparam_opts.sphrad);

switch method
    
    %% Classic Multidimensional Scaling (MDS)
    case 'cmdscale'

        narginchk(3, 6);
        
        % Classical MDS requires a full distance matrix
        if (issparse(d))
            error('Classical MDS does not accept sparse distance matrices')
        end
        
        % classical MDS parametrization. This will produce something
        % similar to a sphere, if the d matrix is not too far from being
        % Euclidean
        y = cmdscale(d);
        y = y(:, 1:3);
        
        % project the MDS solution on the sphere
        [lat, lon] = proj_on_sphere(y);
        [y(:, 1), y(:, 2), y(:, 3)] ...
            = sph2cart(lon, lat, sphparam_opts.sphrad);
        
        % signed volume of tetrahedra formed by sphere triangles and origin
        % of coordinates
        vol = sphtri_signed_vol(tri, y);
        
        % if more than half triangles have negative areas, we mirror the
        % parametrization (because MDS is invariant to "inside-out" sphere
        % transformations)
        if (nnz(vol<0) > length(vol)/2)
            [lon, lat, sphrad] = cart2sph(y(:, 1), y(:, 2), y(:, 3));
            [y(:, 1), y(:, 2), y(:, 3)] = sph2cart(-lon, lat, sphrad);
        end
        
        % Classic MDS always produces a global optimum
        stopCondition = 'Global optimum';
        
        % stress of output parametrization
        sigma = sum(sum((d - dmatrix(y')).^2));
        
        % initial guess is empty, so its stress is empty too
        sigma0 = [];
       
        % time for initial parametrization
        t = toc;
        
    %% SMACOF algorithm ("Scaling by majorizing a convex function"
    case 'smacof'
        
        narginchk(3, 7);
        
        % compute SMACOF parametrization
        [y, stopCondition, sigma, t] = smacof(d, y, [], smacof_opts);
    
        % project the SMACOF solution on the sphere, and recompute the
        % sphere radius as the median of the radii of all vertices
        % projected on the sphere
        [lat, lon, sphparam_opts.sphrad] = proj_on_sphere(y);
        [y(:, 1), y(:, 2), y(:, 3)] ...
            = sph2cart(lon, lat, sphparam_opts.sphrad);
        sigma0 = [];
        
    %% Constrained SMACOF, local optimization
    case 'consmacof-local'
        
        %% Find tangled vertices and group in clusters of connected tangled vertices
        
        % spherical coordinates of the points
        [lon, lat] = cart2sph(y(:, 1), y(:, 2), y(:, 3));
        
        % reorient all triangles so that all normals point outwards
        [~, tri] = meshcheckrepair(x, tri, 'deep');
        
        % signed volume of tetrahedra formed by sphere triangles and origin
        % of coordinates
        vol = sphtri_signed_vol(tri, y);
        
        % we mark as tangled all vertices from tetrahedra with negative
        % volumes, because they correspond to triangles with normals
        % pointing inwards
        isFree = false(N, 1);
        isFree(unique(tri(vol <= 0, :))) = true;
        
        % find triangles that cause self-intersections
        idx = cgal_check_self_intersect(tri, y);
        
        % we also mark as tangled all vertices from triangles that cause
        % self-intersections in the mesh
        isFree(unique(tri(idx>0, :))) = true;
        
        % mesh connectivity matrix
        dcon = dmatrix_mesh(tri);
        
        % find groups of connected tangled vertices
        [Ncomp, cc] = graphcc(dcon(isFree, isFree));
        
        % the vertices in cc refer to the smaller dcon(isFree, isFree)
        % matrix. We need to rename them so that they refer to the full
        % matrix dcon(:, :)
        map = find(isFree)';
        if (~isempty(map))
            cc = cellfun(@(x) map(x), cc, 'UniformOutput', false);
        end
        
        % initialize outputs, with one element per component
        stopCondition = cell(1, Ncomp);
        sigma = cell(1, Ncomp);
        sigma0 = zeros(1, Ncomp);
        t = cell(1, Ncomp);
        
        %% Untangle parametrization: untangle clusters of tangled vertices, one by
        %% one
        
        % untangle each component separately
        for C = 1:Ncomp
            
            if (strcmp(sphparam_opts.Display, 'iter'))
                fprintf('** Untangling component %d/%d\n', C, Ncomp)
            end
            
            % start the local neighbourhood with the free vertices
            isFreenn = false(N, 1);
            isFreenn(cc{C}) = true;
            
            % add to the local neighbourhood all the neighbours of the free
            % vertices. Note that the neighbours must be fixed, because if
            % they were free, they would have been included in the
            % connected component by graphcc()
            nn = full(isFreenn' | (sum(dcon(isFreenn, :), 1) > 0))';
            
            % Local Convex Hull block:
            % This code snippet converts the local neighbourhood to the
            % convex hull of the local neighbourhood. Untangling a convex
            % local neighbourhood should be easier than a non-convex, but
            % it also involves more vertices
            if (sphparam_opts.LocalConvexHull)
                
                % mean point of the local neighbourhood
                [latnnm, lonnnm] = meanm(lat(nn), lon(nn), 'radians');
                ynnm = [0 0 0]';
                [ynnm(1), ynnm(2), ynnm(3)] = sph2cart(lonnnm, latnnm, sphparam_opts.sphrad);
                
                % rotation matrix to take the centroid to lat=0, lon=0
                rot = vrrotvec2mat([cross(ynnm, [1 0 0]'); ...
                    acos(dot(ynnm/norm(ynnm), [1 0 0]'))]);
                
                % rotate all vertices so that the local neighbourhood is centered
                % around (0,0)
                yrot = (rot * y')';
                
                % convert to spherical coordinates
                [lonrot, latrot] = cart2sph(yrot(:, 1), yrot(:, 2), yrot(:, 3));
                
                % update the local neighbourhood so that the local neighbourhood is
                % the convex hull
                nn = nn | inhull([lonrot, latrot], [lonrot(nn), latrot(nn)]);
                
            end
            
            % triangles that triangulate the local neighbourhood
            idxtrinn = sum(ismember(tri, find(nn)), 2) == 3;
            trinn = tri(idxtrinn, :);
            
            % at this point, it's possible that the local triangulation
            % doesn't contain all the vertices in the local neighbourhood.
            % Thus, we drop isolated vertices that don't have an associated
            % triangle
            nn(:) = false;
            nn(unique(trinn)) = true;
            
            % boundary of the triangulation. We are looking for edges that
            % appear only once in the triangulation. Those edges form the
            % boundary.
            edgenn = sort([trinn(:, 1:2); trinn(:, 2:3); trinn(:, [3 1])], 2);
            [aux, ~, idx] = unique(edgenn, 'rows');
            idx = hist(idx, 1:max(idx));
            vedgenn = unique(aux(idx == 1, :));
            
            % to speed things up, we want to pass to SMACOF a subproblem
            % created only from the local neighbourhood. Here, we create
            % the local neighbourhood variables for convenience
            isFreenn = true(N, 1);
            isFreenn(vedgenn) = false;
            isFreenn = isFreenn(nn);
            [trinn, ynn] = tri_squeeze(trinn, y);
            dnn = d(nn, nn);
            
            % recompute bounds and constraints for the spherical problem
            [con, bnd] ...
                = tri_ccqp_smacof_nofold_sph_pip(trinn, ...
                sphparam_opts.sphrad, sphparam_opts.volmin, ...
                sphparam_opts.volmax, isFreenn, ynn);
            
            % solve MDS problem with constrained SMACOF
            [y(nn, :), stopCondition{C}, sigma{C}, sigma0(C), t{C}] ...
                = cons_smacof_pip(dnn, ynn, isFreenn, bnd, [], con, ...
                smacof_opts, scip_opts);
            
            if (sphparam_opts.TopologyCheck)
                
                % assertion check: after untangling, the local neighbourhood cannot
                % produce self-intersections
                if any(cgal_check_self_intersect(trinn, y(nn,:)))
                    warning(['Component ' num2str(C) ...
                        ' contains self-intersections after untangling'])
                end
                
                % assertion check: after untangling, volumes of all tetrahedra in the
                % local neighbourhood must be positive
                aux = sphtri_signed_vol(trinn,  y(nn, :));
                if any(aux < sphparam_opts.volmin | aux > sphparam_opts.volmax)
                    warning(['Component ' num2str(C) ...
                        ' contains tetrahedra with volumes outside the constraint values'])
                end
                
            end
            
            % update spherical coordinates of new points
            [lon(nn), lat(nn)] = cart2sph(y(nn, 1), y(nn, 2), y(nn, 3));
            
            if (strcmp(sphparam_opts.Display, 'iter'))
                fprintf('... Component %d/%d done. Time: %.4e\n', C, Ncomp, toc)
                fprintf('===================================================\n')
            end
            
        end
        
    %% Constrained SMACOF, global optimization
    case 'consmacof-global'
        
        %% Find tangled vertices
        
        % reorient all triangles so that all normals point outwards
        [~, tri] = meshcheckrepair(x, tri, 'deep');
        
        % signed volume of tetrahedra formed by sphere triangles and origin
        % of coordinates
        vol = sphtri_signed_vol(tri, y0);
        
        % we mark as tangled all vertices from tetrahedra with negative
        % volumes, because they correspond to triangles with normals
        % pointing inwards
        isFree = false(N, 1);
        isFree(unique(tri(vol <= 0, :))) = true;
        
        % find triangles that cause self-intersections
        idx = cgal_check_self_intersect(tri, y0);
        
        % we also mark as tangled all vertices from triangles that cause
        % self-intersections in the mesh
        isFree(unique(tri(idx>0, :))) = true;
        
        %% Untangle parametrization
        
        % recompute bounds and constraints for the spherical problem
        [con, bnd] ...
            = tri_ccqp_smacof_nofold_sph_pip(tri, ...
            sphparam_opts.sphrad, sphparam_opts.volmin, ...
            sphparam_opts.volmax, isFree, y0);
        
        % solve MDS problem with constrained SMACOF
        [y, stopCondition, sigma, sigma0, t] ...
            = cons_smacof_pip(d, y0, isFree, bnd, [], con, ...
            smacof_opts, scip_opts);
            
    otherwise
        error(['Unknown parametrization method: ' method])
end

%% post-processing of the solution

% check whether found solution is valid
yIsValid = all(sphtri_signed_vol(tri, y) > 0);

% if the initial guess was better than the found solution, we have to give
% it as the output
switch method
    case 'cmdscale'
        
        % nothing to do
        
    otherwise
    
        % the best valid solution we have found has this stress
        sigmabest = min(sigma);

        % we need stopCondition to be a cell array, so that we can add more
        % strings to it
        if ~iscell(stopCondition)
            stopCondition = {stopCondition};
        end

        % possible cases comparing initial guess to solutions obtained by algorithm
        
        % both initial guess and best found solution are valid
        if (y0IsValid && yIsValid)
            
            % if the best solution is worse than the initial guess, we return the
            % initial guess
            if (sigma0 <= sigmabest)
                y = y0;
                stopCondition{end+1} = 'Returned solution is initial guess';
            else
                stopCondition{end+1} = 'Returned solution is best solution found by algorithm';
            end
            
            % only the initial guess was valid, but best found solution isn't
        elseif (y0IsValid)
            
            % we return the initial guess
            y = y0;
            yIsValid = true;
            stopCondition{end+1} = 'Returned solution is initial guess';
            
            % the best found solution is valid, but the initial guess wasn't
        elseif (yIsValid)
            
            % we don't need to do anything. The best found solution will be
            % returned
            stopCondition{end+1} = 'Returned solution is best solution found by algorithm';
            
            % neither initial guess was valid nor valid solution was found
        else
            
            % return array of NaNs
            y = nan(size(y));
            stopCondition{end+1} = 'No valid solution found';
            
        end
        
end

if (strcmp(sphparam_opts.Display, 'iter'))
    if (yIsValid)
        fprintf('\t Returned parametrization is valid\n')
    else
        fprintf('\t Returned parametrization is invalid\n')
    end
    fprintf('... Parametrization done. Total time: %.4e (sec)\n', toc)
end

% DEBUG: plot parametrization solution
% hold off
% trisurf(tri, y(:, 1), y(:, 2), y(:, 3))
% axis equal


%% assertion check for self-intersections or negative tetrahedra

if (sphparam_opts.TopologyCheck)
    
    if (strcmp(sphparam_opts.Display, 'iter'))
        fprintf('Checking output parametrization tolopogy\n')
    end
    
    % assertion check: after untangling, the local neighbourhood cannot
    % produce self-intersections
    if (yIsValid && any(cgal_check_self_intersect(tri, y)))
        warning('Assertion fail: Output parametrization has only positive triangles but mesh self-intersects')
    else
        disp('Output parametrization has only positive triangles and mesh does not self-intersect')
    end
    
    if (strcmp(sphparam_opts.Display, 'iter'))
        fprintf('... done checking output parametrization tolopogy\n')
    end
    
end

end


%% Auxiliary functions

% Estimate of the sphere's radius such that the sphere's surface is the
% same as the total surface of the mesh
function sphrad = estimate_sphere_radius(tri, x)

% compute area of all mesh triangles
a = cgal_trifacet_area(tri, x);
atot = sum(a);

% estimate radius of parameterization sphere
sphrad = sqrt(atot/4/pi);

end
