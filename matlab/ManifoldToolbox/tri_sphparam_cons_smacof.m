function [y, t] = tri_sphparam_cons_smacof(tri, x, d, sphparam_opts, smacof_opts, scip_opts)
% TRI_SPHPARAM_CONS_SMACOF Pseudo-isometric parametrization of closed
% triangular mesh using constrained SMACOF.
%
% [Y, T] = tri_sphparam_cons_smacof(TRI, X)
%
%   TRI is a 3-column matrix with a surface mesh triangulation. Each row
%   gives the indices of one triangle. The mesh needs to be a 2D manifold,
%   that can be embedded in 2D or 3D space. The orientation of the
%   triangles does not matter.
%
%   X is a 3-column matrix with the coordinates of the mesh vertices. Each
%   row gives the (x,y,z)-coordinates of one vertex.
%
%   Y is a 3-column matrix with the coordinates of the spherical
%   parametrization of the mesh. Each row contains the (x,y,z)-coordinates
%   of a point on the sphere. To compute the spherical coordinates of the
%   points, run
%
%     [lon, lat, r] = cart2sph(y(:, 1), y(:, 2), y(:, 3));
%
%   T is a vector with the accumulated time at each step of the algorithm.
%   T(1) is the time until the initial MDS parametrization is computed.
%   T(2:NCC) are the times until each of the NCC tangled components is
%   untangled. If the final refinement is selected, T(end) is the time
%   until the final refinement is computed.
%
% ... = tri_sphparam_cons_smacof(..., CONS_OPTS, SMACOF_OPTS, SCIP_CONS)
%
%   SPHPARAM_OPTS is a struct with parameters to tweak the spherical
%   parametrization algorithm.
%
%     'volmin':  (default 0) Minimum volume allowed to the oriented
%                spherical tetrahedra at the output, formed by the
%                triangles and the centre of the sphere. Note that if
%                volmin>0, then all output triangles have outwards-pointing
%                normals.
%
%     'volmax':  (default Inf) Maximum volume of the output tetrahedra (see
%                'volmin').
%
%     'Display': (default = 'off') Do not display any internal information.
%                'iter': display internal information at every iteration.
%
%     'LocalConvexHull': (default true) When a local neighbourhood is
%                tangled, untangle the convex hull that contains it. The
%                reason is that untangling a convex domain is simpler, and
%                it produces better quality solutions.
%
%     'FinalRefinement': (default = false) Run constrained SMACOF on the
%                whole mesh after the initial untangling, to try to further
%                reduce the stress measure. This is a very slow process.
%
%   SMACOF_OPTS is a struct with parameters to tweak the SMACOF algorithm.
%   See cons_smacof_pip for details.
%
%   SCIP_OPTS is a struct with parameters to tweak the SCIP algorithm. See
%   cons_smacof_pip for details.
%
% See also: cons_smacof_pip.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.2.0
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
narginchk(3, 6);
nargoutchk(0, 2);

% start clock
tic;

% defaults
if (nargin < 3 || isempty(d))
    
    % by default, we compute the full distance matrix between vertices in
    % the mesh using Fast Marching
    [~, d] = dmatrix_mesh(tri, x, 'fastmarching');
    
end
if (nargin < 4 || isempty(sphparam_opts))
    
    sphparam_opts.volmin = 0;
    sphparam_opts.volmax = Inf;
    sphparam_opts.Display = 'none';
    sphparam_opts.LocalConvexHull = true;
    sphparam_opts.FinalRefinement = false;
    
end    
if (nargin < 5 || isempty(smacof_opts))
    
    smacof_opts.MaxIter = 100;
    smacof_opts.Epsilon = 1e-4;
    smacof_opts.Display = 'none';
    smacof_opts.TolFun = 1e-6;
    
end
if (nargin < 6 || isempty(scip_opts))

    % from SCIP we only need that it enforces the constraints, and we let
    % SCMACOF optimize the stress
    scip_opts.limits_solutions = 1;
    
    % by default, be silent
    scip_opts.display_verblevel = 0;
    
end

% number of vertices and triangles
N = size(x, 1);
Ntri = size(tri, 1);

% check inputs
if (size(tri, 2) ~= 3)
    error('TRI must have 3 columns')
end
if (size(x, 2) ~= 3)
    error('X must have 3 columns')
end
if ((N ~= size(d, 1)) || (N ~= size(d, 2)))
    error('D must be a square matrix with the same number of rows as X')
end

%% Parametrization initialization: classical MDS

if (strcmp(sphparam_opts.Display, 'iter'))
    fprintf('===================================================\n')
    fprintf('Computing initial MDS parametrization\n')
end

% compute area of all mesh triangles
a = cgal_trifacet_area(tri, x);
atot = sum(a);

% estimate radius of parameterization sphere
sphrad = sqrt(atot/4/pi);

% the distances in the dtot matrix are Fast Marching distances over the
% mesh. We are going to assume these are the arc length distances on
% the sphere. We want to convert them to angle distances, and then to
% chord distances for the classical MDS initial solution
dchord = d / sphrad; % angular distances
dchord = 2 * sphrad * sin(dchord/2); % chord distances

% for the initial parameterization we can use classical MDS
y = cmdscale(dchord);
y = y(:, 1:3);

% project the MDS solution on the sphere, and recompute the sphere radius.
% The sphere radius is the median of the radii of all vertices projected on
% the sphere.
% Using classical MDS can create fold-overs in the parametrization
[lat, lon, sphrad] = proj_on_sphere(y);
[y(:, 1), y(:, 2), y(:, 3)] = sph2cart(lon, lat, sphrad);

% reorient all triangles so that all normals point outwards
[~, tri] = meshcheckrepair(y, tri, 'deep');

% recompute chord distances with the new sphere radius
dchord = d / sphrad; % angular distances
dchord = 2 * sphrad * sin(dchord/2); % chord distances

% time for initial parametrization
t = toc;

if (strcmp(sphparam_opts.Display, 'iter'))
    fprintf('... Initial MDS parametrization done. Time: %.4e (sec)\n', t)
    fprintf('===================================================\n')
end
    
% DEBUG: plot classic MDS solution
% hold off
% trisurf(tri, y(:, 1), y(:, 2), y(:, 3))
% axis equal

%% Find tangled vertices and group in clusters of connected tangled vertices

% signed volume of tetrahedra formed by sphere triangles and origin of
% coordinates
vol = sphtri_signed_vol(tri, y);

% we mark as tangled all vertices from tetrahedra with negative volumes,
% because they correspond to triangles with normals pointing inwards
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

% the vertices in cc refer to the smaller dcon(isFree, isFree) matrix. We
% need to rename them so that they refer to the full matrix dcon(:, :)
map = find(isFree)';
cc = cellfun(@(x) map(x), cc, 'UniformOutput', false);

% resize the time vector so that we can add times taken for each component
t(end+1:end+Ncomp) = zeros(1, Ncomp);

%% Untangle parametrization: untangle clusters of tangled vertices, one by
%% one

% untangle each component separately
for C = 1:Ncomp
    
    if (strcmp(sphparam_opts.Display, 'iter'))
        fprintf('Untangling component %d/%d\n', C, Ncomp)
    end
    
    % start the local neighbourhood with the free vertices
    isFreenn = false(N, 1);
    isFreenn(cc{C}) = true;
    
    % add to the local neighbourhood all the neighbours of the free
    % vertices. Note that the neighbours must be fixed, because if they
    % were free, they would have been included in the connected component
    % by graphcc()
    nn = full(isFreenn' | (sum(dcon(isFreenn, :), 1) > 0))';
    
    % Local Convex Hull block:
    % This code snippet converts the local neighbourhood to the convex hull
    % of the local neighbourhood. Untangling a convex local neighbourhood
    % should be easier than a non-convex, but it also involves more
    % vertices
    if (sphparam_opts.FinalRefinement)
        
        % mean point of the local neighbourhood
        [latnnm, lonnnm] = meanm(lat(nn), lon(nn), 'radians');
        ynnm = [0 0 0]';
        [ynnm(1), ynnm(2), ynnm(3)] = sph2cart(lonnnm, latnnm, sphrad);
        
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
    
    % at this point, it's possible that the local triangulation doesn't
    % contain all the vertices in the local neighbourhood. Thus, we drop
    % isolated vertices that don't have an associated triangle
    nn(:) = false;
    nn(unique(trinn)) = true;
    
    % boundary of the triangulation. We are looking for edges that
    % appear only once in the triangulation. Those edges form the
    % boundary.
    edgenn = sort([trinn(:, 1:2); trinn(:, 2:3); trinn(:, [3 1])], 2);
    [aux, ~, idx] = unique(edgenn, 'rows');
    idx = hist(idx, 1:max(idx));
    vedgenn = unique(aux(idx == 1, :));
    
    % to speed things up, we want to pass to SMACOF a subproblem created
    % only from the local neighbourhood. Here, we create the local
    % neighbourhood variables for convenience
    isFreenn = true(N, 1);
    isFreenn(vedgenn) = false;
    isFreenn = isFreenn(nn);
    [trinn, ynn] = tri_squeeze(trinn, y);
    dchordnn = dchord(nn, nn);
    
    % recompute bounds and constraints for the spherical problem
    [con, bnd] = tri_ccqp_smacof_nofold_sph_pip(trinn, sphrad, ...
        sphparam_opts.volmin, sphparam_opts.volmax, isFreenn, ynn);
    
    % solve MDS problem with constrained SMACOF
    y(nn, :) = cons_smacof_pip(dchordnn, ynn, isFreenn, ...
        bnd, [], con, smacof_opts, scip_opts);
    
    % assertion check: after untangling, the local neighbourhood cannot
    % produce self-intersections
    if any(cgal_check_self_intersect(trinn, y(nn,:)))
        warning(['Assertion fail: Component ' num2str(C) ...
            ' contains self-intersections after untangling'])
    end
    
    % assertion check: after untangling, volumes of all tetrahedra in the
    % local neighbourhood must be positive
    aux = sphtri_signed_vol(trinn,  y(nn, :));
    if any(aux < sphparam_opts.volmin | aux > sphparam_opts.volmax)
        warning(['Assertion fail: Component ' num2str(C) ...
            ' contains tetrahedra with volumes outside the constraint values'])
    end
    
    % update spherical coordinates of new points
    [lon(nn), lat(nn)] = cart2sph(y(nn, 1), y(nn, 2), y(nn, 3));
    
    % update time spent until now
    t(C+1) = toc;
    
    if (strcmp(sphparam_opts.Display, 'iter'))
        fprintf('... Component %d/%d done. Time: %.4e\n', C, Ncomp, t(C+1))
        fprintf('===================================================\n')
    end
    
end

%% Final refinement: starting from a solution that fulfils the constraints,
%% decrease stress if possible

if (sphparam_opts.FinalRefinement)
    
    if (strcmp(sphparam_opts.Display, 'iter'))
        fprintf('===================================================\n')
        fprintf('Final refinement\n')
    end
    
%     % recompute bounds and constraints for the spherical problem
%     [con, bnd] = tri_ccqp_smacof_nofold_sph_pip(tri, sphrad, ...
%         sphparam_opts.volmin, sphparam_opts.volmax, [], y);
%     
%     % solve MDS problem with constrained SMACOF
%     y = cons_smacof_pip(dchord, y, [], ...
%         bnd, [], con, smacof_opts, scip_opts);
    
    % time to compute the refinement
    t(end+1) = toc;
    
    if (strcmp(sphparam_opts.Display, 'iter'))
        fprintf('... Final refinement done. Time: %.4e\n', t(end))
        fprintf('===================================================\n')
    end
    
end

%% assertion check for self-intersections or negative tetrahedra

if (strcmp(sphparam_opts.Display, 'iter'))
    fprintf('===================================================\n')
    fprintf('Checking that the output spherical mesh has no tangled vertices\n')
end

% assertion check: after untangling, the local neighbourhood cannot
% produce self-intersections
if any(cgal_check_self_intersect(tri, y))
    warning('Assertion fail: Mesh contains self-intersections after untangling')
end

% assertion check: after untangling, volumes of all tetrahedra in the
% local neighbourhood must be positive
aux = sphtri_signed_vol(tri,  y);
if any(aux < sphparam_opts.volmin | aux > sphparam_opts.volmax)
    warning('Assertion fail: Mesh contains tetrahedra with volumes outside the constraint values')
end
