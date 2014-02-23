% test_tri_qcqp_smacof_nofold_sph.m

%% toy example: triangulation of a sphere

% approximately uniform triangulation of a sphere
[x, tri] = ParticleSampleSphere('N', 10);
N = size(x, 1);
Ntri = size(tri, 1);

% plot sphere
hold off
trisurf(tri, x(:, 1), x(:, 2), x(:, 3))

% compute chord-distance between points
dtot = dmatrix(x');

% connectivity matrix of the mesh
con = dmatrix_mesh(tri);

% local distance matrix
d = dtot;
d(~con) = 0;
d = sparse(d);

% signed volume of the tetrahedra
vs = zeros(Ntri, 1);
for I = 1:Ntri
    idx = tri(I, :);
    vs(I) = tetra_signedvol(x(:), [idx idx+N idx+2*N]);
    if (abs(1/6*det(x(idx, :)) - vs(I)) > 1e-12)
        warning(['tetra_signedvol() failed to compute signed volume value, I = ' num2str(I)])
    end
end
if any(vs <= 0)
    error('The toy example mesh has triangles with negative orientation')
end
min(vs)
max(vs)

% compute constraints
R = 1;
vmin = 0.10;
vmax = 0.15;
[lb, ub, w, A, rl, ru, qc, nl] ...
    = tri_qcqp_smacof_nofold_sph(tri, R, vmin, vmax);

% initialize vertices randomly
%[lon, lat] = cart2sph(x(:, 1), x(:, 2), x(:, 3));
lon0 = (rand(N, 1) - .5) * pi;
lat0 = (rand(N, 1) - .5) * pi/2;
x0 = zeros(N, 3);
[x0(:, 1), x0(:, 2), x0(:, 3)] = sph2cart(lon0, lat0, 1);

% plot algorithm initialization
hold off
trisurf(tri, x0(:, 1), x0(:, 2), x0(:, 3))

% check that the ideal solution agrees with the quadratic constraints
cons_OK = false(1, length(qc.Q));
for I = 1:length(qc.Q)
    cons_OK(I) = (...
        qc.qrl(I) <= (x(:)' * qc.Q{I} * x(:) + x(:)' * qc.l(:, I))) ...
        && ((x(:)' * qc.Q{I} * x(:) + x(:)' * qc.l(:, I)) <= qc.qru(I));
end
if any(~cons_OK)
    warning('Ideal solution does not fulfill quadratic constraints')
end

% SMACOF algorithm parameters
smacof_opts.MaxIter = 100;
smacof_opts.Epsilon = 1e-3;
smacof_opts.Display = 'iter';
smacof_opts.TolFun = 1e-6;

% SCIP algorithm parameters
scip_opts.display = 'off';
scip_opts.display = 'iter';
scip_opts.warnings = 'on';
scip_opts.maxiter = 1500;
scip_opts.maxnodes = 10000;
scip_opts.maxtime = 1000;
scip_opts.gamsfile = 'rcasero-problem.gams';

% solve MDS problem with QCQP-SMACOF with constraints
% we use the mesh with fold-up as starting point
[y, stopCondition, sigma, t] ...
    = qcqp_smacof(d, x0, lb, ub, w, A, rl, ru, qc, nl, ...
    smacof_opts, scip_opts);
