% test_lmdscale.m

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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check that result for full distance matrix is the same as for classical
%% MDS

% quadrangular mesh split into asymmetric triangles

% create rectangle with triangular mesh
[tri, x] = surface_tridomain('rect', 'step', 1/10, [0 0], [1 1/5*3]);
x(:, end+1) = 0;

% number of vertices in the mesh
N = size(x, 1);

% compute full distance matrix from the mesh
d = zeros(N);
for I = 1:N
    d(:, I) = perform_fast_marching_mesh(x, tri, I);
end

% parametrization using classical MDS
uv0 = cmdscale(d);

% the output parametrization is 2D
uv0 = uv0(:, 1:2);

% register solution to input mesh
[~, uv0] = procrustes(x(:, 1:2), uv0, 'Scaling', false);

% plot original mesh and parametrization
con = dmatrix_mesh(tri);
hold off
gplot(con, x(:, 1:2))
hold on
gplot(con, uv0, 'r')

% parametrization using local MDS
opt.MaxIter = 20;
opt.MaxInc = .002;
[u, v, stopCondition, err, dout] = lmdscale(d, [], [], opt);
uv = [u v];

% register solution to input mesh
[~, uv] = procrustes(x(:, 1:2), uv, 'Scaling', false);

% plot error evolution
subplot(2, 2, 1)
hold off
plot(err.rawstress)
title('Raw stress')
hold off
subplot(2, 2, 2)
plot(err.stress1)
title('Stress-1')
subplot(2, 1, 2)
hold off
plot(err.maxinc)
hold on
plot(err.medinc)
title('Location displacement')

% plot original mesh and parametrizations
subplot(1, 1, 1)
hold off
gplot(con, x(:, 1:2))
hold on
gplot(con, uv0, 'r')
gplot(con, uv, 'g')
