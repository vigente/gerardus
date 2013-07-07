% test_congrow.m

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

%% small triangular mesh of a rectangular grid
[tri, uv] = surface_tridomain('rect', 1:3, 1:2);

% connectivity matrix
con = dmatrix_mesh(tri);

% plot mesh
hold off
gplot(con, uv)
axis ij
axis([.5 3.5 .5 2.5])

% grow local neighbourhood to increase connectivity
con2 = congrow(con);

% plot new connections
hold on
gplot(con2, uv, 'r')

% compute distances of new neighbourhood
d = dmatrix_con(con2, uv);
full(d)

% remove too long connections for local neighbourhood
d(d > sqrt(2)) = 0;

% plot corrected local neighbourhood
hold off
gplot(con, uv)
hold on
gplot(d, uv, 'r')
axis ij
axis([.5 3.5 .5 2.5])

%% bigger triangular mesh of a rectangular grid
[tri, uv] = surface_tridomain('rect', 1:10, 1:7);

% connectivity matrix
con = dmatrix_mesh(tri);

% plot mesh
hold off
gplot(con, uv)
axis ij
axis([.5 10.5 .5 7.5])

% grow local neighbourhood to increase connectivity
con2 = congrow(con);

% plot new connections
hold on
gplot(con2, uv, 'r')

% compute distances of new neighbourhood
d = dmatrix_con(con2, uv);

% remove too long connections for local neighbourhood
d(d > sqrt(2)) = 0;

% plot corrected local neighbourhood
hold off
gplot(con, uv)
hold on
gplot(d, uv, 'r')
axis ij
axis([.5 10.5 .5 7.5])

%% spherical mesh

% Uniformly distribute particles across the surface of the unit sphere 
[xyz,tri,~,Ue]=ParticleSampleSphere('N', 10); % this operation takes ~8 sec on my machine (6GB RAM & 2.8 GHz processor)

% plot mesh
hold off
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3), 'FaceColor', 'green')

% connectivity matrix
con = dmatrix_mesh(tri);

% grow local neighbourhood to increase connectivity
con2 = congrow(con);

% plot new connections
hold off
gplot3d(con2, xyz, 'r')
hold on
gplot3d(con, xyz, 'b')

