% test_cgal_trifacet_area.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013-2014 University of Oxford
% Version: 0.1.1
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
%% Tetrahedron
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a surface mesh that is a tetrahedron
x = [
    0 0 0
    1 0 0
    0 1 0
    0 0 1
    ];

% plot points
subplot(1, 1, 1)
hold off
plot3(x(:, 1), x(:, 2), x(:, 3), '.')
axis equal

% apply alphashape
[vol, as] = alphavol(x, 1.3);

% plot surface
hold off
trisurf(as.bnd, x(:, 1), x(:, 2), x(:, 3), 'FaceColor', 'blue', 'FaceAlpha', .5)

% compute area of the triangles in the mesh
a = cgal_trifacet_area(as.bnd, x)

% expected result
[.5 sqrt(3)/2 .5 .5]'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Heart ventricles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load segmentation
scimat = scimat_load('data/008-lvhull-downsampled-4.mha');
scimat = scimat_load('data/008-rvhull-downsampled-4.mha');

% compute surface mesh from segmentation
opt = .004;
method = 'simplify';

tic
[x, tri] = v2s(single(scimat.data), 1, opt, method);
tri = tri(:, 1:3);
toc

% plot mesh
subplot(1, 1, 1)
hold off
plotmesh(x, tri)
axis equal

% compute area of the triangles in the mesh
a = cgal_trifacet_area(tri, x);
