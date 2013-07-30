% test_cgal_closest_trifacet.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
% Version: 0.2.2
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

xi = [
    .6 .6 .6
    .25 .25 .35
    .10 .40 .09
    -.25 .5 .25
    ];

% plot points
hold on
plot3(xi(:, 1), xi(:, 2), xi(:, 3), 'r*')


% find closest facet and distance to a point
[f, d, p] = cgal_closest_trifacet(as.bnd, x, xi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Right ventricle mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load alpha-shape convex hull binary mask
aux = scinrrd_load('data/008-rvhull-downsampled-4.mha');

opt = 10;
method = 'cgalsurf';

opt = .0045;
method = 'simplify';

tic
[node,elem]=v2s(single(aux.data), 1, opt, method);
toc
hold off
plotmesh(node,elem(:,1:3))

% rename the mesh variables so that they are compatible with the rest of
% our code
tri = elem(:, 1:3);
x = node(:, [2 1 3]);
x(:, 1) = x(:, 1) * aux.axis(2).spacing;
x(:, 2) = x(:, 2) * aux.axis(1).spacing;
x(:, 3) = x(:, 3) * aux.axis(3).spacing;

% some test points
xi = [
    40 50 80
    20 20 20
    80 80 160
    60 40 140
    ];
xi(:, 1) = xi(:, 1) * aux.axis(2).spacing;
xi(:, 2) = xi(:, 2) * aux.axis(1).spacing;
xi(:, 3) = xi(:, 3) * aux.axis(3).spacing;

% plot points
hold on
plot3(xi(:, 1), xi(:, 2), xi(:, 3), 'r*')

% find closest facet and distance to a point
tic
[f, d, p] = cgal_closest_trifacet(tri, x, xi);
toc
