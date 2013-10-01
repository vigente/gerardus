% test_cgal_surfsubdivision.m

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
%% Basic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% empty mesh
[tri, x] = cgal_surfsubdivision([], [], 'method', 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tetrahedron subdivision, all approximating schemes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a surface mesh that is a tetrahedron
x = [
    0 0 0
    1 0 0
    0 1 0
    0 0 1
    ];

% apply alphashape
[vol, as] = alphavol(x, 1.3);
tri = as.bnd;
clear as vol

% plot mesh
subplot(1, 2, 1)
hold off
trimesh(tri, x(:,1), x(:,2), x(:,3))
axis([0 1 0 1 0 1])

% surface subdivision
[tri2, x2] = cgal_surfsubdivision(tri, x, 'CatmullClark', 1);

% plot subdivided surface
subplot(1, 2, 2)
hold off
trimesh(tri2, x2(:,1), x2(:,2), x2(:,3))
axis([0 1 0 1 0 1])

% surface subdivision
[tri2, x2] = cgal_surfsubdivision(tri, x, 'Loop', 1);

% plot subdivided surface
subplot(1, 2, 2)
hold off
trimesh(tri2, x2(:,1), x2(:,2), x2(:,3))
axis([0 1 0 1 0 1])

% surface subdivision
[tri2, x2] = cgal_surfsubdivision(tri, x, 'DooSabin', 1);

% plot subdivided surface
subplot(1, 2, 2)
hold off
trimesh(tri2, x2(:,1), x2(:,2), x2(:,3))
axis([0 1 0 1 0 1])

% surface subdivision
[tri2, x2] = cgal_surfsubdivision(tri, x, 'Sqrt3', 1);

% plot subdivided surface
subplot(1, 2, 2)
hold off
trimesh(tri2, x2(:,1), x2(:,2), x2(:,3))
axis([0 1 0 1 0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cube subdivision, all approximating schemes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a surface mesh that is a cube
x = [
    0 0 0
    0 0 1
    0 1 0
    0 1 1
    1 0 0
    1 0 1
    1 1 0
    1 1 1
    ];

% apply alphashape
[vol, as] = alphavol(x, 1.3);
tri = as.bnd;
clear as vol

% plot mesh
subplot(1, 2, 1)
hold off
trimesh(tri, x(:,1), x(:,2), x(:,3))
axis([0 1 0 1 0 1])

% surface subdivision
[tri2, x2] = cgal_surfsubdivision(tri, x, 'CatmullClark', 1);

% plot subdivided surface
subplot(1, 2, 2)
hold off
trimesh(tri2, x2(:,1), x2(:,2), x2(:,3))
axis([0 1 0 1 0 1])

% surface subdivision
[tri2, x2] = cgal_surfsubdivision(tri, x, 'Loop', 1);

% plot subdivided surface
subplot(1, 2, 2)
hold off
trimesh(tri2, x2(:,1), x2(:,2), x2(:,3))
axis([0 1 0 1 0 1])

% surface subdivision
[tri2, x2] = cgal_surfsubdivision(tri, x, 'DooSabin', 1);

% plot subdivided surface
subplot(1, 2, 2)
hold off
trimesh(tri2, x2(:,1), x2(:,2), x2(:,3))
axis([0 1 0 1 0 1])

% surface subdivision
[tri2, x2] = cgal_surfsubdivision(tri, x, 'Sqrt3', 1);

% plot subdivided surface
subplot(1, 2, 2)
hold off
trimesh(tri2, x2(:,1), x2(:,2), x2(:,3))
axis([0 1 0 1 0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cube subdivision, one scheme, different number of iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a surface mesh that is a cube
x = [
    0 0 0
    0 0 1
    0 1 0
    0 1 1
    1 0 0
    1 0 1
    1 1 0
    1 1 1
    ];

% apply alphashape
[vol, as] = alphavol(x, 1.3);
tri = as.bnd;
clear as vol

% plot mesh
subplot(2, 2, 1)
hold off
trimesh(tri, x(:,1), x(:,2), x(:,3))
axis([0 1 0 1 0 1])

levels = 1:2:7;

for I = 1:length(levels)
    
    % surface subdivision
    [tri2, x2] = cgal_surfsubdivision(tri, x, 'CatmullClark', levels(I));
    
    % plot subdivided surface
    subplot(2, 2, I)
    hold off
    trimesh(tri2, x2(:,1), x2(:,2), x2(:,3))
    axis([0 1 0 1 0 1])
    title(['Iter = ' num2str(levels(I))])
    
end
