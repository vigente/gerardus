% test_cgal_meshseg.m

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

% empty segmentation 
[tri, x] = cgal_meshseg([], 1, 30, .25, .25, [3, 2.5, 1.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh small rectangular segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create segmentation
im = zeros(7,6,4, 'double');
im(2:6,2:5,2:3) = 1;

% defaults
[tri, x] = cgal_meshseg(im, 0.5);

% plot segmentation
hold off
trimesh(tri, x(:,1), x(:,2), x(:,3))
axis equal

% mesh the segmentation with other input parameters
[tri, x] = cgal_meshseg(im, 1, 30, .25, .25);

% plot segmentation
hold off
trimesh(tri, x(:,1), x(:,2), x(:,3))
axis equal

% reduce the iso value towards 0, so that the isosurface will grow. This
% test is to check that the bounding sphere computed internally will be
% large enough to contain the isosurface for any kind of iso value between
% 0+epsilon and 1
isoval = linspace(0.01, 1, 4);
for I = 1:4
    
    % compute segmentation
    [tri, x] = cgal_meshseg(im, isoval(I), 30, .25, .25);
    
    % plot result
    subplot(2, 2, I)
    hold off
    trimesh(tri, x(:,1), x(:,2), x(:,3))
    axis equal
    title(['isoval = ' num2str(isoval(I))])
    
end

% plot with a 2-D view
for I = 1:4
    
    % compute segmentation
    [tri, x] = cgal_meshseg(im, isoval(I), 30, .25, .25);
    
    % plot result
    subplot(2, 2, I)
    hold off
    trimesh(tri, x(:,1), x(:,2), x(:,3))
    view(2)
    axis equal
    title(['isoval = ' num2str(isoval(I))])
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Try input segmentation types other than 'double'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create segmentation
im = zeros(7,6,4, 'double');
im(2:6,2:5,2:3) = 1;

% boolean and int64 segmentation: we throw an error because the CGAL mesher
% enters an infinite loop, for some reason
[tri, x] = cgal_meshseg(boolean(im), 0.5, 30, .25, .25);
[tri, x] = cgal_meshseg(int64(im), 0.5, 30, .25, .25);

% mesh the segmentation
[tri, x] = cgal_meshseg(uint8(im), 0.5, 30, .25, .25);

% plot segmentation
subplot(2, 2, 1)
hold off
trimesh(tri, x(:,1), x(:,2), x(:,3))
axis equal
title('type = uint8')

% mesh the segmentation
[tri, x] = cgal_meshseg(uint16(im), 0.5, 30, .25, .25);

% plot segmentation
subplot(2, 2, 2)
hold off
trimesh(tri, x(:,1), x(:,2), x(:,3))
axis equal
title('type = uint16')

% mesh the segmentation
[tri, x] = cgal_meshseg(int32(im), 0.5, 30, .25, .25);

% plot segmentation
subplot(2, 2, 3)
hold off
trimesh(tri, x(:,1), x(:,2), x(:,3))
axis equal
title('type = int32')

% mesh the segmentation
[tri, x] = cgal_meshseg(single(im), 0.5, 30, .25, .25);

% plot segmentation
subplot(2, 2, 4)
hold off
trimesh(tri, x(:,1), x(:,2), x(:,3))
axis equal
title('type = single')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Effect of bounding sphere centre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mesh the segmentation: default centre at centre of mass
[tri, x] = cgal_meshseg(im, 0.5, 30, .25, .25, []);

% plot segmentation
subplot(1, 2, 1)
hold off
trimesh(tri, x(:,1), x(:,2), x(:,3))
axis equal
title('Centre of mass')

% mesh the segmentation: sphere centre close to the surface
[tri, x] = cgal_meshseg(im, 0.5, 30, .25, .25, [2.5 3 2.5]);

% plot segmentation
subplot(1, 2, 2)
hold off
trimesh(tri, x(:,1), x(:,2), x(:,3))
axis equal
title('Close to the top surface')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Manifold vs non-manifold meshes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mesh the segmentation: Manifold
[tri, x] = cgal_meshseg(im, 0.5, 30, .25, .25, [], true);

% plot segmentation
subplot(1, 2, 1)
hold off
trimesh(tri, x(:,1), x(:,2), x(:,3))
axis equal
title('Manifold mesh')

% mesh the segmentation: Non-manifold
[tri, x] = cgal_meshseg(im, 0.5, 30, .25, .25, [], false);

% plot segmentation
subplot(1, 2, 2)
hold off
trimesh(tri, x(:,1), x(:,2), x(:,3))
axis equal
title('Non-manifold mesh')

