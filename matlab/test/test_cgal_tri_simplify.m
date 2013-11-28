% test_cgal_tri_simplify.m

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
%% Sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create meshing of sphere
rng(5)
[x, tri] = ParticleSampleSphere('N', 50);

% plot sphere mesh
subplot(2, 2, 1)
hold off
trisurf(tri, x(:, 1), x(:, 2), x(:, 3));
title('Original')
axis equal

% number of edge ratios for mesh simplification
r = [.75 .50 .25];

for I = 1:length(r)

    % simplify the mesh at several ratios
    [tri2, x2] = cgal_tri_simplify(tri, x, r(I));
    
    % plot simplified sphere mesh
    subplot(2, 2, I + 1)
    hold off
    trisurf(tri2, x2(:, 1), x2(:, 2), x2(:, 3));
    axis equal
    title(['Edge ratio = ' num2str(r(I))])

end

% overlap both results
subplot(1, 1, 1)
hold off
trisurf(tri, x(:, 1), x(:, 2), x(:, 3));
hold on
trisurf(tri2, x2(:, 1), x2(:, 2), x2(:, 3));
axis equal
