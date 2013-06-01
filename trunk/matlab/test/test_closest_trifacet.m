% test__closest_trifacet.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
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
    ];

% plot points
hold on
plot3(xi(:, 1), xi(:, 2), xi(:, 3), 'r*')

% find closest facet and distance to a point
[f, d] = closest_trifacet(as.bnd, x, xi);

% plot results
for I = 1:size(xi, 1)

    % plot tetrahedron
    hold off
    trisurf(as.bnd, x(:, 1), x(:, 2), x(:, 3), 'FaceColor', 'blue', 'FaceAlpha', .5)

    % plot point
    hold on
    plot3(xi(I, 1), xi(I, 2), xi(I, 3), 'r*')
    
    % plot closest facet
    trisurf(as.bnd(f(I), :), x(:, 1), x(:, 2), x(:, 3), 'FaceColor', 'red', 'FaceAlpha', .5)
    
    pause
    
end
