% test_itk_icp_registration.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
% Version: 0.0.1
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
%% empty input

% register y onto x
itk_icp_registration([], [])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% register points from a rectangular parallelepiped to points from a 
%% cylinder

x = [
    cosd(linspace(0, 360, 8)') sind(linspace(0, 360, 8)') zeros(8, 1)
    cosd(linspace(0, 360, 8)') sind(linspace(0, 360, 8)') zeros(8, 1)-.5
    cosd(linspace(0, 360, 8)') sind(linspace(0, 360, 8)') zeros(8, 1)+.5
    ];

y = [
    0 0 0
    0 1 0
    1 0 0
    1 1 0
    0 0 -.5
    0 1 -.5
    1 0 -.5
    1 1 -.5
    0 0 .5
    0 1 .5
    1 0 .5
    1 1 .5
    ];
y(:, 1) = y(:, 1) + 1;
y(:, 2) = y(:, 2) + 1;

% plot points
hold off
plot3(x(:, 1), x(:, 2), x(:, 3), '.')
hold on
plot3(y(:, 1), y(:, 2), y(:, 3), 'r*')
axis equal

% register y onto x
[y2, t] = itk_icp_registration(x, y);

% plot points
hold off
plot3(x(:, 1), x(:, 2), x(:, 3), 'o')
hold on
plot3(y2(:, 1), y2(:, 2), y(:, 3), 'r*')
axis equal

