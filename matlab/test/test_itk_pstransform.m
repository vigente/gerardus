% test_itk_pstransform.m
%
% Script to test MEX function itk_pstransform().

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
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
%% Test 1:
%% 2D grid with 4 landmarks

% source landmarks
x = [...
    .25 .25;...
    .25 .75;...
    .75 .25;...
    .75 .75...
    ];

% target landmarks
y = [...
    .25+.1 .25+.1;...
    .25-.2 .75+.1;...
    .75-.1 .25-.15;...
    .75+.15 .75+.2...
    ];

% grid points to be transformed
[gu, gv] = meshgrid(linspace(0, 1, 11), linspace(0, 1, 11));
xi = [gu(:) gv(:)];

% loop every available transform
subplot(1, 1, 1)
hold off
for TRANSF = {'elastic', 'elasticr', 'tps', 'tpsr2', 'volume', 'bspline'}
    
    % apply transform to grid
    yi = itk_pstransform(TRANSF{:}, x, y, xi);
    
    % plot points
    hold off
    plot(x(:, 1), x(:, 2), 'or')
    hold on
    plot(y(:, 1), y(:, 2), 'xg')
    for J = 1:size(x, 1)
        plot([x(J, 1) y(J, 1)], [x(J, 2) y(J, 2)], 'r')
    end
    plot(xi(:, 1), xi(:, 2), '.')
    plot(yi(:, 1), yi(:, 2), '.k')
    title(TRANSF{:})
    
    pause
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test 2:
%% 3D segmentation of bent vessel

% load segmentation points 'xi', skeleton points 'x', and skeleton
% parameterization 't'
load('test/ps-bent-vessel-3d-3.mat')

% choose by hand a few points in the segmentation, and use them as source
% landmarks
x = [...
    0.0094    0.0153    0.0142; ...
    0.0092    0.0139    0.0127; ...
    0.0094    0.0137    0.0108; ...
    0.0102    0.0127    0.0088; ...
    0.0115    0.0118    0.0085];

% move the source landmarks a bit to create the target landmarks
y = [...
    0.0094-.001    0.0153+.001    0.0142+.001; ...
    0.0092+.0005   0.0139-.001    0.0127-.0005; ...
    0.0094-.0004    0.0137-.0002    0.0108; ...
    0.0102+.0004    0.0127+.0002    0.0088; ...
    0.0115    0.0118    0.0085-.001];

% loop every available transform
for TRANSF = {'elastic', 'elasticr', 'tps', 'tpsr2', 'volume', 'bspline'}
    
    % apply transform
    yi = itk_pstransform(TRANSF{:}, x, y, xi);
    
    % plot points
    subplot(1, 2, 1)
    hold off
    plot3(xi(:, 1), xi(:, 2), xi(:, 3), '.')
    hold on
    plot3(x(:, 1), x(:, 2), x(:, 3), '*r')
    plot3(y(:, 1), y(:, 2), y(:, 3), '*g')
    for J = 1:size(x, 1)
        plot3([x(J, 1) y(J, 1)], [x(J, 2) y(J, 2)], [x(J, 3) y(J, 3)], 'r')
    end
    title(TRANSF{:})
    axis xy equal
    view(-40, 20)
    
    subplot(1, 2, 2)
    hold off
    plot3(yi(:, 1), yi(:, 2), yi(:, 3), '.')
    hold on
    plot3(x(:, 1), x(:, 2), x(:, 3), '*r')
    plot3(y(:, 1), y(:, 2), y(:, 3), '*g')
    for J = 1:size(x, 1)
        plot3([x(J, 1) y(J, 1)], [x(J, 2) y(J, 2)], [x(J, 3) y(J, 3)], 'r')
    end
    title(TRANSF{:})
    axis xy equal
    view(-40, 20)
    
    pause
    
end
