% test_elastix_bspline_grid2param.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2015 University of Oxford
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

% create a B-spline transform in elastix format
t.Transform = 'BSplineTransform';
t.NumberOfParameters = 84;
t.TransformParameters = zeros(1, 84);
t.InitialTransformParametersFileName = 'NoInitialTransform';
t.HowToCombineTransforms = 'Compose';
t.FixedImageDimension = 2;
t.MovingImageDimension = 2;
t.FixedInternalImagePixelType = 'float';
t.MovingInternalImagePixelType = 'float';
t.Size = [1331 1101];
t.Index = [0 0];
t.Spacing = [1 1];
t.Origin = [0 0];
t.Direction = [1 0 0 1];
t.UseDirectionCosines = 'true';
t.GridSize = [7 6];
t.GridIndex = [0 0];
t.GridSpacing = [400 400];
t.GridOrigin = [-535 -450];
t.GridDirection = [1 0 0 1];
t.BSplineTransformSplineOrder = 3;
t.UseCyclicTransform = 'false';
t.ResampleInterpolator = 'FinalBSplineInterpolator';
t.FinalBSplineInterpolationOrder = 3;
t.Resampler = 'DefaultResampler';
t.DefaultPixelValue = 0;
t.ResultImageFormat = 'png';
t.ResultImagePixelType = 'unsigned char';
t.CompressResultImage = 'true';

% coordinates of control points
[gx0, gy0] = elastix_bspline_grid(t);

% plot grid
hold off
plot(gx0(:), gy0(:), 'o')
box = [
    t.Origin
    t.Origin(1)+(t.Size(1)-1)*t.Spacing(1), t.Origin(2)
    t.Origin+(t.Size-1).*t.Spacing
    t.Origin(1), t.Origin(2)+(t.Size(2)-1)*t.Spacing(2)
    t.Origin
    ];
hold on
plot(box(:, 1), box(:, 2), 'r', 'LineWidth', 2)

% randomly assign displacement values to the control points
rng(5)
t.TransformParameters = 200 * (rand(1, t.NumberOfParameters)-0.5);

% coordinates of displaced control points
[gx, gy] = elastix_bspline_grid(t);

% plot displaced points
plot(gx(:), gy(:), 'x')

% convert the displacements back from grid to vector format
v = elastix_bspline_grid2param(gx - gx0, gy - gy0);

% check that we have recovered the original transform parameters vector
any(abs(v - t.TransformParameters) > 1e-12)
