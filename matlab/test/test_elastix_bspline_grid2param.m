% test_elastix_bspline_grid2param.m

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
