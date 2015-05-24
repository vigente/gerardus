function varargout = elastix_bspline_grid(t)
% ELASTIX_BSPLINE_GRID  Control point grid of B-spline transformation.
%
% [GX, GY] = ELASTIX_BSPLINE_GRID(T)
%
%   T is a struct with the 2D B-spline transform parameters, or a string
%   with the path to a file with the transforms.
%
%   GX, GY are matrices with the coordinates of the B-spline coefficients.
%   This is a similar format to the output of NDGRID. They follow the
%   Matlab convention that X->columns and Y->rows. Note that Elastix
%   parameters follow the convention that X->rows, Y->columns, so they are
%   internally transposed to produced GX, GY.
%
% [GX, GY, ..., GN] = ELASTIX_BSPLINE_GRID(T)
% 
%   If T is an N-dimensional transform, the function produces the
%   corresponding N coordinates.
%
% You can visualize the control point grid using
%
%   plot(gx(:), gy(:), 'o')
%
% If you want to see the position of the image with respect to the control
% points:
%
%   hold on
%   box = [
%       t.Origin
%       t.Origin(1)+(t.Size(1)-1)*t.Spacing(1), t.Origin(2)
%       t.Origin+(t.Size-1).*t.Spacing
%       t.Origin(1), t.Origin(2)+(t.Size(2)-1)*t.Spacing(2)
%       t.Origin
%       ];
%   plot(box(:, 1), box(:, 2), 'r', 'LineWidth', 2)
%
% See also: elastix_bspline_grid2param.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2015 University of Oxford
% Version: 0.1.4
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

% check input arguments
narginchk(1, 1);

% read parameters if provided as a file
if (ischar(t))
    t = elastix_read_file2param(t);
end

% check output arguments
D = t.MovingImageDimension;
nargoutchk(0, D);

% coordinates of grid points, separated by dimension
x = cell(D, 1);
for I = 1:D
    % axis ticks along each dimension
    x{I} = (0:t.GridSize(I)-1) * t.GridSpacing(I) + t.GridOrigin(I);
end

% generate the grid (switch X<->Y cells so that they suit the Matlab
% convention of X=columns and Y=rows)
[varargout{[2 1 3:D]}] = ndgrid(x{[2 1 3:end]});

% split control point displacements into one coordinate per cell
dx = mat2cell(t.TransformParameters, 1, ...
    t.NumberOfParameters / D * ones(1, D));

for I = 1:D
    
    % reshape the vector of coefficients to the same size as the grid
    dx{I} = reshape(dx{I}, t.GridSize);
    
    % transpose cols<->rows to adapt to Matlab's convention
    dx{I} = permute(dx{I}, [2 1 3:D]);
    
    % add displacements to control points
    varargout{I} = varargout{I} + dx{I};
end
