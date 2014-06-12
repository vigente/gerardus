function [dx, dy, dxx, dyy, dxy] = scimat_gradients(scimat, method)
% SCIMAT_GRADIENTS  Compute 1st and 2nd order image gradients
%
% [DX, DY, DXX, DYY, DXY] = scimat_gradients(SCIMAT)
%
%   This function computes the first and second order gradients of an image
%   provided in SCIMAT format. The voxel size in each dimension is taken
%   into account when computing the gradients.
%
%   SCIMAT is the struct with the image (see "help scimat" for details).
%
%   DX, DY are the first order gradients of the image I in the X and Y
%   directions, respectively, dI/dx, dI/dy (Note: the X direction
%   corresponds to *columns*, and the Y direction corresponds to *rows*).
%
%   DXX, DYY, DXY are the second order gradients, d^2(I)/dx^2, d^2(I)/dy^2,
%   d^2(I)/dxdy. Second order gradients are computed with an explicit
%   formula using a 3x3 pixel neighbourhood.
%
% ... = scimat_gradients(SCIMAT, METHOD)
%
%   The second order gradients are computed using a explicit formula by
%   default (METHOD='default'). However, it is also possible to compute it
%   running function gradient() twice (METHOD='approx'). Note that 'approx'
%   uses a 5x5 neighbourhood to compute the gradient, instead of 3x3, so it
%   is less accurate. It's also slower.

% Author: Ramon Casero <rcasero@gmail.com>, Vicente Grau
% Copyright © 2010,2014 University of Oxford
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
narginchk(1, 2);
nargoutchk(0, 5);

% defaults
if (nargin < 2)
    method = 'default';
end

% compute first order gradients (beware: spacing(1) is for rows, the
% Y-coordinate)
[dx, dy] = gradient(scimat.data, scimat.axis(2).spacing, ...
    scimat.axis(1).spacing, scimat.axis(3).spacing);

% compute second order gradients (we assume dxy = dyx)
dxx = zeros(size(dx));
dyy = zeros(size(dx));
dxy = zeros(size(dx));

if (strcmp(method, 'default'))
    dxx(:, 2:end-1, :) = (-2 * scimat.data(:, 2:end-1, :) ...
        + scimat.data(:, 1:end-2, :) + scimat.data(:, 3:end, :)) ...
        / scimat.axis(2).spacing^2;
    dyy(2:end-1, :, :) = (-2 * scimat.data(2:end-1, :, :) ...
        + scimat.data(1:end-2, :, :) + scimat.data(3:end, :, :)) ...
        / scimat.axis(1).spacing^2;
    dxy(2:end-1, 2:end-1, :) = (...
        -scimat.data(1:end-2, 3:end, :) - ...
        scimat.data(3:end, 1:end-2, :) + ...
        scimat.data(1:end-2, 1:end-2, :) + ...
        scimat.data(3:end, 3:end, :) ...
        ) ...
        / scimat.axis(1).spacing / scimat.axis(2).spacing;

elseif (strcmp(method, 'approx'))
    % instead of using the explicit formula for the second order gradients, we
    % can use function gradient() twice
    warning('Using inaccurate approximation')
    dxx = gradient(dx, scimat.axis(2).spacing, ...
        scimat.axis(1).spacing, scimat.axis(3).spacing);
    [dxy, dyy] = gradient(dy, scimat.axis(2).spacing, ...
        scimat.axis(1).spacing, scimat.axis(3).spacing);
end
