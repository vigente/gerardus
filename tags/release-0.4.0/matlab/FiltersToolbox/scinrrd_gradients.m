function [dx, dy, dxx, dyy, dxy] = scinrrd_gradients(nrrd, method)
% SCINRRD_GRADIENTS  Compute 1st and 2nd order image gradients
%
% [DX, DY, DXX, DYY, DXY] = SCINRRD_GRADIENTS(NRRD)
%
%   This function computes the first and second order gradients of an image
%   provided in SCI NRRD format. The voxel size in each dimension is taken
%   into account when computing the gradients.
%
%   NRRD is the struct with the data.
%
%   DX, DY are the first order gradients of the image I in the X and Y
%   directions, respectively, dI/dx, dI/dy (Note: the X direction
%   corresponds to *columns*, and the Y direction corresponds to *rows*).
%
%   DXX, DYY, DXY are the second order gradients, d^2(I)/dx^2, d^2(I)/dy^2,
%   d^2(I)/dxdy. Second order gradients are computed with an explicit
%   formula using a 3x3 pixel neighbourhood.
%
% ... = SCINRRD_GRADIENTS(NRRD, METHOD)
%
%   The second order gradients are computed using a explicit formula by
%   default (METHOD='default'). However, it is also possible to compute it
%   running function gradient() twice (METHOD='approx'). Note that 'approx'
%   uses a 5x5 neighbourhood to compute the gradient, instead of 3x3, so it
%   is less accurate. It's also slower.
%
%   Note on SCI NRRD: Software applications developed at the University of
%   Utah Scientific Computing and Imaging (SCI) Institute, e.g. Seg3D,
%   internally use NRRD volumes to store medical data.
%
%   When label volumes (segmentation masks) are saved to a Matlab file
%   (.mat), they use a struct called "scirunnrrd" to store all the NRRD
%   information:
%
%   >>  scirunnrrd
%
%   scirunnrrd = 
%
%          data: [4-D uint8]
%          axis: [4x1 struct]
%      property: []

% Author: Ramon Casero <rcasero@gmail.com>, Vicente Grau
% Copyright © 2010 University of Oxford
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

% check arguments
error( nargchk( 1, 2, nargin, 'struct' ) );
error( nargoutchk( 0, 5, nargout, 'struct' ) );

% defaults
if (nargin < 2)
    method = 'default';
end

% compute first order gradients (beware: spacing(1) is for rows, the
% Y-coordinate)
[dx, dy] = gradient(nrrd.data, nrrd.axis(2).spacing, ...
    nrrd.axis(1).spacing, nrrd.axis(3).spacing);

% compute second order gradients (we assume dxy = dyx)
dxx = zeros(size(dx));
dyy = zeros(size(dx));
dxy = zeros(size(dx));

if (strcmp(method, 'default'))
    dxx(:, 2:end-1, :) = (-2 * nrrd.data(:, 2:end-1, :) ...
        + nrrd.data(:, 1:end-2, :) + nrrd.data(:, 3:end, :)) ...
        / nrrd.axis(2).spacing^2;
    dyy(2:end-1, :, :) = (-2 * nrrd.data(2:end-1, :, :) ...
        + nrrd.data(1:end-2, :, :) + nrrd.data(3:end, :, :)) ...
        / nrrd.axis(1).spacing^2;
    dxy(2:end-1, 2:end-1, :) = (...
        -nrrd.data(1:end-2, 3:end, :) - ...
        nrrd.data(3:end, 1:end-2, :) + ...
        nrrd.data(1:end-2, 1:end-2, :) + ...
        nrrd.data(3:end, 3:end, :) ...
        ) ...
        / nrrd.axis(1).spacing / nrrd.axis(2).spacing;

elseif (strcmp(method, 'approx'))
    % instead of using the explicit formula for the second order gradients, we
    % can use function gradient() twice
    warning('Using inaccurate approximation')
    dxx = gradient(dx, nrrd.axis(2).spacing, ...
        nrrd.axis(1).spacing, nrrd.axis(3).spacing);
    [dxy, dyy] = gradient(dy, nrrd.axis(2).spacing, ...
        nrrd.axis(1).spacing, nrrd.axis(3).spacing);
end
