function nrrd = scinrrd_estimate_bias_field(nrrd, x, a, ANTIA)
% SCINRRD_ESTIMATE_BIAS_FIELD  Estimate MRI bias field
%
%   This function provides an estimate of the bias field from a magnetic
%   resonance image (MRI).
%
%   The function samples an SCI NRRD image in locations provided by the
%   user (typically, corresponding to background voxels), and then creates
%   another image that interpolates the sampled intensity values.
%
%   The locations can be selected e.g. using our Spline Tool extension to
%   the Seg3D platform (click the points, and then export the control
%   points).
%
%   Note that because the interpolant is an interpolating thin-plate
%   spline (TPS), this function is not robust in the presence of noise
%   (this can be solved e.g. by low-pass filtering the image). Ideally, in
%   the future this function should use an approximating TPS instead of an
%   interpolating one.
%
% NRRD2 = SCINRRD_ESTIMATE_BIAS_FIELD(NRRD, X)
%
%   NRRD is an image provided in SCI NRRD format.
%
%   X is a 3-colum matrix with the real world coordinates of the sampling
%   points. Note that each point is rounded to the closest voxel centre.
%
%   NRRD2 is the estimated bias field in SCI NRRD format.
%
% NRRD2 = SCINRRD_ESTIMATE_BIAS_FIELD(NRRD, X, A, ANTIA)
%
%   A is a scaling factor. Because using the TPS to interpolate all voxels
%   can be rather slow, and the bias field is anyway a slow varying field,
%   it's convenient to first quickly reduce the image size by A (using
%   bilinear interpolation), interpolate the bias field in the smaller
%   image with the TPS, and then expand to the original size. By default, A
%   = 1.0 and no rescaling is used.
%
%   ANTIA is a boolean flag to run an anti-aliasing filter before
%   downsampling the volume. By default, ANTIA = true.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
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
error(nargchk(2, 4, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% defaults
if (nargin < 3 || isempty(a))
    a = 1.0;
end
if (nargin < 4 || isempty(ANTIA))
    ANTIA = true;
end
    
% squeeze volume
nrrd = scinrrd_squeeze(nrrd);

% save input parameters for the end
size_in = num2cell([nrrd.axis.size]);
spacing_in = num2cell([nrrd.axis.spacing]);

% recompute the scaling factor so that we obtain an integer number of
% voxels in the output volume
a = a([1 1 1]);
size_out = num2cell(round([nrrd.axis.size] .* a));
a = [size_out{:}] ./ [nrrd.axis.size];
[nrrd.axis.size] = deal(size_out{:});

% recompute the voxel spacing according to the adjusted scaling factor
spacing_out = num2cell([nrrd.axis.spacing] ./ a);
[nrrd.axis.spacing] = deal(spacing_out{:});

if ANTIA
    % compute low-pass anti-aliasing filter
    sigma = 1 ./ a;
    h = fspecial3('gaussian', sigma * 4, sigma);
    
    % low pass filtering of the image
    nrrd.data = imfilter(nrrd.data, h);
end

% % compute size of image in the frequency domain (*4 to avoid border effects
% % due to cylic convolution)
% sz = 2.^ceil(log2([size_in{:}])) * 4;
% 
% % low pass filtering of the image (frequency domain)
% nrrd.data = real(ifftn(fftn(nrrd.data, sz) .* fftn(h, sz)));
% nrrd.data = nrrd.data(1:size_in{1}, 1:size_in{2}, 1:size_in{3});

% reduce volume size
nrrd.data = tformarray(nrrd.data, ...
    maketform('affine', [a(1) 0 0 0; 0 a(2) 0 0; 0 0 a(3) 0; 0 0 0 1]), ...
    makeresampler('linear', 'replicate'), ...
    1:length(nrrd.axis), 1:length(nrrd.axis), ...
    [nrrd.axis.size], [], []);

% extract size of the volume
n = [nrrd.axis.size];

% compute coordinates of the two extreme voxels that define the image
% volume
cmin = scinrrd_index2world([1 1 1], nrrd.axis);
cmax = scinrrd_index2world(n, nrrd.axis);

% if we have point coordinates with values like 100, 400, matrix L for the
% thin-plate spline weight computation is badly scaled. Thus, we make use
% of the invariability of thin-plate splines to scaling and make all the
% coordinate values <= 5.0 (if we make them <= 1.0, interpn gives an error
% saying that they grid is not monotonic)
K = max(cmax)/5;

% Seg3D provides continuous point coordinates over the whole image volume.
% But in fact, intensity values correspond to the voxel centre, so we have
% to convert the coordinates to indices and round them
x = scinrrd_world2index(x, nrrd.axis);
x = round(x);

% % get intensity values corresponding to the sampling points
idx = sub2ind(n, x(:, 1), x(:, 2), x(:, 3));
v = nrrd.data(idx);

% convert back to coordinates
x = scinrrd_index2world(x, nrrd.axis);

% the warp will be defined as from the xyz-points to the corresponding
% intensity values

% compute weights for thin-plate spline interpolation
w = pts_tps_weights( x/K, v );

% interpolate intensity values for each point in the grid
[gx, gy, gz] = meshgrid(...
    linspace(cmin(1)/K, cmax(1)/K, n(1)), ...
    linspace(cmin(2)/K, cmax(2)/K, n(2)), ...
    linspace(cmin(3)/K, cmax(3)/K, n(3)));
y = pts_tps_map( x/K, v, [ gx(:) gy(:) gz(:) ], w, true, false );

% create output volume
nrrd.data = reshape(single(y), n);

% recover original size
nrrd.data = tformarray(nrrd.data, ...
    maketform('affine', [1/a(1) 0 0 0; 0 1/a(2) 0 0; 0 0 1/a(3) 0; 0 0 0 1]), ...
    makeresampler('linear', 'replicate'), ...
    1:length(nrrd.axis), 1:length(nrrd.axis), ...
    [size_in{:}], [], []);

[nrrd.axis.size] = deal(size_in{:});
[nrrd.axis.spacing] = deal(spacing_in{:});
