function nrrd = scinrrd_estimate_bias_field(nrrd, x, a)
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
%   Note that the interpolant is an interpolating thin-plate spline (TPS),
%   but instead of just sampling the selected voxels, a neighbourhood is
%   sampled, and used to locally low-pass filter the image with a Gaussian
%   filter. Thus, this function is somehow robust in the presence of noise.
%   Alternatively, this function could use an approximating TPS instead of
%   an interpolating one.
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
% NRRD2 = SCINRRD_ESTIMATE_BIAS_FIELD(NRRD, X, A)
%
%   A is a scaling factor. Because using the TPS to interpolate all voxels
%   can be rather slow, and the bias field is anyway a slow varying field,
%   it's convenient to first quickly reduce the image size by A (using
%   bilinear interpolation), interpolate the bias field in the smaller
%   image with the TPS, and then expand to the original size. By default, A
%   = 1.0 and no rescaling is used.

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
error(nargchk(2, 3, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% defaults
if (nargin < 3 || isempty(a))
    a = 1.0;
end
    
% squeeze volume
nrrd = scinrrd_squeeze(nrrd);

% size of input volume
nin = [nrrd.axis.size];

% recompute the scaling factor so that we obtain an integer number of
% voxels in the output volume
a = a([1 1 1]);
nout = round(nin .* a);
a = nout ./ nin;

% Seg3D provides continuous point coordinates over the whole image volume.
% But in fact, intensity values correspond to the voxel centre, so we have
% to convert the coordinates to indices and round them
x = scinrrd_world2index(x, nrrd.axis);
x = round(x);

% compute low-pass anti-aliasing filter
sigma = 1 ./ a;
halfsz = round(sigma * 4 / 2);
h = fspecial3('gaussian', halfsz*2+1, sigma);

% loop each sampling point. For each sampling point we want to sample a
% whole cube around it, so that we can apply a low-pass filter without
% having to filter the whole image volume
v = zeros(size(x, 1), 1);
for I = 1:size(x, 1)
    
    % get indices of the area around the sampling point to sample. If the
    % sampling point is too close to the edge, we have to be careful to not
    % overflow
    idxr = max(1, x(I, 1)-halfsz(1)):min(nin(1), x(I, 1)+halfsz(1));
    idxc = max(1, x(I, 2)-halfsz(2)):min(nin(2), x(I, 2)+halfsz(2));
    idxs = max(1, x(I, 3)-halfsz(3)):min(nin(3), x(I, 3)+halfsz(3));
    
    % if necessary, crop low-pass filter so that it has the same size as
    % the sampling area
    hbox = h(...
        max(1, idxr - x(I, 1) + halfsz(1) + 1), ...
        max(1, idxc - x(I, 2) + halfsz(2) + 1), ...
        max(1, idxs - x(I, 3) + halfsz(3) + 1));
    
    % extract image area around the sampling point
    im = nrrd.data(idxr, idxc, idxs);
    
    v(I) = sum(hbox(:) .* im(:)) / sum(hbox(:));
%     v(I) = sum(hbox(:) .* im(:)); % this is how Matlab's imfilter() does
%     it, which is incorrect near the edges
    
end

% convert back to coordinates
x = scinrrd_index2world(x, nrrd.axis);

% compute coordinates of the two extreme voxels that define the image
% volume
cmin = scinrrd_index2world([1 1 1], nrrd.axis);
cmax = scinrrd_index2world(nin, nrrd.axis);

% if we have point coordinates with values like 100, 400, matrix L for the
% thin-plate spline weight computation is badly scaled. Thus, we make use
% of the invariability of thin-plate splines to scaling and make all the
% coordinate values <= 5.0 (if we make them <= 1.0, interpn gives an error
% saying that they grid is not monotonic)
K = max(cmax)/5;

% the warp will be defined as from the xyz-points to the corresponding
% intensity values

% compute what would be the size of the image if we downsample it
nmid = round(nin .* a);

% compute weights for thin-plate spline interpolation
w = pts_tps_weights( x/K, v );

% interpolate intensity values for each point in the grid
[gx, gy, gz] = meshgrid(...
    linspace(cmin(1)/K, cmax(1)/K, nmid(1)), ...
    linspace(cmin(2)/K, cmax(2)/K, nmid(2)), ...
    linspace(cmin(3)/K, cmax(3)/K, nmid(3)));
nrrd.data = pts_tps_map( x/K, v, [ gx(:) gy(:) gz(:) ], w, true, false );

% reshape the interpolated values to go from a vector to an image volume
nrrd.data = reshape(single(nrrd.data), nmid);

% recover original size
nrrd.data = tformarray(nrrd.data, ...
    maketform('affine', [1/a(1) 0 0 0; 0 1/a(2) 0 0; 0 0 1/a(3) 0; 0 0 0 1]), ...
    makeresampler('linear', 'replicate'), ...
    1:length(nrrd.axis), 1:length(nrrd.axis), ...
    nin, [], []);
