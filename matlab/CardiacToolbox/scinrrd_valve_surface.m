function nrrd = scinrrd_valve_surface(nrrd, x, PARAM, INTERP, KLIM)
% SCINRRD_VALVE_SURFACE  Inter atrio-ventricular surface
%
% NRRD = SCINRRD_VALVE_SURFACE(NRRD0, X)
%
%   NRRD0 is the SCI NRRD struct that contains the cardiac Magnetic
%   Resonance Image (MRI).
%
%   X is a 3-row matrix. Each column has the coordinates of a point from
%   one of the valve annula.
%
%   NRRD is a SCI NRRD struct with a segmentation of the surface that
%   interpolates the valve annula points.
%
% NRRD = SCINRRD_VALVE_SURFACE(NRRD0, X, PARAM, INTERP, KLIM)
%
%   PARAM is a string with the method used to parametrize the surface and
%   X:
%
%     'xy' (default): No change, the X coordinates are kept the same.
%
%     'pca': X points are rotated according to their eigenvectors to make
%     the valve surface as horizontal as possible before interpolating.
%
%     'isomap': Use the Isomap method by [1] to "unfold" the curved surface
%     defined by the valves before interpolating.
%
%   INTERP is a string with the interpolation method:
%
%      'tps' (default): Thin-plate spline.
%
%   KLIM is a scalar factor for the extension of the interpolation domain.
%   By default, KLIM=1 and the interpolation domain is a rectangle that
%   tightly contains X. Sections of the interpolated surface that protude
%   from the image volume are removed.
%
%
% [1] J.B. Tenenbaum, V. de Silva and J.C. Langford, "A Global Geometric
% Framework for Nonlinear Dimensionality Reduction", Science 290(5500):
% 2319-2323, 2000.
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

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010-2011 University of Oxford
% Version: 0.1
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
error(nargchk(2, 5, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% defaults
if (nargin < 3 || isempty(PARAM))
    PARAM='xy';
end
if (nargin < 4 || isempty(INTERP))
    INTERP='tps';
end
if (nargin < 5 || isempty(KLIM))
    KLIM=1;
end

% extract data volume and convert to boolean to save space
nrrd.data = boolean(nrrd.data);

%% map the 3D points (x,y,z) to a 2D domain (u,v)

% this is analogous to computing the knot vector for a curve interpolation.
% The idea is that (x,y)->z is not necessarily a function, due to the valve
% folding over. However, we hope that (x,y,z)->(u,v) is a function
switch PARAM
    
    case 'xy'
        
        % (u,v) is simply (x,y)
        em = x(1:2, :);
        
    case 'pca'
        
        % rotate valve points to make the valve surface as horizontal as
        % possible
        m = mean(x, 2);
        em = x - m(:, ones(1, size(x, 2)));
        eigv = pts_pca(em);
        em = eigv' * em;
        em = em(1:2, :);
        
    case 'isomap'
    
    otherwise
        error('Parametrization method not implemented')
end

%% compute interpolation domain

% get voxel size
res = [nrrd.axis.spacing];

% find box that contains embedded coordinates
emmin = min(em, [], 2);
emmax = max(em, [], 2);

% box size and centroid
delta = emmax - emmin;
boxm = mean([emmax emmin], 2);

% extend the box
emmin = boxm - delta/2*KLIM;
emmax = boxm + delta/2*KLIM;

% generate grid for the embedding box
[gy, gx] = ndgrid(emmin(2):res(1):emmax(2), emmin(1):res(2):emmax(1));


%% compute interpolating surface

% source and target points that will define the warp
%s = em; % don't duplicate data in memory
%t = x; % don't duplicate data in memory

% interpolate
switch INTERP
    case 'tps'
        y = pts_tps_map(em', x', [gx(:) gy(:)])';
    otherwise
        error('Interpolation method not implemented')
end

%% map interpolated surface points to voxels

% convert real world coordinates to indices
idx = round(scinrrd_world2index(y', nrrd.axis));

% remove points outside the volume
badidx = isnan(sum(idx, 2));
idx = idx(~badidx, :);

% construct new segmentation mask. Note that zz has type double.
nrrd.data = nrrd.data * 0;
nrrd.data(sub2ind(size(nrrd.data), ...
    idx(:, 1), idx(:, 2), floor(idx(:, 3)))) = 1;
nrrd.data(sub2ind(size(nrrd.data), ...
    idx(:, 1), idx(:, 2), ceil(idx(:, 3)))) = 1;


end
