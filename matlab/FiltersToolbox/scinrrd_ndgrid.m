function [xg, yg, zg] = scinrrd_ndgrid(nrrd)
% SCINRRD_NDGRID  Generation of arrays for 3D SCI NRRD image volumes
%
% [XG, YG, ZG] = SCINRRD_NDGRID(NRRD)
%
%   NRRD is a SCI image volume.
%
%   XG, YG, ZG are arrays with the x-, y-, z-coordinates of the voxels in
%   NRRD.
%
%   Note that XG values change with columns, and YG values change with
%   rows, to accommodate the usual coordinate frame convention.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010-2012 University of Oxford
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
error(nargchk(1, 1, nargin, 'struct' ));
error(nargoutchk(0, 3, nargout, 'struct'));

% squeeze the non-used first dimension of data
nrrd = scinrrd_squeeze(nrrd);

% get image vertices
vmin = [nrrd.axis.min];
vmax = [nrrd.axis.max];

% get voxel size
res = [nrrd.axis.spacing];

% point to voxel centres instead of voxel left edge
vmin = vmin + res/2;
vmax = vmax + res/2;

% generate 3D grid of coordinates: note the inversion of coordinates,
% necessary so that xg will change with columns, and yg with rows
[yg, xg, zg] = ndgrid(vmin(1):res(1):vmax(1), ...
    vmin(2):res(2):vmax(2), vmin(3):res(3):vmax(3));
