function [xg, yg, zg] = scimat_ndgrid(scimat, ri, ci, si)
% SCIMAT_NDGRID  Generation of arrays for 3D SCIMAT image volumes.
%
% [XG, YG, ZG] = scimat_ndgrid(SCIMAT)
%
%   SCIMAT is a struct with an image volume (see "help scimat" for
%   details).
%
%   XG, YG, ZG are arrays with the x-, y-, z-coordinates of the voxels in
%   SCIMAT.
%
%   Note that XG values change with columns, and YG values change with
%   rows, to accommodate the usual coordinate frame convention.
%
% [XG, YG, ZG] = scimat_ndgrid(SCIMAT, RI, CI, SI)
%
%   RI, CI, SI are vectors of voxel indices (rows, columns and slices). The
%   output grid will be generated only for the corresponding image block.
%   For example, RI=1:4, CI=3:6, SI=5:7 means that the grid of coordinates
%   will be generated only for the block of voxels (1:4, 3:6, 5:7).
%
% See also: ndgrid.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010,2014 University of Oxford
% Version: 0.3.0
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
narginchk(1, 4);
nargoutchk(0, 3);

% squeeze the non-used first dimension of data
scimat = scimat_squeeze(scimat);

% defaults
if (nargin < 2 || isempty(ri))
    ri = 1:size(scimat.data, 1);
end
if (nargin < 3 || isempty(ci))
    ci = 1:size(scimat.data, 2);
end
if (nargin < 4 || isempty(si))
    si = 1:size(scimat.data, 3);
end

% local variables
res = [scimat.axis.spacing];

% convert indices to real world coordinates
r = (ri - 1) * res(1) + scimat.axis(1).min + res(1)/2;
c = (ci - 1) * res(2) + scimat.axis(2).min + res(2)/2;
s = (si - 1) * res(3) + scimat.axis(3).min + res(3)/2;

% generate 3D grid of coordinates: note the inversion of coordinates,
% necessary so that xg will change with columns, and yg with rows
[yg, xg, zg] = ndgrid(r, c, s);
