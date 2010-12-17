function x = scinrrd_index2world(idx, ax)
% SCINRRD_INDEX2WORLD  Convert data volume indices to real world
% coordinates for NRRD volumes created by SCI applications (e.g. Seg3D)
%
%   Function SCINRRD_INDEX2WORLD() maps between the indices of the 
%   image volume used to store the voxel intensity values, and the
%   real world coordinates of points within the NRRD data volume.
%
%      [r, c, s] -> [x, y, z]
%
%   This function assumes that input indices are in 
%   (row, column, slice)-order, corresponding to (y, x, z). However, output
%   world coordinates are given in (x, y, z)-order.
%
%   For points that are not within the data volume, the returned
%   coordinates are "NaN".
%
% X = SCINRRD_INDEX2WORLD(IDX, AXIS)
%
%   X is a 3-column matrix where each row contains the real world
%   (x, y, z)-coordinates of a point.
%
%   IDX has the same size as X, and the voxel indices in 
%   (row, column, slice)-order, that corresponds to (y, x, z)-order.
%
%   AXIS is the 4x1 struct array nrrd.axis from an SCI NRRD struct. It
%   contains the following fields:
%
%   >> nrrd.axis
%
%   4x1 struct array with fields:
%       size
%       spacing
%       min
%       max
%       center
%       label
%       unit
%
% Example:
%
% >> x = scinrrd_index2world([55, 189, 780], nrrd.axis)
%
% x =
%
%     0.0100    0.0110    0.0200
%
%
%   Note on SCI NRRD: Software applications developed at the University of
%   Utah Scientific Computing and Imaging (SCI) Institute, e.g. Seg3D,
%   internally use NRRD volumes to store medical data.
%
%   When data or label volumes are saved to a Matlab file (.mat), they use
%   a struct called "scirunnrrd" to store all the NRRD information:
%
%   >>  scirunnrrd
%
%   scirunnrrd = 
%
%          data: [4-D uint8]
%          axis: [4x1 struct]
%      property: []
%
%
% See also: scinrrd_world2index.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2009-2010 University of Oxford
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
error( nargchk( 2, 2, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

if ( size( idx, 2 ) ~= 3 )
    error( 'IDX must be a 3-column matrix, so that each row has the 3 indices of a voxel' )
end

% init output
x = zeros( size( idx ) );

% extract parameters
xmin = [ ax.min ];
xmax = [ ax.max ];
dx = [ ax.spacing ];
% remove dummy dimension, if present
if ( length( xmin ) == 4 )
    xmin = xmin( 2:end );
    xmax = xmax( 2:end );
    dx = dx( 2:end );
end

% number of dimensions (we expect D=3, but in case this gets more general)
D = length( dx );

% convert indices to real world coordinates
for I = 1:D
    x( :, I ) = ( idx( :, I ) - 1 ) * dx( I ) + xmin( I );
end

% find coordinates that are outside the volume
for I = 1:D
    x( x( :, I ) < xmin( I ) | x( :, I ) > xmax( I ), I ) = NaN;
end

% (y, x, z) => (x, y, z)
x = x( :, [ 2 1 3 ] );
