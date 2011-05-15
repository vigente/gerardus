function [im, zp, xp, yp] = scinrrd_intersect_plane(nrrd, m, v, xg, yg, zg)
% SCINRRD_INTERSECT_PLANE  Compute intersection of a plane with an SCI
% NRRD image volume
%
% IM = SCINRRD_INTERSECT_PLANE(NRRD, M, V)
%
%   IM is an image that displays the intersection of the plane defined by
%   M, V with the image volume NRRD.
%
%   M is a 3-vector with the coordinates of a point contained in the plane.
%   
%   V is a 3-vector that represents a normalized vector orthogonal to the
%   plane.
%   
%   Together, M and V define a unique plane in 3D space.
%
% [IM, ZP, XP, YP] = SCINRRD_INTERSECT_PLANE(NRRD, M, V, XG, YG, ZG)
%
%   XG, YG, ZG are grid coordinates for NRRD computed with NDGRID. If not
%   provided, they will be computed internally. Providing them is just a
%   way of saving time if many intersection planes have to be computed.
%
%   ZP is a matrix with the z-coordinates of the pixels in IM.
%
%   Because the XP and YP matrices are easily computed (e.g. they are each
%   of the identical slices of the XG, YG volumes), and they are the same
%   for every plane, they are provided as output arguments after the ZP
%   matrix.
%
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
error( nargchk( 3, 6, nargin, 'struct' ) );
error( nargoutchk( 0, 4, nargout, 'struct' ) );

% remove dummy dimension of data if necessary
nrrd = scinrrd_squeeze( nrrd, true );

% this function has a singularity if the plane is vertical
if ( v(3) == 0 )
    error( 'This function doesn''t accept vertical planes' )
end

% generate 3D grid of coordinates if they are not already provided
if ( (nargin < 6) || isempty( xg ) || isempty( yg ) || isempty( zg ) )
    [ xg, yg, zg ] = scinrrd_ndgrid( nrrd );
end

% extract one slice from the grid for covenience to express the x-,
% y-coordinates of plane points
xp = xg(:,:,1);
yp = yg(:,:,1);

% A plane with orthogonal vector n that goes through the centroid m
% can be expressed by the formula
% 
% nx(x-mx) + ny(y-my) + nz(z-mz) = 0
% 
% We know the horizontal limits of the image. So in order to
% compute the intersection height of the rotated plane at each one
% of 2 vertices that delimitate the image horizontally, we use the
% derived expression
% 
% z = nx/nz(mx-x) + ny/nz(my-y) + mz
% 
% This expression is invalid when nz=0, i.e. the plane has been
% rotated to make it "vertical"

% compute the corresponding heights at the horizontal grid, i.e.
% the rotated plane
zp = v(1)/v(3)*(m(1)-xp) + v(2)/v(3)*(m(2)-yp) + m(3);

% compute intersection of image with the 2D plane
% (if you want to visualize the image as in Seg3D, you need to do
% 'axis xy')
% note: the swapping of xg and yg is because interpn requires the matrices
% in the same order as they are produced by ndgrid()
im = interpn( yg, xg, zg, nrrd.data, yp, xp, zp, 'linear', 0 );
