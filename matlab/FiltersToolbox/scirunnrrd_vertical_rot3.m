function [m, v] = scirunnrrd_vertical_rot3(im)
% VERTICAL_ROTMAT  Compute the 3D rotation matrix to make a 3D segmented
% object vertical
%
% [M, A] = VERTICAL_ROTMAT3(IM)
%
%   This function computes a rotation matrix and centroid so that the input
%   SCI NRRD segmentation mask can be rotated to make the object vertical.
%
%   IM is the SCI NRRD struct.
%
%   M is a 2-vector with the coordinates of the segmentation mask centroid.
%
%   A is a (3,3)-rotation matrix. To make the object outlined by the
%   segmentation mask vertical, first substract the centroid from their
%   coordinates, so that the segmentation is centered around (0,0,0).
%
%   Second, to follow the ITK (Insight Toolkit) convention, A is the
%   transformation from output to input voxel coordinates. So if you want
%   to transform the coordinates of the *input* voxels, you need to use the
%   inverse of the rotation matrix, i.e. the transpose A'.
%
%   Finally, translate the rotated segmentation to the centroid.
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
%

% Copyright © 2010 University of Oxford
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
error( nargchk( 1, 1, nargin, 'struct' ) );
error( nargoutchk( 0, 2, nargout, 'struct' ) );

% compute image size
sz = size( im.data );

% get linear indices of segmented voxels
idx = find( im.data );

% convert the linear indices to volume indices
[ix, iy, iz] = ind2sub( sz( 2:end ), idx );

% compute real world coordinates for those indices
x = scirunnrrd_index2world( [ ix, iy, iz ], im.axis );

% compute centroid
m = mean( x );

% compute principal components so that we know the orientation of the
% points
v = pts_pca( x' );

% if we use the eigenvector matrix to rotate the data, then the first
% eigenvector will be projected on the x-axis, i.e. the object will be
% turned horizontal. So we need to reorganise the vectors. What we are
% going to do is make the object vertical (1st eigenvector projects on
% z-axis), and also that the major horizontal axis is parallel to the
% screen (2nd eigenvector projects on x-axis)
v = v( :, [ 2 3 1 ] );

% to comply with the ITK convention, we need the inverse rotation, i.e. the
% transpose
v = v';
