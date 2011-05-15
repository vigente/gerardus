function [m, v] = scinrrd_vertical_orientation_pca_basis(nrrd)
% SCINRRD_VERTICAL_ORIENTATION_PCA_BASIS  Compute a basis from the
% Principal Components of a set of voxels, such that the vertical axis is
% assigned to the maximum variability, and the basis is right-hand oriented
%
% [M, A] = SCINRRD_VERTICAL_ORIENTATION_PCA_BASIS(NRRD)
%
%   This function computes the centroid of an SCI NRRD segmentation, and a
%   basis where the axes give the Principal Components of variability of
%   the voxels. The vertical axis is assigned  to the maximum variability,
%   and the basis is right-hand oriented
%
%   NRRD is the SCI NRRD struct.
%
%   M is a 3-vector with the coordinates of the segmentation mask centroid.
%
%   A is a (3,3)-matrix. The correspondence between vectors and Principal
%   Components is:
%
%     * A(:,1): Middle variability, i.e. the major axis of the "horizontal
%               " ellipse
%     * A(:,2): Smallest variability
%     * A(:,3): Maximum variability, i.e. the object is considered to be
%               "vertical" in its longest axis
%
%   Also, the basis A is right-hand oriented, i.e. the cross product
%   cross(A(:,1), A(:,2)) = A(:,3).
%
%   Note that you can use M, A to make the segmentation object vertical, if
%   you rotate the voxel coordinates with the transpose matrix A' around
%   the centroid M.
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

% Author: Ramon Casero.
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
error( nargchk( 1, 1, nargin, 'struct' ) );
error( nargoutchk( 0, 2, nargout, 'struct' ) );

% compute Principal Component Analysis of the segmented voxels
[v, d, m] = scinrrd_pca( nrrd );

% if we use the eigenvector matrix to rotate the data, then the first
% eigenvector would be projected on the x-axis, i.e. the object would be
% turned horizontal. So we need to reorganise the vectors. What we are
% going to do is make the object vertical (1st eigenvector projects on
% z-axis), and also that the major horizontal axis is parallel to the
% screen (2nd eigenvector projects on x-axis)
v = v( :, [ 2 3 1 ] );

% note that v is also valid if we take the negative of any vector. In fact,
% we cannot be sure that e.g. the vectors in v fulfill X x Y = Z or 
% X x Y = -Z (where "x" is cross product). In the latter case, we would
% have a vertical reflection if we use v' to rotate the data. To avoid it,
% we impose two constraints...

% ... we want the Z axis to be pointing up in all circumstances
v(:, 3) = v(:, 3) * sign(v(3, 3));

% ... and we make sure that X x Y = Z, i.e. v has the same orientation as
% the Cartesian system. The Cartesian system fulfills that 
% cross([1 0 0],[0 1 0]) = [0 0 1]
aux = cross(v(:, 1), v(:, 2));
v(:, 1) = v(:, 1) * sign(aux(3));

