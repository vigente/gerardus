function [v, d, m] = scinrrd_pca(nrrd)
% SCINRRD_PCA Principal Principal Component Analysis of the selected
% points in a SCI NRRD segmentation mask
%
% [V, D, M] = SCINRRD_PCA(NRRD)
%
%   This function computes Principal Component Analysis (PCA) on the
%   collection of points in a SCI NRRD segmentation mask.
%
%   NRRD is the SCI NRRD struct.
%
%   V is a matrix with the column eigenvectors ordered in decreasing order
%   of the corresponding eigenvalues in D.
%
%   M is a 3-vector with the coordinates of the dataset centroid. The
%   centroid is the centre of the coordinate system defined by the
%   eigenvectors, and the centre of rotation if the data is rotated.
%
%   Assuming a 3D volume, V is a (3, 3)-matrix and D a 3-vector.
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
error( nargchk( 1, 1, nargin, 'struct' ) );
error( nargoutchk( 0, 3, nargout, 'struct' ) );

% squeeze NRRD variables, if necessary
nrrd = scinrrd_squeeze( nrrd );

% extract linear indices of voxels in the segmentation
idx = find( nrrd.data );

% get volume size
sz = size( nrrd.data );

% convert linear index to multiple subscripts
[ir, ic, iz] = ind2sub( sz, idx );

% convert indices to real world coordinates
x = scinrrd_index2world( [ ir, ic, iz ], nrrd.axis );

% compute centroid
m = mean( x );

% compute PCA
[ v, d ] = pts_pca( x' );
