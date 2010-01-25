function [v, d] = scirunnrrd_pca(scirunnrrd)
% SCIRUNNRRD_PCA Principal Principal Component Analysis of the selected
% points in a SCI NRRD segmentation mask
%
% [V, D] = SCIRUNNRRD_PCA(X)
%
%   This function computes Principal Component Analysis (PCA) on the
%   collection of points in a SCI NRRD segmentation mask.
%
%   X is the SCI NRRD struct.
%
%   V is a matrix with the column eigenvectors ordered in decreasing order
%   of the corresponding eigenvalues in D.
%
%   Assuming a 3D volume, V is a (3, 3)-matrix and D a 3-vector.
%
%   Software applications developed at the University of Utah Scientific
%   Computing and Imaging (SCI) Institute, e.g. Seg3D, internally use NRRD
%   volumes to store medical data.
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

% Copyright © 2009 University of Oxford
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

% extract linear indices of voxels in the segmentation
idx = find( scirunnrrd.data );

% get volume size (first dimension is dummy, e.g. [1    62    81   206])
sz = size( scirunnrrd.data );

% convert linear index to multiple subscripts
[ix, iy, iz] = ind2sub( sz( 2:end ), idx );

% convert indices to real world coordinates
x = scirunnrrd_index2world( [ ix, iy, iz ], scirunnrrd.axis );

% compute PCA
[ v, d ] = pts_pca( x' );
