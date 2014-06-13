function [v, d, m] = scimat_pca(scimat)
% SCIMAT_PCA Principal Principal Component Analysis of the selected points
% in a SCIMAT segmentation mask.
%
% [V, D, M] = scimat_pca(SCIMAT)
%
%   This function computes Principal Component Analysis (PCA) on the
%   collection of points in a SCIMAT segmentation mask.
%
%   SCIMAT is a struct with the segmentation (see "help scimat" for
%   details).
%
%   V is a matrix with the column eigenvectors ordered in decreasing order
%   of the corresponding eigenvalues in D.
%
%   M is a 3-vector with the coordinates of the dataset centroid. The
%   centroid is the centre of the coordinate system defined by the
%   eigenvectors, and the centre of rotation if the data is rotated.
%
%   Assuming a 3D volume, V is a (3, 3)-matrix and D a 3-vector.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010,2014 University of Oxford
% Version: 0.2.1
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
narginchk(1, 1) ;
nargoutchk(0, 3) ;

% squeeze SCIMAT variables, if necessary
scimat = scimat_squeeze(scimat);

% extract linear indices of voxels in the segmentation
idx = find(scimat.data);

% get volume size
sz = size(scimat.data);

% convert linear index to multiple subscripts
[ir, ic, iz] = ind2sub(sz, idx);

% convert indices to real world coordinates
x = scimat_index2world([ir, ic, iz], scimat);

% compute centroid
m = mean(x);

% compute PCA
[v, d] = pts_pca(x');
