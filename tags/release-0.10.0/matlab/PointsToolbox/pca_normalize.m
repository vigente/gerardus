function [x, xm, v, d] = pca_normalize(x)
% PCA_NORMALIZE  Normalize a set of points so that its centroid is zero,the
% Principal Components correspond to the Cartesian axes, and the variance
% (eigenvalue) along each Principal Component is one
%
% Y = pca_normalize(X)
%
%   X is a matrix where each row corresponds to a point and each column to
%   a coordinate.
%
%   Y is a matrix that corresponds to X normalized in a PCA sense. That is,
%   X is centered on zero, rotated so that its Principal Components
%   coincide with the Cartesian axes, and scaled so that the variance along
%   each axis is one.
%
% [Y, XM, V, D] = pca_normalize(X)
%
%   XM is the centroid of X.
%
%   V is matrix of Principal Components (eigenvectors) of X.
%
%   D is a vector with the variance (eigenvalues) along the Principal
%   Compoments of X.
%
%   Note that you can recover the original data set running
%
%     y = y .* repmat(sqrt(d)', size(x, 1), 1);  % unscale
%     y = y * v';                                % undo rotation
%     y = y + repmat(xm, size(x, 1), 1);         % undo centering
%
% See also: pts_cn.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

% check arguments
narginchk(1, 1);
nargoutchk(0, 4);

% compute PCA on the points so that we get the main axes of the data
[v, d] = pts_pca(x');

% center the points
xm = mean(x);
x = x - repmat(xm, size(x, 1), 1);

% rotate the points to align their eigenvectors with the XYZ axes (the main
% eigenvectors with the X axis, the second with Y, and the third with Z)
x = x * v;

% equalise the axes so that the new eigenvalues are all one
x = x ./ repmat(sqrt(d)', size(x, 1), 1);
