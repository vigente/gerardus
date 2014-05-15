% test_pca_normalize.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013-2014 University of Oxford
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

% load segmentation
scimat = scimat_load('data/008-lvhull-downsampled-4.mha');

% compute surface mesh from segmentation
opt = .004;
method = 'simplify';

tic
[x, tri] = v2s(single(scimat.data), 1, opt, method);
tri = tri(:, 1:3);
toc

% plot mesh
subplot(2, 1, 1)
hold off
plotmesh(x, tri)
axis equal

% normalize mesh vertices
[y, xm, v, d] = pca_normalize(x);

% plot normalized points
subplot(2, 1, 2)
plotmesh(y, tri)
axis equal

% recover the original data set
y = y .* repmat(sqrt(d)', size(x, 1), 1);
y = y * v';
y = y + repmat(xm, size(x, 1), 1);

max(abs(x(:) - y(:)))
