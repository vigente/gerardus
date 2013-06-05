% test_im2imat.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010-2013 University of Oxford
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

% minimal image volume with 2 voxels in each direction, and all voxels = 1
% so that we can check that the function is computing distances between
% voxels right

d = im2dmatrix(ones(2,2));
full(d)

%          0    1.0000    1.0000    1.4142
%     1.0000         0    1.4142    1.0000
%     1.0000    1.4142         0    1.0000
%     1.4142    1.0000    1.0000         0

d = im2dmatrix(ones(2,2,2));
full(d)

%          0    1.0000    1.0000    1.4142    1.0000    1.4142    1.4142    1.7321
%     1.0000         0    1.4142    1.0000    1.4142    1.0000    1.7321    1.4142
%     1.0000    1.4142         0    1.0000    1.4142    1.7321    1.0000    1.4142
%     1.4142    1.0000    1.0000         0    1.7321    1.4142    1.4142    1.0000
%     1.0000    1.4142    1.4142    1.7321         0    1.0000    1.0000    1.4142
%     1.4142    1.0000    1.7321    1.4142    1.0000         0    1.4142    1.0000
%     1.4142    1.7321    1.0000    1.4142    1.0000    1.4142         0    1.0000
%     1.7321    1.4142    1.4142    1.0000    1.4142    1.0000    1.0000         0

% check averaging of voxel values

d = im2dmatrix([1 2; 3 4]);
full(d)
