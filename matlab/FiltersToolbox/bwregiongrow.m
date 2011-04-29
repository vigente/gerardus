function bwregiongrow
% BWREGIONGROW  Region grow labelling of binary image from multiple
% seeds
%
% LAB = BWREGIONGROW(IM, TODO)
%
%   IM is a 2D matrix or 3D array with a multi-label segmentation. That is,
%   IM contains a background and several objects, each object represented
%   by all the connected voxels with the same label.
%
%   Numerical voxel values in IM are interpreted in the following way:
%
%     0:               background voxel, don't label
%     TODO:            this voxel needs to be labelled using the region
%                      grow algorithm
%     Any other value: this voxel is a seed, and its value will be
%                      propagated using the region grow algorithm
%
%   LAB has the same size as IM, and has the label values that partition
%   the original binary image IM. Each partition has a different label.
%   Partitions are computed with a region grow algorithm that expands the
%   labels from the seeds.
%
%   At each iteration, partitions grow 1 voxel until the whole IM is
%   labelled.
%
% LAB = BWREGIONGROW(..., RES, MAXITER)
%
%   RES is a 2-vector (in 2D) or 3-vector (in 3D) with the voxel size in
%   each dimension. By default, it is assumed that RES=[1, 1, 1]. Voxel
%   size is used to compute distances between voxels in the labelling
%   process.
%
%   MAXITER is a scalar to tell the algorithm to stop after a number of
%   region grow iterations. If MAXITER < 0, the algorithm iterates until
%   all TODO voxels have been labelled. By default, MAXITER = -1.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.2.0
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

error('MEX function not found')
