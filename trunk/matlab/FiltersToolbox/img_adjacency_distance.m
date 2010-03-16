function a = img_adjacency_distance(im)
% ADJACENCY_DISTANCE  Convert image to adjacency-distance list in sparse matrix form
%
% A = IMG_ADJACENCY_DISTANCE(IM)
%
%   IM is an image volume with dimensions (R, C, S).
%
%   A is a sparse matrix with dimensions (R*C*S, R*C*S), where element
%   (i,j) is the mean intensity between voxels with linear indices i and j.
%
%   Voxels with an Inf intensity are skipped.

% Author: Ramon Casero <rcasero@gmail.com>
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

error('Compiled MEX function has not been found')
