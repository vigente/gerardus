function nrrd = scinrrd_squeeze(nrrd, todouble)
% SCINRRD_SQUEEZE  Remove dummy dimension and convert data to double type
%
% NRRD = SCINRRD_SQUEEZE(NRRD, TODOUBLE)
%
%   NRRD is an SCI NRRD image volume struct. Even if the data is 3D, the
%   struct adds a dummy dimension. This function removes it, so that
%   indices are more intuitive (e.g., x-coordinates are data(:, , ) instead
%   of data( , :, , )
%
%   TODOUBLE is a flag to convert the image data from uint8 to double. This
%   is useful for image processing, but the data will be 8 times larger. By
%   default, TODOUBLE=false.
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

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010-2011 University of Oxford
% Version: 0.1.2
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
error(nargchk(1, 2, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% defaults
if (nargin < 2 || isempty(todouble))
    todouble = false;
end

% remove dummy dimension
nrrd.data = squeeze(nrrd.data);
if (length(nrrd.axis) > 3)
    nrrd.axis = nrrd.axis(2:end);
end

% convert data to double
if (todouble)
    nrrd.data = double(nrrd.data);
end
