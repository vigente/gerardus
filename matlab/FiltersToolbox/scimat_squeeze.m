function scimat = scimat_squeeze(scimat, todouble)
% SCIMAT_SQUEEZE  Remove dummy dimension from SCIRUNNRRD struct so that it
% becomes a SCIMAT struct.
%
% SCIMAT = scimat_squeeze(SCIRUNNRRD)
%
%   When the application Seg3D saves images to Matlab format, it saves them
%   to a SCIRUNNRRD struct, which is basically SCIMAT with a dummy
%   dimension.
%
%   SCIMAT is a struct that contains a 3D volume. We use SCIMAT structs
%   widely in Gerardus, because that way we have the image data and
%   metainformation (e.g. voxel size) together in the same
%   variable. For details on SCIMAT structs, see "help scimat".
%
%   This function removes the dummy dimension.
%
% SCIMAT = scimat_squeeze(..., TODOUBLE)
%
%   TODOUBLE is a flag to convert the image data from uint8 to double. This
%   is useful for image processing, but the data will be 8 times larger. By
%   default, TODOUBLE=false.
%
% See also: scimat_load.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010-2014 University of Oxford
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
narginchk(1, 2);
nargoutchk(0, 1);

% defaults
if (nargin < 2 || isempty(todouble))
    todouble = false;
end

% remove dummy dimension
scimat.data = squeeze(scimat.data);
if (length(scimat.axis) > 3)
    scimat.axis = scimat.axis(2:end);
end

% convert data to double
if (todouble)
    scimat.data = double(scimat.data);
end
