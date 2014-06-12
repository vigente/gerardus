function scimat = scimat_unsqueeze(scimat, touint8)
% SCIMAT_UNSQUEEZE  Add dummy dimension to SCIMAT struct so that it becomes
% a SCIRUNNRRD struct.
%
% SCIRUNNRRD = scimat_unsqueeze(SCIMAT)
%
%   SCIMAT is a struct that contains a 3D volume. We use SCIMAT structs
%   widely in Gerardus, because that way we have the image data and
%   metainformation (e.g. voxel size) together in the same
%   variable. For details on SCIMAT structs, see "help scimat".
%
%   The application Seg3D can read Matlab files that contain a SCIRUNNRRD
%   struct, which is basically SCIMAT with a dummy dimension. This function
%   adds that dummy dimension.
%
% SCIRUNNRRD = scimat_unsqueeze(..., TOUINT8)
%
%   TOUINT8 is a flag to convert the image data from double to uint8. This
%   will make the volume 8 times smaller. By default, TOUINT8=false.
%
% See also: scimat_write.

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
if (nargin < 2 || isempty(touint8))
    touint8 = false;
end

% add dummy dimension
if (ndims(scimat.data) == 3)
    scimat.data = reshape(scimat.data, [1 size(scimat.data)]);
end
if (length(scimat.axis) == 3)
    % add a new element at the beginning of scimat.axis that will be used for
    % the dummy dimension
    scimat.axis = [scimat.axis(1); scimat.axis];
    % add dummy dimension values
    scimat.axis(1).size = 1;
    scimat.axis(1).spacing = NaN;
    scimat.axis(1).min = 0;
    scimat.axis(1).max = 1;
    scimat.axis(1).center = 1;
    scimat.axis(1).label = 'axis 0';
    scimat.axis(1).unit = 'no unit';
end

% convert data to uint8
if (touint8)
    scimat.data = uint8(scimat.data);
end
