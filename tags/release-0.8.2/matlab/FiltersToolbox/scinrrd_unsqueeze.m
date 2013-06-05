function nrrd = scinrrd_unsqueeze(nrrd, touint8)
% SCINRRD_UNSQUEEZE  Add dummy dimension and convert data to uint8 data
% type
%
% NRRD = SCINRRD_UNSQUEEZE(NRRD, TOUINT8)
%
%   NRRD is an SCI NRRD image volume struct. Even if the data is 3D, Seg3D
%   adds a dummy dimension to the struct. This dummy dimension is removed
%   by function SCINRRD_SQUEEZE() in this toolbox for covenience. Function
%   SCINRRD_UNSQUEEZE() does the opposite and recovers the dummy dimension.
%
%   TOUINT8 is a flag to convert the image data from double to uint8. This
%   will make the volume 8 times smaller. By default, TOUINT8=false.
%
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
% Copyright © 2010 University of Oxford
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

% check arguments
error( nargchk( 1, 2, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

% defaults
if ( nargin < 2 || isempty( touint8 ) )
    touint8 = false;
end

% add dummy dimension
if ( ndims( nrrd.data ) == 3 )
    nrrd.data = reshape( nrrd.data, [1 size(nrrd.data)] );
end
if ( length( nrrd.axis ) == 3 )
    % add a new element at the beginning of nrrd.axis that will be used for
    % the dummy dimension
    nrrd.axis = [nrrd.axis(1); nrrd.axis];
    % add dummy dimension values
    nrrd.axis(1).size = 1;
    nrrd.axis(1).spacing = NaN;
    nrrd.axis(1).min = 0;
    nrrd.axis(1).max = 1;
    nrrd.axis(1).center = 1;
    nrrd.axis(1).label = 'axis 0';
    nrrd.axis(1).unit = 'no unit';
end

% convert data to uint8
if ( touint8 )
    nrrd.data = uint8( nrrd.data );
end
