function scirunnrrd = scinrrd_save(file, scirunnrrd, touint8)
% SCINRRD_SAVE  Save a NRRD struct to a Matlab format that can be imported
% by Seg3D
%
% SCINRRD_SAVE(FILE, NRRD)
%
%   This function formats the NRRD volume and saves it to a .mat file that
%   can be imported as a segmentation or opened as a volume by Seg3D.
%
%      Note: If you volume is very large, Matlab gives a spurious warning
%      message
%
%      "Warning: Variable 'scirunnrrd' cannot be saved to a MAT-file whose
%      version is older than 7.3.
%      To save this variable, use the -v7.3 switch.
%      Skipping..."
%
%      This has nothing to do with versions, it's some kind of memory bug.
%
%   FILE is a string with the path and name of the .mat file.
%
%   NRRD is the SCI NRRD struct.
%
% SCINRRD_SAVE(FILE, NRRD, TOUINT8)
%
%   TOUINT8 is a flag to convert the image data from double to uint8. This
%   will make the volume 8 times smaller. By default, TOUINT8=false.
%
% NRRD = SCINRRD_SAVE(...)
%
%   NRRD as an output argument gives the reformatted data volume actually
%   saved.
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
%

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010-2011 University of Oxford
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
error(nargchk(2, 3, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% defaults
if (nargin < 3 || isempty(touint8))
    touint8 = false;
end

% make x-,y-coordinates compatible with the Seg3D convention
scirunnrrd = scinrrd_seg3d2matlab(scirunnrrd);

% add dummy dimension, if necessary, and convert data to uint8
scirunnrrd = scinrrd_unsqueeze(scirunnrrd, touint8);

% save data
save(file, 'scirunnrrd');
