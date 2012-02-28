function nrrd = scinrrd_save(file, nrrd, touint8, v73)
% SCINRRD_SAVE  Save a NRRD struct to a Matlab/MetaImage format that can be
% imported by Seg3D
%
% SCINRRD_SAVE(FILE, NRRD)
%
%   This function formats the NRRD volume and saves it to a .mat or .mha
%   file that can be imported as a segmentation or opened as a volume by
%   Seg3D.
%
%   FILE is a string with the path and name of the .mat/.mha file.
%
%   NRRD is the SCI NRRD struct.
%
% SCINRRD_SAVE(FILE, NRRD, TOUINT8)
%
%   TOUINT8 is a flag to convert the image data from double to uint8. This
%   will make the volume 8 times smaller. By default, TOUINT8=false.
%
% SCINRRD_SAVE(FILE, NRRD, TOUINT8, V73)
%
%   V73 is a flag to save the data to Matlab's HDF5-based MAT file format
%   v7.3. The reason is that if you volume is >= 2GB, Matlab skips saving
%   it and gives the warning message:
%
%      "Warning: Variable 'scirunnrrd' cannot be saved to a MAT-file whose
%      version is older than 7.3.
%      To save this variable, use the -v7.3 switch.
%      Skipping..."
%
%   This is a limitation of the v7 format. Read more here
%
%      http://www.mathworks.com/help/techdoc/rn/bqt6wls.html
%
%   However, note that Seg3D 1.x cannot read v7.3 files. By default,
%   V73=false.
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
% Version: 0.3.2
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
error(nargchk(2, 4, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% defaults
if (nargin < 3 || isempty(touint8))
    touint8 = false;
end
if (nargin < 4 || isempty(v73))
    v73 = false;
end

% get extension of output filename
[~, ~, ext] = fileparts(file);

switch lower(ext)
    case '.mat'
        % make x-,y-coordinates compatible with the Seg3D convention
        nrrd = scinrrd_seg3d2matlab(nrrd);
        
        % add dummy dimension, if necessary, and convert data to uint8, if
        % requested
        nrrd = scinrrd_unsqueeze(nrrd, touint8);
        
        % save data
        if (v73)
            save(file, 'scirunnrrd', '-v7.3');
        else
            save(file, 'scirunnrrd');
        end
    case '.mha'
        % swap rows and columns so that we have x-coordinates in the first
        % column, and y-coordinates in the second column, as expected by
        % the MetaImage format
        nrrd.data = permute(nrrd.data, [2 1 3:ndims(nrrd.data)]);
        
        % save data, doing the same permutation of the axis values
        writemetaimagefile(file, nrrd.data, ...
            [nrrd.axis([2 1 3]).spacing], ...
            [nrrd.axis([2 1 3]).min]+[nrrd.axis([2 1 3]).spacing]/2);
        
    otherwise
        error('Unrecognised output file format')
end
