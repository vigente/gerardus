function nrrd = scinrrd_load(file)
% SCINRRD_LOAD  Load a NRRD struct saved to Matlab format
%
% NRRD = SCINRRD_LOAD(FILE)
%
%   This function loads the NRRD volume, removes the dummy dimension and
%   corrects the row/column order so that the Seg3D image matches Matlab's.
%
%   FILE is a string with the path and name of the .mat file that contains
%   the NRRD volume.
%
%    NRRD is the SCI NRRD struct.
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

% check arguments
error( nargchk( 1, 1, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

% load data
nrrd = load(file);

% rename NRRD volume for convenience
nrrd = nrrd.scirunnrrd;

% remove dummy dimension
nrrd = scinrrd_squeeze(nrrd);

% correct x-,y-coordinates
nrrd = scinrrd_seg3d2matlab(nrrd);
