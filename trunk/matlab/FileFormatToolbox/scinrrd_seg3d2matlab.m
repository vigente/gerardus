function nrrd = scinrrd_seg3d2matlab(nrrd)
% SCINRRD_SEG3D2MATLAB  Correct dimensions of data loaded from SCI NRRD so
% that rows and columns will follow Matlab's coordinates
%
% NRRD = SCINRRD_SEG3D2MATLAB(NRRD)
%
%   Seg3D displays images with x-coordinates spanning columns, and
%   y-coordinates spanning rows. However, when segmentations are exported
%   to Matlab files, it swaps this convention. Thus, when you load a
%   segmentation mask in Matlab, coordinates are inverted: x-coords span
%   rows and y-coords span columns.
%
%   This function undoes the swapping, so that the image loaded in Matlab
%   will follow the convention x <-> columns, y <-> rows.
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
% Copyright © 2010-2012 University of Oxford
% Version: 0.1.1
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
error(nargchk(1, 1, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% axis info may be missing, due to a bug in Seg3D2. In that case, the image
% will have 3 dimensions
if (~isfield(nrrd, 'axis') || isempty(nrrd.axis) ...
        || ~isfield(nrrd.axis, 'size'))
    warning('Missing axis information')
    % indices to swap rows and columns
    idx = [2 1 3];

else
    
    % the volume may have a dummy dimension or not (if scinrrd_squeeze()
    % has been used)
    if (length([nrrd.axis.size]) == 4) % dummy dimension present
        % indices to swap rows and columns
        idx = [1 3 2 4];
    else % dummy dimension removed
        % indices to swap rows and columns
        idx = [2 1 3];
    end
    
    nrrd.axis = nrrd.axis(idx);
end

nrrd.data = permute(nrrd.data, idx);
