function scimat = scimat_seg3d2matlab(scimat)
% SCIMAT_SEG3D2MATLAB  Correct dimensions of data loaded from a *.mat file
% created with Seg3D so that it will follow Matlab's convention.
%
% SCIMAT = scimat_seg3d2matlab(SCIRUNNRRD)
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
%   SCIRUNNRRD is a struct loaded from a .mat file created with the
%   software application Seg3D.
%
%   SCIMAT is a struct with the order of rows and columns corrected so that
%   it follows Matlab's convention.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010-2014 University of Oxford
% Version: 0.2.0
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
narginchk(1, 1);
nargoutchk(0, 1);

% axis info may be missing, due to a bug in Seg3D2. In that case, the image
% will have 3 dimensions
if (~isfield(scimat, 'axis') || isempty(scimat.axis) ...
        || ~isfield(scimat.axis, 'size'))
    warning('Missing axis information')
    % indices to swap rows and columns
    idx = [2 1 3];

else
    
    % the volume may have a dummy dimension or not (if scimat_squeeze()
    % has been used)
    if (length([scimat.axis.size]) == 4) % dummy dimension present
        % indices to swap rows and columns
        idx = [1 3 2 4];
    else % dummy dimension removed
        % indices to swap rows and columns
        idx = [2 1 3];
    end
    
    scimat.axis = scimat.axis(idx);
end

scimat.data = permute(scimat.data, idx);
