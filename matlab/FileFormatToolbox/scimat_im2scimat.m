function scimat = scimat_im2scimat(im, res, offset, rotmat)
% SCIMAT_IM2SCIMAT  Create SCIMAT struct from scratch.
%
% This function creates a struct with the scimat format (see "help scimat"
% for details).
%
% SCIMAT = scimat_im2scimat(IM, RES, OFFSET, ROTMAT)
%
%   IM is a Matlab array with the image or segmentation. IM can be of class
%   logical, (u)int8, (u)int16, (u)int32, (u)int64, single or double. IM
%   can have 3 to 4 dimensions.
%
%   RES is a vector with the voxel size. By default, RES is 1.0 in each
%   dimension.
%
%   OFFSET is a 3-vector with the coordinates of the *centre* of the first
%   voxel in the image.
%
%   Note that the order of RES and OFFSET is the same as the axes in IM.
%   E.g. for a 4D image
%
%     RES(1) --> rows     (y axis)
%     RES(2) --> columns  (x axis)
%     RES(3) --> slices   (z axis)
%     RES(4) --> frames   (t axis)
%
%   ROTMAT is a rotation matrix. By default, ROTMAT=[], which means that no
%   rotation is applied. This is a better solution than using an identity
%   matrix, e.g. because a 3D image could represent a 2D+t or 3D image,
%   which requires a (2,2)- or (3,3)-rotation matrix, respectively.
%
% See also: scimat.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011-2015 University of Oxford
% Version: 0.3.0
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
narginchk(1, 4);
nargoutchk(0, 1);

% defaults
if (nargin < 2 || isempty(res))
    res = ones(1, ndims(im));
end
if (nargin < 3 || isempty(offset))
    offset = zeros(1, ndims(im));
end
if (nargin < 4 || isempty(rotmat))
    rotmat = [];
end

% create NRRD struct

% data volume
scimat.data = im;

% loop some of the fields
for I = 1:3
    
    % data volume size
    scimat.axis(I).size = size(im, I);
    
    % image resolution
    scimat.axis(I).spacing = res(I);
    
    % left edge of first voxel
    scimat.axis(I).min = offset(I) - res(I) / 2;
    
end

% we need scimat.axis to be a column vector
scimat.axis = scimat.axis';

% add rotation matrix
scimat.rotmat = rotmat;
