function ra = scimat_scimat2imref(scimat)
% SCIMAT_SCIMAT2IMREF  Matlab's image reference frame from SCIMAT metadata.
%
% SCIMAT_SCIMAT2IMREF produces an image frame (i.e. imref2d object) as used
% by Matlab to display images in real world coordinates. For example,
% Matlab's imshow() function uses imref2d objects as inputs.
%
% Note that imref2d objects do not support rotations, so the rotation
% matrix SCIMAT.rotmat will be ignored with a warning if it's not the
% identity.
%
% RA = SCIMAT_SCIMAT2IMREF(SCIMAT)
%
%   SCIMAT is a struct with a 2D image and its metadata (see "help
%   scimat"). The image can have several frames (2D+t), but it cannot be
%   3D. The image can also be multi-channel (e.g. RGB).
%
%   RA is the output imref2d object (see "help imref2d").
%
% Example:
%
% % create a SCIMAT image with some offset and non-isotropic pixel spacing
% offset = [2, 3];
% spacing = [0.3, 0.6];
% scimat = scimat_im2scimat(imread('cameraman.tif'), spacing, offset);
% 
% % alternative 1: display it using scimat_imagesc
% subplot(2, 1, 1)
% hold off
% scimat_imagesc(scimat)
% axis equal
% 
% % alternative 2: display it using a reference frame from scimat_scimat2imref
% ra = scimat_scimat2imref(scimat);
% subplot(2, 1, 2)
% hold off
% imshow(scimat.data, ra);
%
%
% See also: scimat, imref2d, scimat_imagesc, imshow.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2015 University of Oxford
% Version: 0.1.0
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

% check input dimensions
if (size(scimat.data, 3) > 1)
    error('Implemented only for 2D images, not for 3D images')
end

% Matlab's image frame doesn't support rotation matrices
if (isfield(scimat, 'rotmat') && any(scimat.rotmat(:) ~= [1 0 0 1]'))
    warning('Matlab''s image frame doesn''t support rotation matrices. Ignoring rotation matrix')
end

% get image limits
ImageSize = [size(scimat.data, 1) size(scimat.data, 2)];
xmin = scimat_index2world([1, 1], scimat);
xmax = scimat_index2world(ImageSize, scimat);
XWorldLimits = [xmin(1) xmax(1)];
YWorldLimits = [xmin(2) xmax(2)];

ra = imref2d(ImageSize, XWorldLimits, YWorldLimits);
