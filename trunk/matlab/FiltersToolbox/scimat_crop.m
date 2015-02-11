function scimat = scimat_crop(scimat, from, to)
% SCIMAT_CROP  Crop a SCIMAT image or segmentation volume.
%
% SCIMAT2 = scimat_crop(SCIMAT, FROM, TO)
%
%   SCIMAT is a struct with a metainfo-enriched image or segmentation (see
%   "help scimat" for details).
%
%   FROM, TO are vectors with the index coordinates that define the
%   cropping box. The number of indices depend on the dimensionality of the
%   image. For example, for a 3D volume, FROM=[2 3 7], TO=[15 20 22].
%
%   If FROM or TO don't have all the dimensions, those dimensions are not
%   cropped. E.g. in an image with size [256 512 128 10],
%
%     FROM=[50 75];
%     TO=[100 89];
%
%   is the same as
%
%     FROM=[50 75 1 1];
%     TO=[100 89 128 10];
%
%   SCIMAT2 is the cropped volume.
%
% Example:
%
%   % create example scimat struct
%   scimat = scimat_im2scimat(zeros(256, 512, 128), [2 3 5], [-3.0 7.0 11.0]);
%
%   % crop struct
%   scimat2 = scimat_crop(scimat, [50 75], [100 89]);

% Authors: Ramon Casero <rcasero@gmail.com>, Benjamin Villard <b.016434@gmail.com>
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
narginchk(3, 3);
nargoutchk(0, 1);

% we only work with up to 4D images
MAXDIM = 4;
if (ndims(scimat.data) > MAXDIM)
    error('Image in SCIMAT must be 1D, 2D, 3D or 4D')
end

% FROM and TO can be shorter than the number of dimensions of the image,
% but not longer
if (length(from) > ndims(scimat.data))
    error('FROM has more dimensions than the image')
end
if (length(to) > ndims(scimat.data))
    error('TO has more dimensions than the image')
end

% default: extend "from" and "to" with dummy dimensions, if necessary, so
% that they are 4D
% E.g. in an image with size [256 512 128 10],
%
%     FROM=[50 75];
%     TO=[100 89];
%
%   is the same as
%
%     FROM=[50 75 1 1];
%     TO=[100 89 128 10];
%
% In addition, a 2D image with size [256 512] can be cropped as if it were
% a 4D image with size [256 512 1 1]
sz = size(scimat.data);
sz(length(sz)+1:MAXDIM) = 1;
from(length(from)+1:MAXDIM) = 1;
to(length(to)+1:MAXDIM) = sz(length(to)+1:MAXDIM);

% crop the image
scimat.data = scimat.data(from(1):to(1), from(2):to(2), from(3):to(3), ...
    from(4):to(4));

% correct the metainformation in the scimat volume
for I = 1:ndims(scimat.data)
    
    % size
    scimat.axis(I).size = size(scimat.data, I);
    
    % "left" edge of first voxel
    scimat.axis(I).min = scimat.axis(I).min ...
        + (from(I) - 1) * scimat.axis(I).spacing;
    
end
