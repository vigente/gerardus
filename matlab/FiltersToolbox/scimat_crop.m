function scimat = scimat_crop(scimat, from, to)
% SCIMAT_CROP  Crop a SCIMAT image or segmentation volume of any dimension.
%
% SCIMAT2 = scimat_crop(SCIMAT, FROM, TO)
%
%   SCIMAT is a struct with a metainfo-enriched image or segmentation (see
%   "help scimat" for details). The image can have any number of
%   dimensions.
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
%   If FROM=[] and TO=[], or they are not given as input argument, then the
%   image is not cropped.
%
%   SCIMAT2 is the output cropped volume.
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
% Version: 0.4.0
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
narginchk(1, 3);
nargoutchk(0, 1);

% defaults
if (nargin < 2)
    from = [];
end
if (nargin < 3)
    to = [];
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
% that they are the same dimension as the image
% E.g. in an image with size [256 512 128 10],
%
%     FROM=[50 75];
%     TO=[100 89];
%
%   is the same as
%
%     FROM=[50 75 1 1];
%     TO=[100 89 128 10];
sz = size(scimat.data);
sz(length(sz)+1:ndims(scimat.data)) = 1;
from(length(from)+1:ndims(scimat.data)) = 1;
to(length(to)+1:ndims(scimat.data)) = sz(length(to)+1:end);

% "bottom-left" coordinates of what is going to become the first voxel of
% the volume after cropping
%
% swap X and Y coordinates, so that they match the scimat convention for
% the axes: (r, c, s, f) <-> (y, x, y, t)
xmin = scimat_index2world(from, scimat);
xmin = xmin([2 1 3:end]) - [scimat.axis.spacing]/2;

% crop and correct metainformation in the scimat image
for I = 1:length(scimat.axis)
    
    % crop one dimension of the image. We are going to be shifting
    % circularly the dimensions of the image so that we always crop the
    % first dimension
    idx = cell(1, length(scimat.axis));
    idx{1} = from(I):to(I);
    idx(2:end) = {':'};
    scimat.data = scimat.data(idx{:});
    
    % size
    scimat.axis(I).size = size(scimat.data, 1);
    
    % "left" edge of first voxel
    scimat.axis(I).min = xmin(I);
    
    % circular shift of the dimensions of image so that we can crop along
    % the first dimension in the next iteration
    idx = 1:length(scimat.axis);
    scimat.data = permute(scimat.data, idx([2:end 1]));
    
end

% Note: the last circular shift of the image dimension has already returned
% the image to its original shape, so no need to do anything after the loop
