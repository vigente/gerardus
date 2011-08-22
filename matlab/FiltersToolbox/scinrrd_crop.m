function nrrd = scinrrd_crop(nrrd, from, to)
% SCINRRD_CROP  Crop an SCI NRRD image or segmentation volume
%
% NRRD2 = SCINRRD_CROP(NRRD, FROM, TO)
%
%   NRRD is an SCI NRRD image or segmentation volume.
%
%   FROM, TO are vectors with the index coordinates that define the
%   cropping box. For example, FROM=[2 3 7], TO=[15 20 22].
%
%   NRRD2 is the cropped volume.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.0.1
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
error(nargchk(3, 3, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% if we have 2-vectors, make them 3-vectors for uniformity
if (length(from) == 2)
    from(3) = 1;
end
if (length(to) == 2)
    to(3) = 1;
end

% crop the data volume
nrrd.data = nrrd.data(from(1):to(1), from(2):to(2), from(3):to(3));

% correct the metainformation in the nrrd volume
for I = 1:3
    
    % size
    nrrd.axis(I).size = size(nrrd.data, I);
    
    % "left" edge of first voxel
    nrrd.axis(I).min = nrrd.axis(I).min ...
        + (from(I) - 1) * nrrd.axis(I).spacing;
    
    % "left" edge of last voxel
    nrrd.axis(I).max = nrrd.axis(I).min ...
        + (nrrd.axis(I).size - 1) * nrrd.axis(I).spacing;
    
end
