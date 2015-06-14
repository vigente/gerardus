function c = imcmass(im)
% IMCMASS  Centre of mass of an image.
%
% IMCMASS computes the centre of mass of a grayscale image weighting the
% pixel coordinates by the pixel intensities.
%
% C = IMCMASS(IM)
%
%   IM is a grayscale 2D image.
%
%   C is the center of mass of the image in [row, col] format.

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

if (size(im, 3) > 1)
    error('IM must be a grayscale image')
end

% grid for image coordinates
[grow, gcol] = ndgrid(1:size(im, 1), 1:size(im, 2));

% centre of mass
c = single(im(:));
c = sum([c .* grow(:), c .* gcol(:)]) / sum(c);
