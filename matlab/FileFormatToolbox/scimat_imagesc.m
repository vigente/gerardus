function h = scimat_imagesc(scimat, varargin)
% SCIMAT_IMAGESC  Display 2D image in real world coordinate axes.
%
% SCIMAT_IMAGESC displays the 2D image in SCIMAT.data scaling the axes to
% real world coordinates using the metadata in SCIMAT.axis.
%
% This function is a thin interface to imagesc(), with the same syntax,
% except that it saves you from having to compute X, Y for 
% imagesc(X, Y, SCIMAT.data).
%
% SCIMAT_IMAGESC(SCIMAT)
%
%   SCIMAT is a struct with an image and metadata (see "help scimat").
%
% SCIMAT_IMAGESC(SCIMAT, ...)
% 
%   See "help image" for extra options.
%
% See also: image, imagesc.

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
nargoutchk(0, 1);

% image size
N = size(scimat.data, 1);
M = size(scimat.data, 2);

if (size(scimat.data, 3) > 1 || size(scimat.data, 4) > 1)
    error('SCIMAT must contain a 2D image, and no more than one frame')
end

% real world coordinates of the first and last voxels
xmin = scimat_index2world([1, 1], scimat);
xmax = scimat_index2world([N M], scimat);

% display image
if (nargout == 0)
    imagesc(...
        linspace(xmin(1), xmax(1), M), ...
        linspace(xmin(2), xmax(2), N), ...
        squeeze(scimat.data), ...
        varargin{:} ...
        )
else
        h = imagesc(...
        linspace(xmin(1), xmax(1), M), ...
        linspace(xmin(2), xmax(2), N), ...
        squeeze(scimat.data), ...
        varargin{:} ...
        );
end
