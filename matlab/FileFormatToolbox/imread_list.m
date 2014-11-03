function im = imread_list(pathstr, file, TOGRAY)
% IMREAD_LIST  Read list of image files into volume.
%
% IMREAD_LIST reads a list of image files, converts them to grayscale if
% desired, and collates them to form a 3D or 4D image volume.
%
% IM = IMREAD_LIST(PATHSTR, FILE)
%
%   PATHSTR is a string with the path to the files. All images must be the
%   same size and have the same number of channels.
%
%   FILE is an array of structs obtained with dir().
%
%   IM is the resulting image volume with indices:
%
%     IM(rows, columns, slice, channel)
%
% IM = IMREAD_LIST(..., TOGRAY)
%
%   TOGRAY is a boolean to convert colour images to grayscale. By default,
%   TOGRAY='false'.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

% check arguments
narginchk(2, 3);
nargoutchk(0, 1);

% defaults
if (nargin < 3 || isempty(TOGRAY))
    TOGRAY = false;
end

% number of files
N = length(file);

% no files to read
if (N == 0)
    im = [];
    return
end

% get number of rows and columns in the blockface volume, and whether
% images are colour. The easiest way is to just load the first image,
% instead of using imfinfo and then having to convert e.g. BitDepth to a
% Matlab class
im0 = imread([pathstr filesep file(1).name]);

% number of channels
nchan = size(im0, 3);

% initialize output volume
if (nchan == 1 || TOGRAY)
    im = zeros(size(im0, 1), size(im0, 2), N, class(im0));
else
    im = zeros(size(im0, 1), size(im0, 2), N, nchan, class(im0));
end
im(:, :, 1, :) = im0;

% loop image files
for I = 2:N
    
    % load image from file
    if (nchan == 1 || TOGRAY)
        im(:, :, I) = rgb2gray(imread([pathstr filesep file(I).name]));
    else
        im(:, :, I, :) = imread([pathstr filesep file(I).name]);
    end
    
end
