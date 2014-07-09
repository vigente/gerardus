function im = blockface_imread_gray(file_expr)
% BLOCKFACE_IMREAD_GRAY  Read multiple image files into grayscale volume.
%
% blockface_imread_gray reads a list of image files, converts them to
% grayscale if necessary, and collates them to form a 3D image volume.
%
% IM = blockface_imread_gray(FILE_EXPR)
%
%   FILE_EXPR is a string with an expression passed to dir() to list the
%   blockface image files, e.g. '/data/blockface/*.bmp' or
%   'C:\data\blockface\*.bmp'.
%
%   IM is the resulting image volume.

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
narginchk(1, 1);
nargoutchk(0, 1);

% list of blockface files
file = dir(file_expr);

% directory with the files
pathstr = fileparts(file_expr);

% number of files
N = length(file);

if (N == 0)
    im = [];
    return
end

% load the first image, so that we know the number of rows and columns in
% the blockface volume
im = imread([pathstr filesep file(1).name]);

% initialize matrix to keep the blockface volume
im = zeros(size(im, 1), size(im, 2), N, class(im));

% loop image files
for I = 1:N
    
    % load image from file
    imaux = imread([pathstr filesep file(I).name]);
    
    % convert to grayscale if image is in colour
    if (size(imaux, 3) == 3)
        imaux = rgb2gray(imaux);
    end
    
    % add to the volume
    im(:, :, I) = imaux;
    
end
