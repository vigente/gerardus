function writemetaimagefile(filename, img, resolution, offset, orientation)
% WRITEMETAIMAGEFILE  Write a MetaImage file (.mha) with both header and
% data.
%
% WRITEMETAIMAGEFILE(FILENAME, IMG, RESOLUTION, DATA_TYPE, OFFSET, ORIENTATION)
%
%   FILENAME is the path and name of the file to be written, e.g.
%   'foo.mha'.
%
%   IMG is an array that contains the image. The index-coordinate
%   correspondence must be:
%
%     [row, col, slice, frame, channel] -> [x, y, z, t, channel]
%
%   which is different from the Matlab standard for images, where rows
%   correspond to the y-coordinate. This is also different from the
%   Gerardus "scimat" image specification (see "help scimat").
%
%   RESOLUTION is a 3-vector with the voxel size in the 3 directions. The
%   order must be [dx, dy, dz].
%
%   OFFSET is a vector with the real world coordinates in metres (not index
%   coordinates) of the centre of the first voxel in the volume. For
%   example, OFFSET=[0.014552, 0.010486, 0.00142]. By default, 
%   OFFSET=[0 0 0]. The order must be the same as in RESOLUTION.
%
%   ORIENTATION is a vector with a linearized rotation matrix. E.g. for 3
%   dimensions ORIENTATION=[X(1:3) Y(1:3) Z(1:3)] means that the rotation
%   matrix is 
%                      [X(1:3)]
%                      [Y(1:3)]
%                      [Z(1:3)]
%
% See also scimat.

% Author(s): Ramon Casero <rcasero@gmail.com> and Vicente Grau
% Copyright © 2012, 2015 University of Oxford
% Version: 0.2.2
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
narginchk(2, 5);
nargoutchk(0, 0);

% [row, col, slice, frame, channel]
if (ndims(img) < 2 || ndims(img) > 5)
    error('IMG dimension must be 2, 3, 4 or 5')
end

% number of physical dimensions (i.e. excluding channels because MetaImage
% treats them separately, and time frames because we haven't implemented
% writing time series yet)
if (nargin >= 3)
    D = length(resolution);
else
    D = min([ndims(img) 3]);
end

% defaults
if (nargin < 3 || isempty(resolution))
    resolution = ones(1, D);
end
if (nargin < 4 || isempty(offset))
    offset = zeros(1, D);
else
    if (length(offset) ~= D)
        error('OFFSET must have the same dimensions as RESOLUTION')
    end
end
if (nargin < 5 || isempty(orientation))
    orientation = eye(D);
    orientation = orientation(:)';
else
    if (length(orientation) ~= D^2)
        error('ORIENTATION must be a vector with D^2 elements, if D is the number of elements in RESOLUTION')
    end
end

% get image size
img_size = size(img);

% TODO: writing time series not implemented
if (length(img_size) >= 4 && img_size(4) > 1)
    error('Writing of time series not implemented')
end

% number of channels
if (length(img_size) < 5)
    nchannel = 1;
else
    nchannel = img_size(5);
end

% get pixel type
data_type = class(img);

% open file for writing
fid = fopen(filename, 'w');
if(fid <= 0) 
    error(['Could not open file: ' filename]);
end

% number of spatial physical dimensions
fprintf(fid, 'NDims = %d\n', D);

aux = num2cell(img_size(1:D));
fprintf(fid, ['DimSize =' repmat(' %d', 1, D) '\n'], aux{:});

fprintf(fid, 'ElementNumberOfChannels = %d\n', nchannel);

switch data_type
    case {'logical', 'uint8'}
        fprintf(fid, 'ElementType = MET_UCHAR\n');
    case 'int8'
        fprintf(fid, 'ElementType = MET_CHAR\n');
    case 'uint16'
        fprintf(fid, 'ElementType = MET_USHORT\n');
    case 'int16'
        fprintf(fid, 'ElementType = MET_SHORT\n');
    case 'uint32'
        fprintf(fid, 'ElementType = MET_UINT\n');
    case 'int32'
        fprintf(fid, 'ElementType = MET_INT\n');
    case 'single'
        fprintf(fid, 'ElementType = MET_FLOAT\n');
    case 'double'
        fprintf(fid, 'ElementType = MET_DOUBLE\n');
    otherwise
        error('Unrecognized data type')
end

aux = num2cell(offset);
fprintf(fid, ['Offset =' repmat(' %1.12e', 1, D) '\n'], aux{:});

aux = num2cell(orientation);
fprintf(fid, ['Orientation =' repmat(' %1.12e', 1, D^2) '\n'], aux{:});

aux = num2cell(resolution);
fprintf(fid, ['ElementSpacing =' repmat(' %1.12e', 1, D) '\n'], aux{:});

fprintf(fid, 'ElementByteOrderMSB = False\n');

fprintf(fid, 'ElementDataFile = LOCAL\n');

% rearrange dimensions so that channels go from 5th to 1st dimension
if (ndims(img) >= 5)
    img = permute(img, [5 1:4]);
end
        
% write image data to file
fwrite(fid, img, data_type);

fclose(fid);
