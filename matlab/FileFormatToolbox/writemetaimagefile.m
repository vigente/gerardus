function writemetaimagefile(filename, img, resolution, offset)
% WRITEMETAIMAGEFILE  Write a MetaImage file (.mha) with both header and data
%
% WRITEMETAIMAGEFILE(FILENAME, IMG, RESOLUTION, DATA_TYPE, OFFSET)
%
%   FILENAME is the path and name of the file to be written, e.g.
%   'foo.mha'.
%
%   IMG is a 3D matrix that contains the image volume.
%
%   RESOLUTION is a 3-vector with the voxel size in the 3 directions.
%
%   OFFSET is a vector with the real world coordinates in metres (not index
%   coordinates) of the centre of the first voxel in the volume. For
%   example, OFFSET=[0.014552, 0.010486, 0.00142]. By default, 
%   OFFSET=[0 0 0].
%
%   See also: WriteMhaFile, WriteRawFile.

% Author(s): Ramon Casero <rcasero@gmail.com> and Vicente Grau
% Copyright © 2012 University of Oxford
% Version: 0.1.3
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
error(nargchk(2, 4, nargin, 'struct'));
error(nargoutchk(0, 0, nargout, 'struct'));

% defaults
if (nargin < 3 || isempty(resolution))
    resolution = [1.0 1.0 1.0];
end
if (nargin < 4 || isempty(offset))
    offset = [0.0 0.0 0.0];
end

% get image size
img_size = size(img);

% get pixel type
data_type = class(img);

% open file for writing
fid=fopen(filename, 'w');
if(fid<=0) 
    fprintf('Impossible to open file %s\n', filename);
end

ndims=numel(resolution);

if(ndims == 3)
    fprintf(fid, 'NDims = 3\n');

    fprintf(fid, 'DimSize = %d %d %d\n', img_size(2), img_size(1), img_size(3));

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

    fprintf(fid, 'Offset = %1.12e %1.12e %1.12e\n', ...
        offset(1), offset(2), offset(3));

    fprintf(fid, 'ElementSpacing = %1.12e %1.12e %1.12e\n', ...
        resolution(1), resolution(2), resolution(3));

elseif(ndims==4)
    fprintf(fid, 'NDims = 4\n');

    fprintf(fid, 'DimSize = %d %d %d %d\n', ...
        img_size(2), img_size(1), img_size(3), img_size(4));

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

    fprintf(fid, 'Offset = %1.12e %1.12e %1.12e\n', ...
        offset(1), offset(2), offset(3));

    fprintf(fid, 'ElementSpacing = %1.12e %1.12e %1.12e %1.12e\n', ...
        resolution(1), resolution(2), resolution(3), resolution(4));
       
end

fprintf(fid, 'ElementByteOrderMSB = False\n');

fprintf(fid, 'ElementDataFile = LOCAL\n');

% write image data to file
fwrite(fid, img, data_type);

fclose(fid);
