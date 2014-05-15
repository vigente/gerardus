function WriteMhaFile(filename, img_size, resolution, data_type, offset)
% WRITEMHAFILE  Write the header part of a MetaImage file (.mha)
%
% WRITEMHAFILE(FILENAME, IMG_SIZE, RESOLUTION, DATA_TYPE, OFFSET)
%
%   FILENAME is the path and name of the file to be written, e.g.
%   'foo.mha'.
%
%   IMG_SIZE is a 3-vector with the size of the output volume.
%
%   RESOLUTION is a 3-vector with the voxel size in the 3 directions.
%
%   DATA_TYPE is a string with the data type as given my Matlab, e.g.
%   'uint8', 'short', 'uint16'.
%
%   OFFSET is a vector with the real world coordinates in metres (not index
%   coordinates) of the first voxel in the volume. For example,
%   OFFSET=[0.014552, 0.010486, 0.00142]. By default, OFFSET=[0 0 0].
%
%   See also: WriteRawFile, writemetaimagefile.

% Author(s): Ramon Casero <rcasero@gmail.com> and Vicente Grau
% Copyright © 2009-2012 University of Oxford
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


% check arguments
error(nargchk(4, 5, nargin, 'struct'));
error(nargoutchk(0, 0, nargout, 'struct'));

% defaults
if (nargin < 5 || isempty(offset))
    offset = [0.0 0.0 0.0];
end

[path, name] = fileparts(filename);

% open file for writing
fid=fopen(filename, 'w');
if(fid<=0) 
    fprintf('Impossible to open file %s\n', filename);
end

ndims=numel(resolution);

if(ndims == 3)
    fprintf(fid, 'NDims = 3\n');

    fprintf(fid, 'DimSize = %d %d %d\n', img_size(2), img_size(1), img_size(3));

    if(strcmp(data_type, 'uint8'))
        fprintf(fid, 'ElementType = MET_UCHAR\n');
    elseif(strcmp(data_type, 'short'))
        fprintf(fid, 'ElementType = MET_SHORT\n');
    elseif (strcmp(data_type, 'uint16'))
        fprintf(fid, 'ElementType = MET_USHORT\n');
    else
        error('Not implemented data type')
    end

    fprintf(fid, 'Offset = %1.6f %1.6f %1.6f\n', ...
        offset(1), offset(2), offset(3));

    fprintf(fid, 'ElementSpacing = %1.12f %1.12f %1.12f\n', resolution(1), resolution(2), resolution(3));

elseif(ndims==4)
    fprintf(fid, 'NDims = 4\n');

    fprintf(fid, 'DimSize = %d %d %d %d\n', img_size(2), img_size(1), img_size(3), img_size(4));

    if(strcmp(data_type, 'uint8'))
        fprintf(fid, 'ElementType = MET_UCHAR\n');
    elseif(strcmp(data_type, 'short'))
        fprintf(fid, 'ElementType = MET_SHORT\n');
    elseif (strcmp(data_type, 'uint16'))
        fprintf(fid, 'ElementType = MET_USHORT\n');
    else
        error('Not implemented data type')
    end

    fprintf(fid, 'Offset = %1.6f %1.6f %1.6f\n', ...
        offset(1), offset(2), offset(3));

    fprintf(fid, 'ElementSpacing = %1.12f %1.12f %1.12f %1.12f\n', resolution(1), resolution(2), resolution(3), resolution(4));
       
end

fprintf(fid, 'ElementByteOrderMSB = False\n');

fprintf(fid, 'ElementDataFile = %s\n', [name, '.raw']);

fclose(fid);
