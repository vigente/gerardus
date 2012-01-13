function WriteRawFile(filename, img, resolution, data_type)
% WRITERAWFILE  Write the binary part of a MetaImage file (.raw)
%
% WRITERAWFILE(FILENAME, IMG, RESOLUTION, DATA_TYPE)
%
%   FILENAME is a string with the path and name of the output file.
%
%   IMG is a 3D matrix that contains the image volume.
%
%   RESOLUTION is a 3-vector with the voxel size in the X, Y and Z
%   directions.
%
%   DATA_TYPE is a string with the data type, e.g. 'MET_UCHAR',
%   'MET_SHORT', 'MET_USHORT'...
%
%   See also: WRITEMHAFILE to write the header part.

% Author(s): Ramon Casero <rcasero@gmail.com> and Vicente Grau
% Copyright © 2009-2010 University of Oxford
% Version: 0.1.0
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


[pathstr,name, ext]=fileparts(filename);

fid=fopen(filename, 'w');
if(fid<=0) 
    printf('Impossible to open file %s\n', filename);
end

if(ndims(img) == 3)

    % permute rows and columns so that rows are written to Y-coordinates and
    % columns are written to X-coordinates
    img = permute(img, [2 1 3]);

    fprintf(fid, 'NDims = 3\n');

    fprintf(fid, 'DimSize = %d %d %d\n', size(img,1), size(img,2), size(img,3));

    if(strcmp(data_type, 'char') || strcmp(data_type, 'uint8'))
        fprintf(fid, 'ElementType = MET_UCHAR\n');
    elseif(strcmp(data_type, 'short'))
        fprintf(fid, 'ElementType = MET_SHORT\n');
    elseif (strcmp(data_type, 'uint16'))
        fprintf(fid, 'ElementType = MET_USHORT\n');
    else
        error('Not implemented data type')
    end

    fprintf(fid, 'ElementSpacing = %1.12f %1.12f %1.12f\n', resolution(1), resolution(2), resolution(3));

elseif(ndims(img)==4)

    % permute rows and columns so that rows are written to Y-coordinates and
    % columns are written to X-coordinates
    img = permute(img, [2 1 3 4]);

    fprintf(fid, 'NDims = 4\n');

    fprintf(fid, 'DimSize = %d %d %d %d\n', size(img,1), size(img,2), size(img,3), size(img,4));

    if(strcmp(data_type, 'char') || strcmp(data_type, 'uint8'))
        fprintf(fid, 'ElementType = MET_UCHAR\n');
    elseif(strcmp(data_type, 'short'))
        fprintf(fid, 'ElementType = MET_SHORT\n');
    elseif (strcmp(data_type, 'uint16'))
        fprintf(fid, 'ElementType = MET_USHORT\n');
    else
        error('Not implemented data type')
    end

    fprintf(fid, 'ElementSpacing = %1.12f %1.12f %1.12f %1.12f\n', resolution(1), resolution(2), resolution(3), resolution(4));
       

elseif(ndims(img)==2)
    
    % permute rows and columns so that rows are written to Y-coordinates and
    % columns are written to X-coordinates
    img = permute(img, [2 1]);

    fprintf(fid, 'NDims = 2\n');

    fprintf(fid, 'DimSize = %d %d \n', size(img,1), size(img,2));

    if(strcmp(data_type, 'char') || strcmp(data_type, 'uint8'))
        fprintf(fid, 'ElementType = MET_UCHAR\n');
    elseif(strcmp(data_type, 'short'))
        fprintf(fid, 'ElementType = MET_SHORT\n');
    elseif (strcmp(data_type, 'uint16'))
        fprintf(fid, 'ElementType = MET_USHORT\n');
    else
        error('Not implemented data type')
    end

    fprintf(fid, 'ElementSpacing = %1.12f %1.12f\n', resolution(1), resolution(2));
    
else
    
    error('Number of dimensions in image must be in [2, 4]')
       
end

fprintf(fid, 'ElementByteOrderMSB = False\n');

fprintf(fid, 'ElementDataFile = %s\n', strcat(name, '.raw'));

fclose(fid);

fid=fopen(strcat(pathstr, '/', name, '.raw'), 'w');
if(fid<=0) 
    printf('Impossible to open file %s\n', img_filename);
end
fwrite(fid, img, data_type);
fclose(fid);
