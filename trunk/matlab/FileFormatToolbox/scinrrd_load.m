function nrrd = scinrrd_load(file)
% SCINRRD_LOAD  Load a SCINRRD struct from a Matlab, MetaImage or LSM file
%
% NRRD = SCINRRD_LOAD(FILE)
%
%   This function loads the NRRD volume, removes the dummy dimension and
%   corrects the row/column order so that if you do e.g.
%   imagesc(nrrd.data(:,:,50)), it produces the same image as the 50 axial
%   slice in Seg3D.
%
%   FILE is a string with the path and name of the .mat, .mha or .lsm file
%   that contains the 2D or 3D image.
%
%   NRRD is the SCI NRRD struct.
%
%     NRRD.axis(1) ==> rows
%     NRRD.axis(2) ==> columns
%     NRRD.axis(3) ==> slices
%
%   Note on SCI NRRD: Software applications developed at the University of
%   Utah Scientific Computing and Imaging (SCI) Institute, e.g. Seg3D,
%   internally use NRRD volumes to store medical data.
%
%   When label volumes (segmentation masks) are saved to a Matlab file
%   (.mat), they use a struct called "scirunnrrd" to store all the NRRD
%   information:
%
%   >>  scirunnrrd
%
%   scirunnrrd = 
%
%          data: [4-D uint8]
%          axis: [4x1 struct]
%      property: []

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010-2011 University of Oxford
% Version: 0.2.7
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
error(nargchk(1, 1, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% extract extension of filename in lower case
[pathstr, ~, ext] = fileparts(file);
ext = lower(ext);

switch lower(ext)
    
    case '.mat' % Matlab file
        % load data
        nrrd = load(file);
        
        % rename NRRD volume for convenience
        nrrd = nrrd.scirunnrrd;
        
        % remove dummy dimension
        nrrd = scinrrd_squeeze(nrrd);
        
        % correct x-,y-coordinates
        nrrd = scinrrd_seg3d2matlab(nrrd);
        
    case {'.mha', '.mhd'} % MetaImage file
        
        % open file to read
        fid=fopen(file, 'r');
        if (fid<=0)
            error(['Cannot open file: ' file])
        end
        
        % default values for the text header
        N = [];
        sz = [];
        data_type = [];
        offset = [];
        res = [];
        msb = [];
        rawfile = [];
        
        % process text header, and stop if we get to the raw data
        while 1
            % read text header line
            tline = fgetl(fid);
            
            % if end of text header stop reading
            % Warning: this only works if the first voxel=0. Otherwise, we
            % would read some voxels as part of a header line. To avoid it,
            % we now break after finding ElementDataFile
            if (tline(1) == 0), break, end
            
            % pointer to current position in header
            eoh0 = ftell(fid);
            
            % parse text header line
            
            % find location of "=" sign
            idx = strfind(tline, '=');
            if isempty(idx), break, end
            
            switch getlabel(tline, idx)
                case 'ndims'
                    N = getnumval(tline, idx);
                case 'dimsize'
                    sz = getnumval(tline, idx);
                case 'elementtype'
                    switch lower(strtrim(tline(idx+1:end)))
                        case 'met_uchar'
                            data_type = 'uint8';
                        case 'met_char'
                            data_type = 'int8';
                        case 'met_ushort'
                            data_type = 'uint16';
                        case 'met_short'
                            data_type = 'int16';
                        case 'met_uint'
                            data_type = 'uint32';
                        case 'met_int'
                            data_type = 'int32';
                        case 'met_float'
                            data_type = 'single';
                        case 'met_double'
                            data_type = 'double';
                        otherwise
                            error('Unrecognized ElementType')
                    end
                case 'offset'
                    offset = getnumval(tline, idx);
                case 'elementspacing'
                    res = getnumval(tline, idx);
                case 'elementbyteordermsb'
                    msb = strcmpi(strtrim(tline(idx+1:end)), 'true');
                case 'elementdatafile'
                    rawfile = strtrim(tline(idx+1:end));
                    break;
                case 'compresseddata'
                    if strcmpi(strtrim(tline(idx+1:end)), 'true')
                        error('Cannot read compressed MHA data')
                    end
                otherwise
                    warning(['Unrecognized line: ' tline])
            end
            
            % update position of end of header
            eoh = eoh0;

        end
        
        % the raw data can be after the text header, or in a separate file. If
        % there's a pointer to an external file, we assume that the data is
        % there
        if (isempty(rawfile) || strcmp(rawfile, 'LOCAL')) % data after text header
            % move file pointer to the beginning of the raw data
            if (fseek(fid, eoh, 'bof') == -1)
                error('Cannot read file');
            end
        else % data in external file
            % close mha file
            fclose(fid);
            
            % open raw file to read
            fid=fopen([pathstr filesep rawfile], 'r');
            if (fid<=0)
                error(['Cannot open file: ' pathstr filesep rawfile])
            end
        end
        
        % read all the raw data into a vector, because we cannot read it
        % into a 3D volume
        nrrd.data = fread(fid, prod(sz), [data_type '=>' data_type]);
        
        % reshape the data to create the data volume
        nrrd.data = reshape(nrrd.data, sz);
        
        % permute the X and Y coordinates
        nrrd.data = permute(nrrd.data, [2 1 3]);
        
        % close file
        fclose(fid);
        
        % check that we have enough data to create the output struct
        if (isempty(sz) || isempty(res) || isempty(offset))
            error('Incomplete header in .mha file')
        end
        
        % create output struct (we have read sz, res, etc in x, y, z order)
        for I = 1:N
            nrrd.axis(I).size = sz(I);
            nrrd.axis(I).spacing = res(I);
            nrrd.axis(I).min = offset(I) - res(I)/2;
            nrrd.axis(I).max = nrrd.axis(I).min + (sz(I)-1)*res(I);
            nrrd.axis(I).center = 1;
        end
        nrrd.axis(1).label = 'axis 1';
        nrrd.axis(2).label = 'axis 2';
        nrrd.axis(3).label = 'axis 3';
        nrrd.axis(1).unit = 'no unit';
        nrrd.axis(2).unit = 'no unit';
        nrrd.axis(3).unit = 'no unit';
        nrrd.axis = nrrd.axis';
        nrrd.property = [];
        
        % now we need to permute the axis so that we have [row, col,
        % slice]-> [y, x, z]
        nrrd.axis([1 2]) = nrrd.axis([2 1]);
        
    case '.lsm' % Carl Zeiss LSM format
        
        % read TIFF file
        warning('off', 'tiffread2:LookUp')
        stack = tiffread(file);
        warning('on', 'tiffread2:LookUp')
        
        % convert to sci format
        nrrd = scinrrd_tiff2nrrd(stack);
        
    otherwise
        
        error('Invalid file extension')
end

end

function s = getlabel(s, idx)
s = lower(strtrim(s(1:idx-1)));
end
function n = getnumval(s, idx)
n = str2num(s(idx+1:end));
end
