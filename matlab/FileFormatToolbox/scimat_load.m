function scimat = scimat_load(file, varargin)
% SCIMAT_LOAD  Load an image into a SCIMAT struct from a Matlab, MetaImage,
% Carl Zeiss LSM, Hamamatsu VMU or regular graphics file (PNG, TIFF, etc).
%
% SCIMAT = SCIMAT_LOAD(FILE)
%
%   This function loads the image and metainformation into a scimat struct.
%   If necessary, it swaps rows and columns to follow Matlab's convention
%   that (rows, columns) <=> (y, x).
%
%   FILE is a string with the path and name of the file that contains the
%   2D or 3D image:
%
%     .mat: Matlab binary file with a "scirunnrrd" struct (see below for
%           details).
%
%     .mha/.mhd: MetaImage file (developed for the ITK and VTK libraries).
%          The .mha/.mhd file can be pure text, and point to .raw file with
%          the image, or contain both the header and the binary data.
%          Compressed and uncompressed images are supported
%          (http://www.itk.org/Wiki/ITK/MetaIO/Documentation).
%
%     .lsm: Carl Zeiss microscopy image format. Binary file.
%
%     .vmu: Hamamatsu Uncompressed Virtual Microscope Specimen. Text file
%           containing only the image metadata, and a path to the file with
%           the actual binary image data.
%
%     Otherwise the file is assumed to be a regular graphics file, e.g.
%     .tif, .png, .jpg, .bmp, etc. Any file that imread() can read is
%     valid.
%
%   SCIMAT is the struct with the image data and metainformation (see "help
%   scimat" for details).
%
% ... = SCIMAT_LOAD(FILE, parameter, value, ...)
%
%   The function can be followed by parameter/value pairs to modify the
%   behaviour of the function. Currently, the only one implemented is
%
%   'HeaderOnly': (def false) Read only the metainformation from the file,
%       stopping before reading the image itself. This returns
%       SCIMAT.data = [], and can save time when the intensities are not
%       required. (Currently only used for .mha/.mhd, .tif, .png files).
%
% See also: scimat, scimat_save.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010-2016 University of Oxford
% Version: 0.5.10
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

% parse input arguments
parser = inputParser;
addRequired(parser, 'file', @isstr);
addParameter(parser, 'HeaderOnly', false , @islogical);
parser.parse(file, varargin{:});
headerOnly = parser.Results.HeaderOnly;

% extract extension of filename in lower case
[pathstr, ~, ext] = fileparts(file);
ext = lower(ext);

switch lower(ext)
    
    case '.mat' % Matlab file in Seg3D scirunnrrd format
        
        % load data
        scimat = load(file);
        
        % rename SCIMAT volume for convenience
        scimat = scimat.scirunnrrd;
        
        % remove dummy dimension in old files created with Seg3D 1.x
        scimat = scimat_squeeze(scimat);
        
        % correct x-,y-coordinates
        scimat = scimat_seg3d2matlab(scimat);
        
        % remove extra metainformation that is not used
        if isfield(scimat, 'property')
            scimat = rmfield(scimat, 'property');
        end
        if isfield(scimat, 'axis')
            for f = {'max', 'center', 'label', 'unit'}
                if isfield(scimat.axis, f)
                    scimat.axis = rmfield(scimat.axis, f);
                end
            end
        end
        
        % empty rotation matrix
        scimat.rotmat = [];
        
    case {'.mha', '.mhd'} % MetaImage file
        
        % open file to read
        fid = fopen(file, 'r');
        if (fid <= 0)
            error(['Cannot open file: ' file])
        end
        
        % default values for variables that we are trying to read from the
        % text header
        N = [];
        sz = [];
        data_type = [];
        offset = [];
        res = [];
        nchannel = 1;
        rotmat = [];
        rawfile = [];
        dataIsCompressed = [];
        compressedSize = [];
        
        % process text header, and stop if we get to the raw data
        while (true)
            
            % read text header line
            tline = fgetl(fid);
            
            % if end of text header stop reading
            % Warning: this only works if the first voxel=0. Otherwise, we
            % would read some voxels as part of a header line. To avoid it,
            % we now break after finding ElementDataFile
            if (tline(1) == 0), break, end
            
            % pointer to current position in header
            eoh = ftell(fid);
            
            % parse text header line
            
            % find location of "=" sign
            idx = strfind(tline, '=');
            if isempty(idx), break, end
            
            switch getlabel(tline, idx)
                case 'ObjectType'
                    if (~strcmpi(strtrim(tline(idx+1:end)), 'Image'))
                        error('ObjectType ~= Image. Functionality not implemented')
                    end
                case 'NDims'
                    N = getnumval(tline, idx);
                case 'DimSize'
                    sz = getnumval(tline, idx);
                case 'ElementType'
                    switch lower(strtrim(tline(idx+1:end)))
                        case 'met_uchar'
                            data_type = 'uint8';
                            data_size = 1;
                        case 'met_char'
                            data_type = 'int8';
                            data_size = 1;
                        case 'met_ushort'
                            data_type = 'uint16';
                            data_size = 2;
                        case 'met_short'
                            data_type = 'int16';
                            data_size = 2;
                        case 'met_uint'
                            data_type = 'uint32';
                            data_size = 4;
                        case 'met_int'
                            data_type = 'int32';
                            data_size = 4;
                        case 'met_float'
                            data_type = 'single';
                            data_size = 5;
                        case 'met_double'
                            data_type = 'double';
                            data_size = 6;
                        otherwise
                            error('Unrecognized ElementType')
                    end
                case {'Offset', 'Position', 'Origin'}
                    offset = getnumval(tline, idx);
                case 'ElementSpacing'
                    res = getnumval(tline, idx);
                    if (all(res == 0))
                        warning('File provides ElementSpacing(:) == 0. This is a physical impossibility. Changing to ElementSpacing(:) = 1')
                        res(:) = 1;
                    end
                case 'ElementByteOrderMSB'
                    if (getnumval(tline, idx) ~= true)
                        error('ElementByteOrderMSB = False. Functionality not implemented')
                    end
                case 'BinaryDataByteOrderMSB'
                    if (getnumval(tline, idx) ~= true)
                        error('BinaryDataByteOrderMSB = False. Functionality not implemented')
                    end
                case 'ElementNumberOfChannels'
                    nchannel = getnumval(tline, idx);
                case 'ElementDataFile'
                    rawfile = strtrim(tline(idx+1:end));
                    break;
                case 'CompressedData'
                    dataIsCompressed = strcmpi(strtrim(tline(idx+1:end)), 'true');
                case {'Orientation', 'TransformMatrix', 'Rotation'}
                    % orientation of X axis, orientation of Y axis, ...
                    rotmat = getnumval(tline, idx);
                    rotmat = reshape(rotmat, [N, N])';
                case 'CenterOfRotation'
                    rotc = getnumval(tline, idx);
                    if (any(rotc ~= 0))
                        error('Image has center of rotation different from 0. That option has not been implemented yet in this reader')
                    end
                case 'HeaderSize'
                    if (getnumval(tline, idx) ~= -1)
                        error('HeaderSize ~= -1. Functionality not implemented')
                    end
                case 'BinaryData'
                    if (getnumval(tline, idx) ~= true)
                        error('BinaryData = False. Functionality not implemented')
                    end
                case 'AnatomicalOrientation'
                    % ignore
                case 'CompressedDataSize'
                    compressedSize = getnumval(tline, idx);
                otherwise
                    warning(['Unrecognized line: ' tline])
            end
            
        end
        
        % the raw data can be after the text header, or in a separate file. If
        % there's a pointer to an external file, we assume that the data is
        % there
        if (isempty(rawfile))
            error('No pointer to data in header')
        elseif (strcmp(rawfile, 'LOCAL')) % data after text header
            % do nothing
        else % data in external file
            % close mha file
            fclose(fid);
            
            % open raw file to read
            fid = fopen([pathstr filesep rawfile], 'r');
            if (fid <= 0)
                error(['Cannot open file: ' pathstr filesep rawfile])
            end
        end
        
        if (headerOnly) % don't waste time reading the data
            
            scimat.data = [];
            
        else % read the image data
            
            % read data, and decompress if necessary
            if (dataIsCompressed)
                
                % reposition the reading pointer. In principle, it should be
                % after the last \n (ASCII 10, LF new line feed) after
                % "ElementDataFile = LOCAL". However, if the first byte of the
                % image data is ASCII 13 (CR = carriage return), then fgetl
                % will have interpreted LF+CR as the end of line, and the
                % pointer will be one byte too far ahead
                fseek(fid, -compressedSize, 'eof');
                
                % read all the raw data into a vector
                if (isempty(compressedSize))
                    
                    warning('File did not provide CompressedDataSize. Trying to read to the end of file')
                    scimat.data = fread(fid, prod(sz) * nchannel, ...
                        [data_type '=>' data_type]);
                    
                else
                    
                    scimat.data = fread(fid, compressedSize, ...
                        [data_type '=>' data_type]);
                    
                end

                % decompress the data
                scimat.data = zlib_decompress(scimat.data, data_type);
                
            else
                
                % reposition the reading pointer. In principle, it should be
                % after the last \n (ASCII 10, LF new line feed) after
                % "ElementDataFile = LOCAL". However, if the first byte of the
                % image data is ASCII 13 (CR = carriage return), then fgetl
                % will have interpreted LF+CR as the end of line, and the
                % pointer will be one byte too far ahead
                fseek(fid, -prod(sz) * nchannel * data_size, 'eof');
                
                % read all the raw data into a vector
                scimat.data = fread(fid, prod(sz) * nchannel, ...
                    [data_type '=>' data_type]);
                if (length(scimat.data) ~= prod(sz) * nchannel)
                    error(['We read ' num2str(length(scimat.data)) ...
                        ' byte from the file instead of ' ...
                        num2str(prod(sz) * nchannel) ' byte'])
                end

            end
            
            % reshape the data to create the data volume
            scimat.data = reshape(scimat.data, [nchannel sz]);
            
            % rearrange dimensions so that channels are the 5th dimension
            scimat.data = permute(scimat.data, [2 3 4 5 1]);
            
            % permute the X and Y coordinates
            scimat.data = permute(scimat.data, [2 1 3:ndims(scimat.data)]);
            
        end

        % close file
        fclose(fid);
        
        % defaults
        if (isempty(offset))
            offset = zeros(1, N);
        end
        if (isempty(res))
            res = ones(1, N);
        end
        if (isempty(sz))
            sz = size(scimat.data);
        end
        if (isempty(rotmat))
            rotmat = eye(N);
        end

        % create output struct (we have read sz, res, etc in x, y, z order)
        for I = 1:N
            scimat.axis(I).size = sz(I);
            scimat.axis(I).spacing = res(I);
            scimat.axis(I).min = offset(I) - res(I)/2;
        end
        scimat.axis = scimat.axis';
        
        % now we need to permute the axis so that we have [row, col,
        % slice]-> [y, x, z]
        scimat.axis([1 2]) = scimat.axis([2 1]);
        
        % rotation matrix
        scimat.rotmat = rotmat;
        
    case '.lsm' % Carl Zeiss LSM format
        
        % read TIFF file
        warning('off', 'tiffread2:LookUp')
        stack = tiffread(file);
        warning('on', 'tiffread2:LookUp')
        
        % convert to sci format
        scimat = scimat_tiff2scimat(stack);
        
    case '.vmu' % Hamamatsu miscroscope format (Uncompressed Virtual Microscope Specimen)
        
        % init resolution and size vectors
        rawfile = '';
        res = [0 0 NaN];
        sz = [0 0 1];
        offset = [0 0 0];
        bpp = [];
        pixelorder = '';
        
        % open header file to read
        fid = fopen(file, 'r');
        if (fid == -1)
            error('Gerardus:CannotReadFile', ...
                ['Cannot open file for reading:\n' file])
        end
        
        % read metainfo from header
        tline = fgetl(fid);
        while ischar(tline)
            
            % find the '=' in the line
            idx = strfind(tline, '=');
            if (isempty(idx))
                % if this line has no label=value format, skip to next
            end
            
            % split line into label and value
            label = tline(1:idx-1);
            value = tline(idx+1:end);
            
            switch label
                
                case 'ImageFile(0)'
                    
                    % name of file with the image data
                    rawfile = value;

                case 'PixelWidth'
                    
                    sz(2) = str2double(value);
                
                case 'PixelHeight'
                    
                    sz(1) = str2double(value);
                
                case 'PhysicalWidth'
                    
                    % units are nm
                    res(2) = str2double(value) * 1e-9;
                
                case 'PhysicalHeight'
                    
                    % units are nm
                    res(1) = str2double(value) * 1e-9;
                
                case 'XOffsetFromSlideCentre'
                    
                    offset(1) = str2double(value) * 1e-9;
                
                case 'YOffsetFromSlideCentre'
                    
                    offset(2) = str2double(value) * 1e-9;
                    
                case 'BitsPerPixel'
                    
                    bpp = str2double(value);
                    
                case 'PixelOrder'
                    
                    pixelorder = value;
                    
            end
            
            % read next line
            tline = fgetl(fid);
            
        end
        fclose(fid);
        
        % at this point in the code, res, sz -> (rows, cols, slices)
        
        % compute pixel size from the total image size
        res(1) = res(1) / sz(1);
        res(2) = res(2) / sz(2);
        
        % check that the image is RGB
        if (~strcmpi(pixelorder, 'RGB'))
            error('Microscopy image is not RGB')
        else
            numchannels = 3;
        end

        % check that we know the number of bits per pixel
        if (isempty(bpp))
            warning('File does not provide field BitsPerPixel. Assuming 3 bytes per pixel')
            bpp = 24;
        end
        
        % number of bytes per channel per pixel
        numbyte = bpp / 8 / numchannels;
        
        % translate number of bits per pixel to Matlab data type
        switch numbyte
            case 1
                pixeltype = 'uint8';
            otherwise
                error('Unknown data type. File does not contain one byte per channel per pixel')
        end
        
        
        % read image data to an array with the appropriate pixel type
        % (note: '*uint8' is shorthand for 'uint8=>uint8')
        [pathstr, ~, ~] = fileparts(file);
        fid = fopen([pathstr filesep rawfile], 'r');
        if (fid == -1)
            error('Gerardus:CannotReadFile', ...
                ['Cannot open file for reading:\n' file])
        end
        im = fread(fid, numchannels * sz(1) * sz(2), ['*' pixeltype]);
        fclose(fid);
        
        % reshape image array to produce an image R*C*channels that can be
        % visualized with imagesc()
        im = reshape(im, [numchannels sz(2) sz(1)]);
        im = permute(im, [3 2 1]);
        
        % add a dummy dimension, so that it is clear that we don't have 3
        % slices instead of 3 channels
        im = reshape(im, [size(im, 1) size(im, 2) 1 size(im, 3)]);
        
        % create SCI MAT struct
        scimat = scimat_im2scimat(im, res, offset);
        
    otherwise % we assume it's a normal image (.tif, .png, .bmp, etc)
        
        % read image metadata header
        info = imfinfo(file);
        
        % if the image has a thumbnail, we discard the thumbnail
        if (length(info) > 1)
            
            info = info(1);
            
        end
        
        if (headerOnly)
            
            scimat.data = [];
            
        else
        
            % read image pixels
            scimat.data = imread(file);
            
            % shift dimensions so that number of channels is at 5th
            % dimension, to comply with the scimat format
            scimat.data = reshape(scimat.data, ...
                [size(scimat.data, 1) size(scimat.data, 2) 1 1 ...
                size(scimat.data, 3)]);
            
        end

        % spacing units
        if (~isfield(info, 'ResolutionUnit') || isempty(info.ResolutionUnit))
            warning('File does not provide ResolutionUnit. Assuming ''meter''')
            info.ResolutionUnit = 'meter';
        end
        switch (info.ResolutionUnit)
            case 'Centimeter'
                
                unit = 1e-2;
                
            case 'meter'
                
                unit = 1;
                
            otherwise
                
                error(['I do not know how to deal with the ResolutionUnit provided by this file: ' info.ResolutionUnit])
                
        end
        
        % checking missing fields
        if (~isfield(info, 'XResolution') || isempty(info.XResolution))
            warning('File does not provide XResolution. Assuming ''1.0''')
            info.XResolution = 1.0;
        end
        if (~isfield(info, 'YResolution') || isempty(info.YResolution))
            warning('File does not provide YResolution. Assuming ''1.0''')
            info.YResolution = 1.0;
        end
        
        % set scimat axes values
        scimat.axis(1).size = info.Height;
        scimat.axis(2).size = info.Width;
        scimat.axis(1).spacing = unit ./ info.YResolution;
        scimat.axis(2).spacing = unit ./ info.XResolution;
        scimat.axis(1).min = -0.5 * scimat.axis(1).spacing;
        scimat.axis(2).min = -0.5 * scimat.axis(2).spacing;
        
        % we assume that the image is oriented parallel to Cartesian axes
        scimat.rotmat = eye(2);
        
end

end

%% Auxiliary functions

% MHA file headers have lines like e.g. 
%   BinaryData = True
% These two functions getlabel() and getnumval() extract the left and right
% sides of the "=", respectively. The latter also converts the string to
% numerical values. If you don't expect a numerical value, don't use
% getnumval(), and instead just get that part of the string.
function s = getlabel(s, idx)
s = strtrim(s(1:idx-1));
end
function n = getnumval(s, idx)
n = str2num(s(idx+1:end));
end

% Decompress compressed MHA images. Code copied from mha_read_volume() by
% Dirk-Jan Kroon, 10 Nov 2010 (Updated 23 Feb 2011)
% http://uk.mathworks.com/matlabcentral/fileexchange/29344-read-medical-data-3d/content/mha/mha_read_volume.m
function M = zlib_decompress(Z,DataType)

import com.mathworks.mlwidgets.io.InterruptibleStreamCopier
a=java.io.ByteArrayInputStream(Z);
b=java.util.zip.InflaterInputStream(a);
isc = InterruptibleStreamCopier.getInterruptibleStreamCopier;
c = java.io.ByteArrayOutputStream;
isc.copyStream(b,c);
M=typecast(c.toByteArray,DataType);

end
