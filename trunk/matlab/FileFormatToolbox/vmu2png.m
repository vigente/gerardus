function vmu2png(file)
% VMU2PNG  Convert a microscope image file from .vmu to .png format.
%
% vmu2png(FILE)
%
%   FILE is a string with the path and name of a .vmu (Uncompressed Virtual
%   Microscope Specimen) Hamamatsu file. It is expected that the .vmu file
%   contains only the header information, and a tag 'ImageFile(0)' with the
%   name of the .raw file with the actual image data.
%
%   This function reads the image and all the headers from the .vmu file
%   and produces a .png file with the same name. PNG files use lossless
%   Deflate compression.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
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


% check arguments
narginchk(1, 1);
nargoutchk(0, 0);

%% Read input .vmu file

% extract extension of filename in lower case
[pathstr, filestr, ext] = fileparts(file);
ext = lower(ext);

if (~strcmpi(ext, '.vmu'))
    error('Input file must have extension .vmu, .VMU')
end

% init size vector
sz = [0 0];

% initialize cell array to put all the metadata tags that we are going to
% pass to the PNG file
metadata = {};

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
        tline = fgetl(fid);
        continue
    end
    
    % split line into label and value
    label = tline(1:idx-1);
    value = tline(idx+1:end);
    
    switch label
        
        case 'NoLayers'
            
            check_tag_value(label, value, '1');
            metadata{end+1} = label;
            metadata{end+1} = value;
            
        case 'LayerSpacing'
        
            check_tag_value(label, value, '0');
            metadata{end+1} = label;
            metadata{end+1} = value;
            
        case 'LayerOffset'
            
            check_tag_value(label, value, '0');
            metadata{end+1} = label;
            metadata{end+1} = value;
            
        case 'ImageFile(0)'
            
            % name of file with the image data
            rawfile = value;
            
        case 'MapFile'
            
            % this is a file name of a file with a lower-resolution version
            % of all the ImageFiles. We ignore it.
            
        case 'MapScale'
            
            % downsample factor of the MapFile. We ignore it.
            
        case 'PixelWidth'
            
            sz(2) = str2double(value);
%             metadata{end+1} = 'Width';
%             metadata{end+1} = sz(2);
            
        case 'PixelHeight'
            
            sz(1) = str2double(value);
%             metadata{end+1} = 'Height';
%             metadata{end+1} = sz(1);
            
        case 'PhysicalWidth'
            
            % units are nm in the .vmu file
            physicalWidth = str2double(value) * 1e-9;
            
            % resolution in the PNG sense, i.e. pixels / meter
            xResolution = sz(2) / physicalWidth;

            % PNG valid field
            metadata{end+1} = 'XResolution';
            metadata{end+1} = xResolution;
            
        case 'PhysicalHeight'
            
            % units are nm in the .vmu file
            physicalHeight = str2double(value) * 1e-9;
            
            % resolution in the PNG sense, i.e. pixels / meter
            yResolution = sz(1) / physicalHeight;
            
            metadata{end+1} = 'YResolution';
            metadata{end+1} = yResolution;
            
        case 'XOffsetFromSlideCentre'
            
            % units are nm in the .vmu file
            xOffset = str2double(value) * 1e-9 - physicalWidth / 2;

            metadata{end+1} = 'XOffset';
            metadata{end+1} = sprintf('%0.13e', xOffset);
            
            
        case 'YOffsetFromSlideCentre'
            
            % units are nm in the .vmu file
            yOffset = str2double(value) * 1e-9 - physicalHeight / 2;

            metadata{end+1} = 'YOffset';
            metadata{end+1} = sprintf('%0.13e', yOffset);
            
        case 'BitsPerPixel'
            
            % we assume that the input image is RGB. Thus, if
            % BitsPerPixel=24 because of the 3 channels, then the BitDepth
            % in PNG is 24/3 = 8
            numChannels = 3;
            bitDepth = str2double(value) / numChannels;
            metadata{end+1} = 'BitDepth';
            metadata{end+1} = bitDepth;
            numBytesPerPixelPerChannel = bitDepth / 8;
            
        case 'PixelOrder'
            
            check_tag_value(label, value, 'RGB');
            
        case 'Creator'
            
            % software that created the .vmu file, e.g. NDP.scan 2.5.86
            metadata{end+1} = 'Software';
            metadata{end+1} = value;
            
        case 'IlluminationMode'
            
            check_tag_value(label, value, '0');
            
        case 'ExposureMultiplier'
            
            check_tag_value(label, value, '1.000000');
            
        case {'GainRed', 'GainGreen', 'GainBlue'}
            
            check_tag_value(label, value, '1');
            
        case 'NoSkippedLanes'
            
            check_tag_value(label, value, '0');
            
        case {'BlobMapWidth', 'BlobMapHeight', 'BlobMapImage'}
            
            % e.g. R-001 - 2014-03-25 14.39.17_blobmap.raw
            % We ignore it, as we are not going to load the blob data into
            % the PNG file

        otherwise % tag not supported by the PNG format
            
            % 'SourceLens':     magnification factor of the lens
            % 'Reference':      first part of the filename, without the
            %                   date, etc.
            % 'AuthCode':       unknown
            % 'ExposureTime':   unknown
            % 'FocalPlaneTolerance': unknown
            % 'NMP':            unknown
            % 'MacroIllumination': unknown
            % 'MacroLabelObscured': unknown
            % 'FocusOffset':    unknown
            % 'RefocusInterval': unknown
            % 'CubeName':       unknown
            % 'HardwareModel':  microscope model
            % 'HardwareSerial': microscope serial number
            % 'FirmwareVersion': unknown
            % 'CalibrationInfoFile': unknown
            % 'NoFocusPoints':  number of focus points
            % 'FocusPoint0X', 'FocusPoint1X', ... : coordinates of focus points
            % 'FocusPoint0Y', 'FocusPoint1Y', ... : coordinates of focus points
            % 'FocusPoint0Z', 'FocusPoint1Z', ... : coordinates of focus points
            % 'NoBlobPoints':   number of blob points
            % 'BlobPoint0Blob', 'BlobPoint1Blob', ...: unknown
            % 'BlobPoint0FocusPoint', 'BlobPoint1FocusPoint', ...: unknown
            % 'Wavelength':     unknown
            % 'LampAge':        unknown
            % 'FocusTime':      unknown
            % 'ScanTime':       unknown
            % 'WriteTime':      unknown
            
            % we will pass the tag to the PNG file anyway, and it will get
            % stored as an "OtherText" tag
            metadata{end+1} = label;
            metadata{end+1} = value;
            
    end
    
    % read next line
    tline = fgetl(fid);
    
end
fclose(fid);

% translate number of bits per pixel per channel to Matlab data type
switch numBytesPerPixelPerChannel
    case 1
        pixeltype = 'uint8';
    otherwise
        error('Unknown data type. File does not contain one byte per channel per pixel')
end


% read image data to an array with the appropriate pixel type
% (note: '*uint8' is shorthand for 'uint8=>uint8')
fid = fopen([pathstr filesep rawfile], 'r');
if (fid == -1)
    error('Gerardus:CannotReadFile', ...
        ['Cannot open file for reading:\n' file])
end
im = fread(fid, numChannels * sz(1) * sz(2), ['*' pixeltype]);
fclose(fid);

% reshape image array to produce an image R*C*channels that can be
% visualized with imagesc()
im = reshape(im, [numChannels sz(2) sz(1)]);
im = permute(im, [3 2 1]);

%% Wrte output .png file

% output filename
outfile = [pathstr filesep filestr '.png'];

% write output PNG file
imwrite(im, outfile, 'PNG', metadata{:});

end

% auxiliary function to check expected value in a .vmu tag, and give a
% warning if it's different
function check_tag_value(tag, value, expval)
if (~strcmp(value, expval))
    warning(['Input image has ' tag '=' value '. We expected =' expval '.'])
end
end
