function lsmwrite(file, stack)
% LSMWRITE  Write a 3D stack to Zeiss LSM 5 format
%

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.0.0
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
error(nargchk(2, 2, nargin, 'struct'));
error(nargoutchk(0, 0, nargout, 'struct'));

% constants
lsminfolen = 500; % length in bytes of LSM info header

% force output filename to have an '.lsm' extension
[pathstr, file, ext] = fileparts(file);
if ~strcmp(ext, 'lsm')
    ext = 'lsm';
end
file = [pathstr filesep file '.' ext];

% open file to write. We use little endian format
fid = fopen(file, 'w', 'ieee-le');
if (fid < 0)
    error('Cannot open file for writing')
end

% byte order: II = little endian, MM = big endian
count = fwrite(fid, 'II', '*char');
if (count ~= 2)
    error('When writing byte order')
end

% number which identifies file as TIFF format
count = fwrite(fid, 42, 'uint16');
if (count ~= 1)
    error('When writing TIFF file identifier')
end

% byte offset for the first image file directory (IFD)
count = fwrite(fid, 8, 'uint32');
if (count ~= 1)
    error('When writing offset for the first image file directory (IFD)')
end

% offset of the LSM header
%
% the LSM header comes at the end of all the IFDs. The values that
% contribute to the offset value below are:
%
%   8: TIFF file header
%   12 * 11 * 1: 12 bytes each entry, 11 entries in the first frame
%   12 * 10 * (n-1): 12 bytes each entry, 10 entries in the rest of
%                    the frames
%   6 * n: 6 bytes for the header of each IFD
offsetlsm = 8 + 12 * 11 * 1 + 12 * 10 * (length(stack)-1) ...
    + 6 * length(stack);

% offsets for first data strip
offsetdatastrip = zeros(1, length(stack));
offsetdatastrip(1) = offsetlsm + lsminfolen;

% loop every slice
for I = 1:length(stack)
    
    % the first IFD has 10+1 entries: the 10 common to all, plus one for
    % the LSM info. The rest have only 10
    if (I == 1)
        count = fwrite(fid, 10+1, 'uint16');
    else
        count = fwrite(fid, 10, 'uint16');
    end
    if (count ~= 1)
        error(['When writing number of entries in IFD ' num2str(I) '/' num2str(length(stack))])
    end
    
    % tag 254: Subfile type. As we are ignoring thumbnails, value is
    % always 0
    writeifdentry(fid, 254, 0, 'uint32');
    
    % tag 256: ImageWidth (number of columns)
    writeifdentry(fid, 256, stack(I).width, 'uint32');
    
    % tag 257: ImageLength (number of rows)
    writeifdentry(fid, 257, stack(I).height, 'uint32');
    
    % tag 258: BitsPerSample
    writeifdentry(fid, 258, 8, 'uint16');
    
    % tag 259: Compression (1 = no compresion)
    writeifdentry(fid, 259, 1, 'uint16');
    
    % tag 262: PhotoMetricInterpretation (1 = black is zero, 2^N - 1 is white)
    writeifdentry(fid, 262, 1, 'uint16');
    
    % tag 273: StripOffsets
    writeifdentry(fid, 273, offsetdatastrip(I), 'uint32');
    [~, sz] = mat2tifftype(class(stack(I).data)); % size of data type in bytes
    striplen = sz * numel(stack(I).data); % length of data strip in bytes
    if (I < length(stack))
        offsetdatastrip(I+1) = offsetdatastrip(I) + striplen; % offset of next data strip
    end
    
    % tag 277: SamplesPerPixel (1 = 1 channel data; 3 = RGB)
    writeifdentry(fid, 277, 1, 'uint16');
    
    % tag 279: StripByteCounts (number of bytes in the slice)
    writeifdentry(fid, 279, striplen, 'uint32');
    
    % tag 284: PlanarConfiguration (1 = chunky format, component values for
    % each pixel are stored contiguously, e.g. RGBRGBRGB.... Given that we
    % are writing 1 channel data, this doesn't really matter)
    writeifdentry(fid, 284, 1, 'uint16');
    
    % the first frame also contains the LSM info
    if (I == 1)
       
        % tag 34412: Private tag for the LSM info header
        writeifdlsmentry(fid, offsetlsm, lsminfolen);
        
    end
    
    % last 4 bytes of the IFD, with the offset of the next IFD (0 if
    % there's no next IFD)
    if (I < length(stack))
        count = fwrite(fid, ftell(fid)+4, 'uint32');
    else
        count = fwrite(fid, 0, 'uint32');
    end
    if (count ~= 1)
        error(['When writing last 4 bytes of IFD ' num2str(I) '/' num2str(length(stack)) ' with offset to next IFD'])
    end
    
end % loop every slice

disp(['LSM header pos: ' num2str(ftell(fid))])%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% after writing one FID per frame, we should be at the beginning of the LSM
% info header
if (offsetlsm ~= ftell(fid))
    error('We are at the wrong position to start writing the LSM info header')
end

% write LSM info header
writelsminfo(fid, stack, lsminfolen);

disp(['Image data pos: ' num2str(ftell(fid))])%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
offsetdatastrip%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for I = 1:length(stack)

    if (offsetdatastrip(I) ~= ftell(fid))
        error(['We are at the wrong position to start writing the image data for frame ' num2str(I)])
    end

    % transpose the data matrix, so that voxels are written in the right
    % order
    stack(I).data = stack(I).data';
    
    % write the image data
    count = fwrite(fid, stack(I).data(:), class(stack(I).data));
    if (count ~= numel(stack(I).data)), error(['Image frame ' num2str(I)]), end

end

% close file
if (fclose(fid) ~= 0)
    error('Cannot close file')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% map Matlab to TIFF types
%
%  n: size in bytes
%  tifftype: code for TIFF type
function [tifftype, n] = mat2tifftype(matlabtype)

switch (matlabtype)
    case 'uint8'
        n = 1;
        tifftype = 1;
    case 'uchar'
        n = 1;
        tifftype = 2;
    case 'uint16'
        n = 2;
        tifftype = 3;
    case 'uint32'
        n = 4;
        tifftype = 4;
    case 'float32'
        n = 4;
        tifftype = 11;
    case 'float64'
        n = 8;
        tifftype = 12;
    otherwise
        error('matlab type %d not supported', matlabtype)
end

end

% write a single valued entry to the IFD
function writeifdentry(fid, tag, val, valtype)

    disp(['File pos: ' num2str(ftell(fid))])%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % tag ID
    count = fwrite(fid, tag, 'uint16');
    if (count ~= 1), error(['TIFF tag ' num2str(tag) ' (Tag ID)']), end
    
    % TIFF type of the value
    count = fwrite(fid, mat2tifftype(valtype), 'uint16');
    if (count ~= 1), error(['TIFF tag ' num2str(tag) ' (value type)']), end

    % number of values associated to this tag
    %
    % for the 10 common tags, all tags have 1 value. This changes with the
    % LSM info tag, that has more values (it's a struct), but we will deal
    % with it separately in another auxiliary function
    count = fwrite(fid, 1, 'uint32');
    if (count ~= 1), error(['TIFF tag ' num2str(tag) ' (number of values)']), end

    % tag value
    count = fwrite(fid, val, 'uint32');
    if (count ~= 1), error(['TIFF tag ' num2str(tag) ' (value)']), end

end

% write IFD entry in IFD that points to the LSM info somewhere else in the
% file
function writeifdlsmentry(fid, offsetlsm, lsminfolen)

    disp(['File pos: ' num2str(ftell(fid))])%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % constants
    tag = 34412;
    valtype = 'uint8';
    
    % tag ID
    count = fwrite(fid, tag, 'uint16');
    if (count ~= 1), error(['TIFF tag ' num2str(tag) ' (Tag ID)']), end
    
    % TIFF type of the value
    count = fwrite(fid, mat2tifftype(valtype), 'uint16');
    if (count ~= 1), error(['TIFF tag ' num2str(tag) ' (value type)']), end

    % number of values associated to this tag. lsminfolen is the number of
    % bytes we are going to write, although the LSM head provided by the
    % microscope has 500 bytes
    count = fwrite(fid, lsminfolen, 'uint32');
    if (count ~= 1), error(['TIFF tag ' num2str(tag) ' (number of values)']), end

    % offset value that points to the beginning of the LSM info block
    count = fwrite(fid, offsetlsm, 'uint32');
    if (count ~= 1), error(['TIFF tag ' num2str(tag) ' (offset)']), end
    
end

% write LSM info
function writelsminfo(fid, stack, lsminfolen)

    disp(['File pos to write LSM info: ' num2str(ftell(fid))])%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % MagicNumber
    count = fwrite(fid, hex2dec(stack(1).lsm.MagicNumber(3:end)), 'uint32');
    if (count ~= 1), error('LSM info header: MagicNumber'), end
    
    % StructureSize (size of the LSM info header in bytes)
    count = fwrite(fid, lsminfolen, 'uint32');
    if (count ~= 1), error('LSM info header: StructureSize'), end

    % DimensionX
    count = fwrite(fid, stack(1).width, 'uint32');
    if (count ~= 1), error('LSM info header: DimensionX'), end

    % DimensionY
    count = fwrite(fid, stack(1).height, 'uint32');
    if (count ~= 1), error('LSM info header: DimensionY'), end

    % DimensionZ
    count = fwrite(fid, length(stack), 'uint32');
    if (count ~= 1), error('LSM info header: DimensionZ'), end

    % DimensionChannels
    count = fwrite(fid, stack(1).lsm.DimensionChannels, 'uint32');
    if (count ~= 1), error('LSM info header: DimensionChannels'), end

    % DimensionTime
    count = fwrite(fid, stack(1).lsm.DimensionTime, 'uint32');
    if (count ~= 1), error('LSM info header: DimensionTime'), end

    % IntensityDataType
    count = fwrite(fid, stack(1).lsm.IntensityDataType, 'uint32');
    if (count ~= 1), error('LSM info header: IntensityDataType'), end

    % ThumbnailX (number of voxels)
    count = fwrite(fid, 0, 'uint32');
    if (count ~= 1), error('LSM info header: ThumbnailX'), end

    % ThumbnailY (number of voxels)
    count = fwrite(fid, 0, 'uint32');
    if (count ~= 1), error('LSM info header: ThumbnailY'), end

    % VoxelSizeX
    count = fwrite(fid, stack(1).lsm.VoxelSizeX, 'float64');
    if (count ~= 1), error('LSM info header: VoxelSizeX'), end

    % VoxelSizeY
    count = fwrite(fid, stack(1).lsm.VoxelSizeY, 'float64');
    if (count ~= 1), error('LSM info header: VoxelSizeY'), end

    % VoxelSizeZ
    count = fwrite(fid, stack(1).lsm.VoxelSizeZ, 'float64');
    if (count ~= 1), error('LSM info header: VoxelSizeZ'), end

    % OriginX
    count = fwrite(fid, stack(1).lsm.OriginX, 'float64');
    if (count ~= 1), error('LSM info header: OriginX'), end

    % OriginY
    count = fwrite(fid, stack(1).lsm.OriginY, 'float64');
    if (count ~= 1), error('LSM info header: OriginY'), end

    % OriginZ
    count = fwrite(fid, stack(1).lsm.OriginZ, 'float64');
    if (count ~= 1), error('LSM info header: OriginZ'), end

    % ScanType
    count = fwrite(fid, stack(1).lsm.ScanType, 'uint16');
    if (count ~= 1), error('LSM info header: ScanType'), end

    % SpectralScan
    count = fwrite(fid, stack(1).lsm.SpectralScan, 'uint16');
    if (count ~= 1), error('LSM info header: SpectralScan'), end

    % DataType
    count = fwrite(fid, stack(1).lsm.DataType, 'uint32');
    if (count ~= 1), error('LSM info header: DataType'), end

    % DataType
    count = fwrite(fid, stack(1).lsm.DataType, 'uint32');
    if (count ~= 1), error('LSM info header: DataType'), end
    
    % write a blank block to fill up the rest of the header
    count = fwrite(fid, zeros(1, lsminfolen - 100), 'uint8');
    if (count ~= lsminfolen - 100), error('LSM info header: Blank block'), end

%     % OffsetVectorOverlay
%     count = fwrite(fid, stack(1).lsm.OffsetVectorOverlay, 'uint32');
%     if (count ~= 1), error('LSM info header: OffsetVectorOverlay'), end
% 
%     % OffsetInputLut
%     count = fwrite(fid, stack(1).lsm.OffsetInputLut, 'uint32');
%     if (count ~= 1), error('LSM info header: OffsetInputLut'), end
% 
%     % OffsetOutputLut
%     count = fwrite(fid, stack(1).lsm.OffsetOutputLut, 'uint32');
%     if (count ~= 1), error('LSM info header: OffsetOutputLut'), end
% 
%     % OffsetChannelColors
%     count = fwrite(fid, stack(1).lsm.OffsetChannelColors, 'uint32');
%     if (count ~= 1), error('LSM info header: OffsetChannelColors'), end
% 
%     % TimeInterval
%     count = fwrite(fid, stack(1).lsm.TimeInterval, 'float64');
%     if (count ~= 1), error('LSM info header: TimeInterval'), end
% 
%     % OffsetChannelDataTypes
%     count = fwrite(fid, stack(1).lsm.OffsetChannelDataTypes, 'uint32');
%     if (count ~= 1), error('LSM info header: OffsetChannelDataTypes'), end
% 
%     % OffsetScanInformation
%     count = fwrite(fid, stack(1).lsm.OffsetScanInformation, 'uint32');
%     if (count ~= 1), error('LSM info header: OffsetScanInformation'), end
% 
%     % OffsetKsData
%     count = fwrite(fid, stack(1).lsm.OffsetKsData, 'uint32');
%     if (count ~= 1), error('LSM info header: OffsetKsData'), end
%    
%     % OffsetTimeStamps
%     count = fwrite(fid, stack(1).lsm.OffsetTimeStamps, 'uint32');
%     if (count ~= 1), error('LSM info header: OffsetTimeStamps'), end
%    
%     % OffsetEventList
%     count = fwrite(fid, stack(1).lsm.OffsetEventList, 'uint32');
%     if (count ~= 1), error('LSM info header: OffsetEventList'), end
%    
%     % OffsetRoi
%     count = fwrite(fid, stack(1).lsm.OffsetRoi, 'uint32');
%     if (count ~= 1), error('LSM info header: OffsetRoi'), end
%    
%     % OffsetBleachRoi
%     count = fwrite(fid, stack(1).lsm.OffsetBleachRoi, 'uint32');
%     if (count ~= 1), error('LSM info header: OffsetBleachRoi'), end
%    
%     % OffsetNextRecording
%     count = fwrite(fid, stack(1).lsm.OffsetNextRecording, 'uint32');
%     if (count ~= 1), error('LSM info header: OffsetNextRecording'), end

end
