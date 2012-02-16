function nrrd = scinrrd_tiff2nrrd(stack)
% SCINRRD_TIFF2NRRD  Create SCI NRRD struct from TIFF stack
%
% This function creates a struct with the correct format that the Gerardus
% Toolbox uses for NRRD variables. This is the same format you obtain when
% loading a .mat file using scinrrd_load(), and can be saved to a .mat file
% using scinrrd_save().
%
% NRRD = SCINRRD_TIFF2NRRD(STACK)
%
%   STACK is a struct array obtained from loading a TIFF or LSM file with
%   tiffread(). See below for details on both formats.
%
%   Note: We assume that all frames in stack have the same resolution and
%   offset (only values from frame 1 are read).
%
%   When the microscope grayscale data is exported as RGB, only the red
%   channel is used. So if the image data is RGB, we only read the red
%   channel.
%
% =========================================================================
%   TIFF format:
%
%    >> stack = tiffread('file.tif');
%    >> stack = 
% 
%    1x185 struct array with fields:
%        filename
%        width
%        height
%        bits
%        info
%        x_resolution
%        y_resolution
%        resolution_unit
%        cmap
%        colors
%        data
% 
%    >> stack(1)
%
%             filename: 'file.tif'
%                width: 512
%               height: 512
%                 bits: 8
%                 info: [1x82 char]
%         x_resolution: [2x1 double]
%         y_resolution: [2x1 double]
%      resolution_unit: 1
%                 cmap: [768x1 double]
%               colors: 256
%                 data: [512x512 uint8]
%
%    >> stack(1).info
% 
%    ans =
% 
%    ImageJ=1.44i
%    images=185
%    slices=185
%    unit=um
%    spacing=0.6560000000000001
%    loop=false
%
% =========================================================================
%   LSM format:
%
%    >> stack = tiffread('file.lsm');
%    >> stack = 
%
%            filename: 'file.lsm'
%               width: 512
%              height: 512
%                 bits: 8
%                 data: [512x512 uint8]
%                  lsm: [1x1 struct]
%
%    >> stack.lsm
%
%    ans = 
%
%              MagicNumber: '0x00400494C'
%               DimensionX: 512
%               DimensionY: 512
%               DimensionZ: 1
%        DimensionChannels: 1
%            DimensionTime: 1
%        IntensityDataType: 1
%               ThumbnailX: 128
%               ThumbnailY: 128
%               VoxelSizeX: 4.3945e-07
%               VoxelSizeY: 4.3945e-07
%               VoxelSizeZ: 6.5600e-07
%                  OriginX: 0
%                  OriginY: 0
%                  OriginZ: 0
%                 ScanType: 0
%             SpectralScan: 0
%                 DataType: 0
%             TimeInterval: 0
%                TimeStamp: 4.0385e+03
%               TimeOffset: 0
         
% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.2.4
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

% if stack was read from a TIFF file, it will have an 'info' field
if isfield(stack, 'info')
    
    % image volume
    nrrd.data = cat(3, stack.data);

    % image resolution in x and y
    nrrd.axis(1).spacing = 1/stack(1).x_resolution(1);
    nrrd.axis(2).spacing = 1/stack(1).y_resolution(1);
    
    % image resolution in z
    a = strfind(stack(1).info, 'spacing=');
    b = strfind(stack(1).info, 'loop=');
    if (isempty(a) || isempty(b) || b < a)
        error('TIFF stack has no Z-spacing info')
    end
    nrrd.axis(3).spacing = str2double(stack(1).info(a+8:b-2));
    
    b = a;
    a = strfind(stack(1).info, 'unit=');
    unit = stack(1).info(a+5:b-2);
    switch unit
        case 'm'
        case 'dm'
            nrrd.axis(3).spacing = nrrd.axis(3).spacing * 1e-1;
        case 'cm'
            nrrd.axis(3).spacing = nrrd.axis(3).spacing * 1e-2;
        case 'mm'
            nrrd.axis(3).spacing = nrrd.axis(3).spacing * 1e-3;
        case 'um'
            nrrd.axis(3).spacing = nrrd.axis(3).spacing * 1e-6;
        otherwise
            error('Z-axis units not recognised')
    end

    % loop some of the fields
    for I = 1:3
        
        % left edge of first voxel
        nrrd.axis(I).min = -nrrd.axis(I).spacing / 2;
        
        % left edge of last voxel
        nrrd.axis(I).max = nrrd.axis(1).min ...
            + (size(nrrd.data, I) - 1) * nrrd.axis(I).spacing;
        
    end
    
elseif isfield(stack, 'lsm') % stack read from an LSM v5 file

    % image volume
    if (stack(1).lsm.DimensionChannels == 1)
        
        % grayscale data
        nrrd.data = cat(3, stack.data);
        
    elseif (stack(1).lsm.DimensionChannels == 3)
        
        % RGB data
        nrrd.data = zeros([size(stack(1).data{1}) length(stack)], ...
            class(stack(1).data{1}));
        for I = 1:length(stack)
            nrrd.data(:, :, I) = stack(I).data{1};
        end
        
    else
        error('Unrecognized number of channels in LSM file')
    end
    
    % the microscope camera can save to 12-bit, but this is read as uint16
    % by the LSM reading function
    if ((stack(1).lsm.IntensityDataType == 2) ...
            && strcmp(class(nrrd.data), 'uint16')) % LSM data is 12-bit
        
        if (max(nrrd.data(:)) > (2^12-1))
            warning('I think that LSM data is 12-bit, but values are larger than expected')
        end
        
    end
    
    % voxel resolution
    nrrd.axis(1).spacing = stack(1).lsm.VoxelSizeX;
    nrrd.axis(2).spacing = stack(1).lsm.VoxelSizeY;
    nrrd.axis(3).spacing = stack(1).lsm.VoxelSizeZ;
    
    % left edge of first voxel
    nrrd.axis(1).min = stack(1).lsm.OriginX - nrrd.axis(1).spacing / 2;
    nrrd.axis(2).min = stack(1).lsm.OriginY - nrrd.axis(2).spacing / 2;
    nrrd.axis(3).min = stack(1).lsm.OriginZ - nrrd.axis(3).spacing / 2;
    
    for I = 1:3
        
        % left edge of last voxel
        nrrd.axis(I).max = nrrd.axis(1).min ...
            + (size(nrrd.data, I) - 1) * nrrd.axis(I).spacing;
        
    end
    
else
    
    error('Stack struct format not recognised')
    
end

% loop some of the fields
for I = 1:3
    
    % data volume size
    nrrd.axis(I).size = size(nrrd.data, I);
    
    % unused
    nrrd.axis(I).center = 1;
    nrrd.axis(I).unit = 'no unit';
    
end

% other
nrrd.axis(1).label = 'axis 2';
nrrd.axis(2).label = 'axis 1';
nrrd.axis(3).label = 'axis 3';

% we need nrrd.axis to be a column vector
nrrd.axis = nrrd.axis';
