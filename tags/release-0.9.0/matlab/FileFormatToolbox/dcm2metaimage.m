function varargout = dcm2metaimage( str, res, scale, crop, ext, file )
% DCM2METAIMAGE  Read a batch of DICOM files, collate them and save as a
% single MetaImage volume (one .mha and one .raw file)
%
% DCM2METAIMAGE(STR)
%
%   STR is a string with the path and output mha header file name. For
%   example,
%
%     STR='/home/john/data/study01/study01.mha'; % linux
%     STR='C:\data\study01\study01.mha';         % windows
%
%   It is assumed that the MetaImage file (.mha) will be created from the
%   DICOM files found in the target directory.
%
% DCM2METAIMAGE(STR, RES, SCALE, CROP, EXT, FILE)
%
%   RES is a 3-vector with the pixel size in the x-, y- and z-coordinates.
%   If RES is not provided, the function tries to read the values from the
%   DICOM header. If the resolution values cannot be read, by default
%   RES=[1 1 1].
%
%   SCALE is a factor to reduce the size of each frame (the number of
%   frames doesn't change). By default, SCALE=1.
%
%   CROP is a vector with values [xmin, xmax, ymin, ymax, zmin, zmax].
%   These values allow to crop the data volume. It follows the same
%   convention as Seg3D [1]. For example, if xmin=0, xmax=3, then voxels 1,
%   2 and 3 are selected.
%   Note that the output MetaImage file will have an offset header so that
%   the data retains its true coordinates.
%
%   EXT is a string with the file extension of the images. By default,
%   EXT='tif'.
%
%   FILE is the list of image files. By default, all the files in the
%   target directory are read, but FILE allows to change the system listing
%   order, select a subset of files, etc. This is useful, for example, when
%   the system listing order does not correspond to the order in which you
%   want to collate the slices. FILE is expected to have the same struct
%   format as the output of function DIR().
%
%   If you have filenames like
%
%    MISAS_DT_22_10_2009.MR.HEAD_GENERAL.6.10.2009.10.23.16.02.38.562500.9699412.DCM
%    MISAS_DT_22_10_2009.MR.HEAD_GENERAL.6.11.2009.10.23.16.02.38.562500.9699451.DCM
%    MISAS_DT_22_10_2009.MR.HEAD_GENERAL.6.1.2009.10.23.16.02.38.562500.9699056.DCM
%    MISAS_DT_22_10_2009.MR.HEAD_GENERAL.6.12.2009.10.23.16.02.38.562500.9699490.DCM
%
%   where the slice order is given by the field between the "6" and the
%   "2009", then you can use Gerardus function sortdirbynum() to order the
%   slices correctly. Example:
%
%   >> str1 = 'MISAS_DT_22_10_2009.MR.HEAD_GENERAL.6.*.DCM';
%   >> str2 = 'MISAS_DT_22_10_2009.MR.HEAD_GENERAL.6.mha';
%   >> file = sortdirbynum( str1, 5, '.' );
%   >> dcm2metaimage( str2, [], [], [], 'DCM', file );
%
% IM = DCM2METAIMAGE(...);
%
%   IM is the image volume.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2009-2010 University of Oxford
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
error( nargchk( 1, 6, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

% defaults
% if ( nargin < 2 || isempty( res ) )
%     res = [ 1, 1, 1 ];
% end
% "res" default set below
if ( nargin < 3 || isempty( scale ) )
    scale = 1;
end
if ( nargin < 4 || isempty( crop ) )
    crop = [];
end
if ( nargin < 5 || isempty( ext ) )
    ext = 'tif';
end
if ( nargin < 6 || isempty( file ) )
    file = [];
end

% check arguments
if ( scale <= 0 || scale > 1 )
    error( 'SCALE must be in (0, 1]' );
end

% get path to the data files
[ dirdata, name ] = fileparts( str );
if (isempty(dirdata))
    dirdata = '.';
end

% get list of image files
if isempty( file )
    file = dir( [ dirdata filesep '*.' ext ] );
end

% if resolution is not provided, try to read it from DICOM file
if ( nargin < 2 || isempty( res ) )
    % get DICOM header
    info = dicominfo( [ dirdata filesep file( 1 ).name ] );
    
    % if resolution is in the DICOM header, read it...
    if ( isfield( info, 'SliceThickness' ) && isfield( info, 'PixelSpacing' ) )
        res = [ info.PixelSpacing' info.SliceThickness ];
    else % ... if not, then fall back to default resolution
        res = [ 1, 1, 1 ];
    end
end

% adjust the image resolution to the scaling factor
res( 1:2 ) = res( 1:2 ) / scale;

% load first slice
frame = dicomread( [ dirdata filesep file( 1 ).name ] );

% remove alpha channel or other layers if present
frame = frame( :, :, 1 );

% init volume to load data
im = zeros( [ size( frame )*scale length( file ) ], class( frame ) );

% load slices of whole volume
im( :, :, 1 ) = imresize( frame, scale, 'bilinear' );
for I = 2:length( file )
    frame = dicomread( [ dirdata filesep file( I ).name ] );
    frame = frame( :, :, 1 );
    % resize frames
    im( :, :, I ) = imresize( frame, scale, 'bilinear' );
end

% crop volume (this is less efficient than initialising the volume as
% cropped, and then cropping each frame, but it makes the code simpler)
%
% The "+1" is necessary to follow the Seg3D convention
if ~isempty( crop )
    im = im( ...
        crop( 1 )+1:crop( 2 ), ...
        crop( 3 )+1:crop( 4 ), ...
        crop( 5 )+1:crop( 6 ) ...
        );
end

% compute offset in the coordinates of the cropped volume
if isempty( crop )
    offset = [ 0.0 0.0 0.0 ];
else
    offset = crop( [ 1 3 5 ] ) .* res;
end

% write mha file
writemetaimagefile([dirdata filesep name '.mha'], im, res, offset);

% avoid outputing the image volume unless the user has requested it
% explicitly
if ( nargout == 0 )
    varargout = [];
else
    varargout = { im };
end
