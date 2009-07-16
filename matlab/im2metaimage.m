function varargout = im2metaimage( str, res, scale, crop )
% IM2METAIMAGE  Read a batch of image files, collate them and save as a
% single MetaImage "file" (actually, one .mha and one .raw file)
%
% IM2METAIMAGE(STR)
%
%   STR is a string with the path and output mha header file name. For
%   example,
%
%     STR='/home/john/data/study01/study01.mha'; % linux
%     STR='C:\data\study01\study01.mha';         % windows
%
%   It is assumed that the MetaImage files (.mha and .raw) will be created
%   from the TIFF files found in the target directory.
%
% IM2METAIMAGE(STR, RES, SCALE, CROP)
%
%   RES is a 3-vector with the pixel size in the x-, y- and z-coordinates.
%   By default RES=[1 1 1];
%
%   SCALE is a scalar factor to reduce the size of each frame (the number
%   of frames doesn't change). By default, SCALE=1.
%
%   CROP is a vector with values [xmin, xmax, ymin, ymax, zmin, zmax].
%   These values allow to crop the data volume. It follows the same
%   convention as Seg3D [1]. For example, if xmin=0, xmax=3, then voxels 1,
%   2 and 3 are selected.
%   Note that the output MetaImage file will have an offset header so that
%   the data retains its true coordinates.
%
% IM = IM2METAIMAGE(...);
%
%   IM is the image volume.
%
% [1] http://www.sci.utah.edu/cibc/software/seg3d.html


% Copyright © 2009 University of Oxford
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
error( nargchk( 1, 4, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

% defaults
if ( nargin < 2 )
    res = [ 1, 1, 1 ];
end
if ( nargin < 3 )
    scale = 1;
end
if ( nargin < 4 )
    crop = [];
end

% check arguments
if ( scale <= 0 || scale > 1 )
    error( 'SCALE must be in (0, 1]' );
end

% get path to the data files
[ dirdata, name ] = fileparts( str );

% get list of image files
file = dir( [ dirdata filesep '*.tif' ] );

% adjust the image resolution to the scaling factor
res( 1:2 ) = res( 1:2 ) / scale;

% load first slice
frame = imread( [ dirdata filesep file( 1 ).name ] );

% remove alpha channel or other layers if present
frame = frame( :, :, 1 );

% init volume to load data
im = zeros( [ size( frame )*scale length( file ) ], class( frame ) );

% load slices of whole volume
im( :, :, 1 ) = imresize( frame, scale, 'bilinear' );
for I = 2:length( file )
    frame = imread( [ dirdata filesep file( I ).name ] );
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
WriteMhaFile( [ dirdata filesep name '.mha' ], ...
    [ size( im, 1 ), size( im, 2 ), size( im, 3 ) ], res, class( im ), ...
    offset )

% write raw file
WriteRawFile( [ dirdata filesep name '.raw' ], ...
    im, res, class( im ) )

% avoid outputing the image volume unless the user has requested it
% explicitly
if ( nargout == 0 )
    varargout = [];
else
    varargout = { im };
end
