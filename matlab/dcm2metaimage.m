function varargout = dcm2metaimage( str, res, scale )
% DCM2METAIMAGE  Read a batch of DICOM files, collate them and save as a
% single MetaImage volume (one .mha and one .raw file)
%
% DCM2METAIMAGE(STR)
%
%   STR is a string with the path and file name expression. For example
%
%     STR='/home/john/data/study01/img*.tif'; % linux
%     STR='C:\data\study01\img*.tif';         % windows
%
%   will create two files
%
%     /home/john/data/study01/img.mha
%     /home/john/data/study01/img.raw
%
%     C:\data\study01\img.mha
%     C:\data\study01\img.raw
%
% DCM2METAIMAGE(STR, RES, SCALE)
%
%   RES is a 3-vector with the pixel size in the x-, y- and z-coordinates.
%   By default RES=[1 1 1];
%
%   SCALE is a factor to reduce the size of each frame (the number of
%   frames doesn't change). By default, SCALE=1.
%
% IM = DCM2METAIMAGE(...);
%
%   IM is the image volume.


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
error( nargchk( 1, 3, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

% defaults
if ( nargin < 2 )
    res = [ 1, 1, 1 ];
end
if ( nargin < 3 )
    scale = 1;
end

% check arguments
if ( scale <= 0 || scale > 1 )
    error( 'SCALE must be in (0, 1]' );
end

% get path to the data files
[ dirdata, name ] = fileparts( str );

% get list of image files
file = dir( str );

% remove wildcards from file name
name = strrep( name, '*', '' );
name = strrep( name, '?', '' );

% adjust the image resolution to the scaling factor
res( 1:2 ) = res( 1:2 ) / scale;

% load first slice
frame = dicomread( [ dirdata filesep file( 1 ).name ] );

% init volume to load data
im = zeros( [ size( frame )*scale length( file ) ], ...
    class( frame ) );

% load slices of whole volume
im( :, :, 1 ) = dicomresize( frame, scale, 'bilinear' );
for I = 2:length( file )
    frame = dicomread( [ dirdata filesep file( I ).name ] );
    % resize frames
    im( :, :, I ) = imresize( frame, scale, 'bilinear' );
end

% write mha file
WriteMhaFile( [ dirdata filesep name '.mha' ], ...
    [ size( im, 1 ), size( im, 2 ), size( im, 3 ) ], res, class( im ) )

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
