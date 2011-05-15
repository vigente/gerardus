function matlab2vox_seg3d( path_filename )
% MATLAB2VOX_SEG3D  Convert Matlab file with Seg3D segmentation to vox
% format (for the Tarantula meshing application)
%
% MATLAB2VOX_SEG3D( PATH_FILENAME )
%
%   PATH_FILENAME is a string with the full path and filename for a Matlab
%   file created with Seg3D that contains a 3D segmentation.
%
%   This function creates a file with extension .vox and the same name in
%   a format that can be read by the Tarantula application to create
%   meshes.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2009 University of Oxford
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
error( nargchk( 1, 1, nargin, 'struct' ) );
error( nargoutchk( 0, 0, nargout, 'struct' ) );

% output filename
[ datadir, filename ] = fileparts( path_filename );
fileout = [ filename '.vox' ];

% load the segmentation
im = load( path_filename );

% get resolution for later
res = [ im.scirunnrrd.axis.spacing ];
res = res( 2:end );

% get just image data
im = squeeze( im.scirunnrrd.data );

% open output file and write header
fid = fopen( [ datadir filesep fileout ], 'wt' );
fprintf( fid, num2str( size( im ) ) );
fprintf( fid, '\n' );
fprintf( fid, num2str( res ) );
fprintf( fid, '\n' );

% write image data and close file
fprintf( fid, '%d\n', im( : ) );
fclose( fid );
