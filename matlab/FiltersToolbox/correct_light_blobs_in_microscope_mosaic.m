function im2 = correct_light_blobs_in_microscope_mosaic( im, tilesz, N )
% CORRECT_LIGHT_BLOBS_IN_MICROSCOPE_MOSAIC  Correct the colour blob created
% by the microscope's light in each tile of a mosaic, e.g. for histology.
%
% B = CORRECT_LIGHT_BLOBS_IN_MICROSCOPE_MOSAIC(A, TILESZ, N)
%
%   A is the RGB microscope image given as a 3D array.
%
%   TILESZ is a 2-vector with the (rows, cols) size of each mosaic tile.
%   Note that sometimes tiles at the bottom row and right column of the
%   mosaic are cropped. This function takes care of that automatically.
%
%   N is a 2-vector with the number of tiles in (rows, cols).

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
error( nargchk( 3, 3, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

% get image size
imsz = [ size( im, 1 ), size( im, 2 ) ];

% extract top left tile
tile = im( 1:tilesz( 1 ), 1:tilesz( 2 ), : );

% create volume for all tiles
vol = zeros( tilesz( 1 ), tilesz( 2 ), 3, prod( N ), 'double' ) * nan;

% extract all tiles; it first gets 1 column, and goes row by row, then next
% colum, etc
for J = 1:N( 2 ) % cols
    for I = 1:N( 1 ) % rows
        
        % get tile corners
        ymin = ( I - 1 ) * tilesz( 1 ) + 1;
        ymax = I * tilesz( 1 );
        xmin = ( J - 1 ) * tilesz( 2 ) + 1;
        xmax = J * tilesz( 2 );
        
        % adjust, if the tile is cropped
        if ymax > imsz( 1 )
            ymax = imsz( 1 );
        end
        if xmax > imsz( 2 )
            xmax = imsz( 2 );
        end
        
        % get tile
        aux = im( ymin:ymax, xmin:xmax, : );
        
        % stack tile in volume
        vol( 1:size( aux, 1 ), 1:size( aux, 2 ), :, ( J - 1 ) * N( 1 ) + I ) = aux;
        
    end
end

% for each tile pixel and RGB level, remove the NaNs and compute the median
tilem = double( tile ) * nan;
for I = 1:size( tile, 1 ) % rows
    for J = 1:size( tile, 2 ) % cols
        for C = 1:3 % rgb component
            
            % remove the nans
            aux = vol( I, J, C, : );
            aux = squeeze( aux( ~isnan( aux ) ) );
            
            % compute median component intensity for that pixel
            tilem( I, J, C ) = median( double( aux ) );
            
        end
    end
end

% filter image removing the median tile from each tile position
im2 = double( im ) * 0;
for J = 1:N( 2 ) % cols
    for I = 1:N( 1 ) % rows
        
        % get tile corners
        ymin = ( I - 1 ) * tilesz( 1 ) + 1;
        ymax = I * tilesz( 1 );
        xmin = ( J - 1 ) * tilesz( 2 ) + 1;
        xmax = J * tilesz( 2 );
        
        % adjust, if the tile is cropped
        if ymax > imsz( 1 )
            ymax = imsz( 1 );
        end
        if xmax > imsz( 2 )
            xmax = imsz( 2 );
        end
        
        % get tile from original image
        aux = im( ymin:ymax, xmin:xmax, : );
        
        % correct tile, and save to output image
        im2( ymin:ymax, xmin:xmax, : ) = double( aux ) ...
            - tilem( 1:size( aux, 1 ), 1:size( aux, 2 ), : );
        
    end
end

% rescale output image if necessary to make it uint8
switch class( im )
    case 'uint8'
        im2 = uint8( ( im2 + 255 ) / 2 );
end
