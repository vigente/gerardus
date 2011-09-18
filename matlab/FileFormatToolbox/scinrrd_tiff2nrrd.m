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
%   tiffread(), e.g.
%
%    >> stack - tiffread('file.tif');
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

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
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
error(nargchk(1, 1, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

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
    
    % data volume size
    nrrd.axis(I).size = size(nrrd.data, I);
    
    % left edge of first voxel
    nrrd.axis(I).min = nrrd.axis(I).spacing / 2;
    
    % left edge of last voxel
    nrrd.axis(I).max = (size(nrrd.data, I) - 1) * nrrd.axis(I).spacing;
    
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
