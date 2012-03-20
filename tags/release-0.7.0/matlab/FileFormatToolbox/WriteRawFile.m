function WriteRawFile(filename, img)
% WRITERAWFILE  Write the binary part of a MetaImage file (.raw)
%
% WRITERAWFILE(FILENAME, IMG, RESOLUTION)
%
%   FILENAME is a string with the path and name of the output file.
%
%   IMG is a 3D matrix that contains the image volume. Note that it is
%   expected that x-coordinates are in the first dimension (rows),
%   y-coordinates are in the second dimension (columns). Meanwhile,
%   the Matlab convention is that x-coordinates are in columns and
%   y-coordinates are in rows.
%
%   See also: WriteMhaFile, writemetaimagefile.

% Author(s): Ramon Casero <rcasero@gmail.com> and Vicente Grau
% Copyright © 2009-2010 University of Oxford
% Version: 0.2.0
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
error(nargchk(2, 3, nargin, 'struct'));
error(nargoutchk(0, 0, nargout, 'struct'));

% data type to be written
data_type = class(img);

% write image to binary file
fid=fopen(filename, 'w');
if(fid<=0) 
    printf('Impossible to open file %s\n', img_filename);
end
fwrite(fid, img, data_type);
fclose(fid);
