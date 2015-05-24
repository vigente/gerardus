function elastix_delete_param_file(filename)
% elastix_delete_param_file  Delete nested transform files.
%
% elastix_delete_param_file deletes a transform parameter file, and all
% nested transforms it points to.
%
% elastix_delete_param_file(FILENAME)
%
%   FILENAME is a string with the path and filename to a transform
%   parameter file. This parameter file may point to another previous
%   transform in field InitialTransformParametersFileName.
%
% See also: elastix, transformix, elastix_write_param2file,
% elastix_read_file2param.

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

if (~ischar(filename))
    error('FILENAME must be a string')
end

% open file
fid = fopen(filename);
if (fid == -1)
    error(['File cannot be open: ' filename])
end

% find InitialTransformParametersFileName line in parameters file 
% (InitialTransformParametersFileName "NoInitialTransform")
while (~feof(fid))
   s = fgetl(fid);
   isInitialTransformFound ...
       = ~isempty(strfind(s, 'InitialTransformParametersFileName'));
   if (isInitialTransformFound)
       break
   end
end

% close file
fclose(fid);

% delete the file
delete(filename)

% if no nested transforms, return
if (~isInitialTransformFound)
    return
end

% get path to the nested transform
filename = regexp(s, '".*"', 'match');
filename = filename{1}(2:end-1);

% if nested transform, delete that one too
if (~strcmp(filename, 'NoInitialTransform'))
    elastix_delete_param_file(filename);
end
