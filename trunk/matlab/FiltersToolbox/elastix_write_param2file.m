function filename = elastix_write_param2file(filename, param)
% elastix_write_param2file  Write struct with elastix parameters to file.
%
% FILENAME = elastix_write_param2file([], PARAM)
% FILENAME = elastix_write_param2file('', PARAM)
%
%   PARAM is a struct with parameters for elastix (the result of a
%   transform or the parameters for registration). E.g.
%
%                             Transform: 'SimilarityTransform'
%                    NumberOfParameters: 4
%                   TransformParameters: [1.0001 -7.2000e-05 4.6345 -0.2104]
%    InitialTransformParametersFileName: 'NoInitialTransform'
%                HowToCombineTransforms: 'Compose'
%                   FixedImageDimension: 2
%                  MovingImageDimension: 2
%           FixedInternalImagePixelType: 'float'
%                                     ...
%
%    A temp file created and the struct is saved in a format understood by
%    elastix. For the example above,
%
%   (Transform "SimilarityTransform")
%   (NumberOfParameters 4)
%   (TransformParameters 1.0001 -7.2000e-05 4.6345 -0.2104)
%   (InitialTransformParametersFileName "NoInitialTransform")
%   (HowToCombineTransforms "Compose")
%   (FixedImageDimension 2)
%   (MovingImageDimension 2)
%   (FixedInternalImagePixelType "float")
%   ...
%
%   FILENAME is the path and name of the text file with the parameters. The
%   filename is generated with tempname to make sure that it's unique.
%
% ... = elastix_write_param2file(FILENAME, PARAM)
%
%   This syntax allows to decide the output filename, rather than letting
%   the function produce a random name.

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
narginchk(2, 2);
nargoutchk(0, 1);

if (~isstruct(param))
    error('PARAM must be a struct with the transformation parameters')
end

% if no filename is provided, then we create a temp filename
if (isempty(filename))
    [pathstr, name] = fileparts(tempname);
    filename = [pathstr filesep 'ElastixParameters-' name '.txt'];
end

% open file for writing
fid = fopen(filename, 'w');
if (fid == -1)
    error(['Cannot open file for writing: ' filename])
end

% loop fields in the struct
fieldname = fieldnames(param);
for I = 1:length(fieldname)
    
    val = param.(fieldname{I});
    
    % write field and value to a line in the textfile
    if (ischar(val))
        fprintf(fid, '(%s "%s")\n', fieldname{I}, val);
    else
        % if field value is numeric, we need to convert to string
        val = num2str(val);
        fprintf(fid, '(%s %s)\n', fieldname{I}, val);
    end
    
end

% close file
st = fclose(fid);
if (st == -1)
    error(['Cannot close file for writing: ' filename])
end
