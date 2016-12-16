function im_resize(filein, fileout, K)
% IM_RESIZE Resize a list of image files (.tif, .png, etc) with ImageMagick
% convert, directly file to file
%
% This function implements a fast way of resizing a lot of image files
% (.tif, .png, etc), by using the linux command line program "convert"
% running in parallel with the bash "parallel" program.
%
% IM_RESIZE(FILEIN, FILEOUT, K)
%
%   FILEIN, FILEOUT are cell vectors of the same length. Each element
%   contains the name of an input file that will be resized, and the
%   filename the result is saved to.
%
%   K is a scalar or 2-vector with the resizing factor. E.g. K=0.4 is a 40%
%   resize. If it's a two-vector, K=[K_width, K_height].

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2016 University of Oxford
% Version: 0.1.0
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
narginchk(3, 3);
nargoutchk(0, 0);

% if input is a simple filename, put it into a cell
if (~iscell(filein))
    filein = {filein};
end
if (~iscell(fileout))
    filein = {fileout};
end
if (isscalar(K))
    K = [K K];
end

% number of input files
N = length(filein);
if (N ~= length(fileout))
    error('Number of input and output files must be the same')
end

% open file to pass commands to parallel
fileParam = [tempname '.txt'];
fid = fopen(fileParam, 'w');
if (fid == -1)
    error(['Cannot open file to save commands for "parallel": ' fileParam])
end

Kreal = zeros(N, 2); % actual scaling parameter (width, height)
Dout = zeros(N, 2); % density (width, height)
for I = 1:N
    
    % check that the input and out filenames are not the same
    if (strcmp(filein{I}, fileout{I}))
        error(['Input and output filenames are the same: ' filein{I}])
    end
    
    % read header info of the input image
    info = imfinfo(filein{I});
    
    % discard thumbnail if present
    info = info(1);
    
    % compute actual scaling factor for width and height of the image
    szin = [info.Width, info.Height];
    szout = round(szin .* K);
    Kreal(I, :) = szout ./ szin;
    
    % density of input image
    Din = [info.XResolution info.YResolution];
    
    % density of output image
    Dout = round(Din .* Kreal(I, :));
    
    % units of density
    if (strcmpi(info.ResolutionUnit, 'centimeter'))
        ResolutionUnit = 'PixelsPerCentimeter';
    else
        error('ResolutionUnit not implemented')
    end
    
    % write command to parameters file
    fprintf(fid, '%dx%d,%s,%15.13f%%x%15.13f%%,%s,%s\n', ...
        Dout(1), Dout(2), ResolutionUnit, ...
        Kreal(I, 1)*100, Kreal(I, 2)*100, ...
        filein{I}, fileout{I});
    
end

% close parameter file
status = fclose(fid);
if (status == -1)
    error(['Cannot close commands file for "parallel": ' fileParam])
end

% run the resizing command
status = system(['parallel -a ' fileParam ' --colsep '','' convert -density {1} -units {2} -resize {3} {4} {5}']);
if (status ~= 0)
    error('parallel convert could not resize the images')
end

% delete temp file
delete(fileParam)
