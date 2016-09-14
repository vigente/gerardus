function scimat2 = scimat_resize3(scimat, sz, sigma)
% SCIMAT_RESIZE3 Resize a 3D scimat image.
%
% This function requires the Gerardus binary resize3DImage, built into
% directory "programs" by the Gerardus CMake build system. The binary has
% to be in the path too, which can be done running script
% add_gerardus_paths.m.
%
% Note that the function writes the scimat struct to a temp file for
% resize3DImage. If you have the image already in a file format that ITK
% can read (e.g. .mha), it will be faster to run resize3DImage directly on
% the file, rather than using this Matlab interface.
%
% SCIMAT2 = SCIMAT_RESIZE3(SCIMAT, SZ, SIGMA)
%
%   SCIMAT, SCIMAT2 are the input and output scimat structs that contain
%   the image to be resized (see "help scimat" for details).
%
%   SZ is a 3-vector with the size of the output SCIMAT2 in [row col slice]
%   order.
%
%   SIGMA is a 3-vector with the standard deviation (in pixel units, [row
%   col slice] order) of the Gaussian low-pass filter used for
%   anti-aliasing. The more decimation, the larger SIGMA should be to avoid
%   aliasing artifacts. On the other hand, too much blurring will degrade
%   the result unnecessarily. By default, SIGMA = size(SCIMAT.data) ./
%   size(SCIMAT2.data), but often it's better to have less blurring.
%   
% 
% See also: scimat_resample, scimat_downsample.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2016 University of Oxford
% Version: 0.1.1
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
narginchk(2, 3);
nargoutchk(0, 1);

% check that the input image is 3D
if (ndims(scimat.data) ~= 3)
    error('SCIMAT must contain a 3D image')
end
if (length(sz) ~= 3)
    error('SZ must be a 3-vector')
end

% defaults
if (nargin < 3 || isempty(sigma))
    sigma = size(scimat.data) ./ sz;
end

% convert from row/column to x/y to match the format of resize3DImage
sz(1:2) = sz([2 1]);
sigma(1:2) = sigma([2 1]);
    
% temp names for internal use files to resize the image
infile = [tempname '.mha'];
outfile = [tempname '.mha'];

% save input image to file
scimat_save(infile, scimat);

% resize the image
err = system(['resize3DImage --sigmaInVoxels', ...
    ' --sigx ' num2str(sigma(1)) ...
    ' --sigy ' num2str(sigma(2)) ...
    ' --sigz ' num2str(sigma(3)) ...
    ' ' num2str(sz(1)) ... % sx
    ' ' num2str(sz(2)) ... % sy
    ' ' num2str(sz(3)) ... % sz
    ' ' infile ...
    ' -o ' outfile ...
    ]);
if (err)
    error('resize3DImage failed')
end

% read output
scimat2 = scimat_load(outfile);

% delete temp images
delete(infile)
delete(outfile)
