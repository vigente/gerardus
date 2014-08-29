function blockface_correct_perspective(indir, files55, tform, outdir)
% blockface_correct_perspective  Correct perspective of Bewster angle
% blockface images.
%
% blockface_correct_perspective(INDIR, FILES55, TFORM, OUTDIR)
%
%   INDIR is a string with the directory where the input files are kept.
%
%   FILES55 is the result of a dir() command with a list of the files to
%   correct.
%
%   TFORM is a projective2d object with the transform correction for the
%   frames.
%
%   OUTDIR is a string with the directory for output images. INDIR and
%   OUTDIR must be different to avoid overwriting the input files.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

% check arguments
narginchk(4, 4);
nargoutchk(0, 0);

% safety to prevent us from overwriting the source files
if (strcmp(indir, outdir))
    error('Input and output directory are the same, and files would be overwritten')
end

% if destination directory doesn't exist, create it
if (isempty(dir(outdir)))
    ok = mkdir(outdir);
    if (~ok)
        error(['Cannot create output directory: ' outdir])
    end
end

% correct perspective
for I = 1:length(files55)
    
    % load image
    im55 = imread([indir filesep files55(I).name]);
    
    % make output image the same size and coordinates as input image
    ra = imref2d(size(im55));
    
    % correct perspective
    im55to90 = imwarp(im55, tform, 'cubic', ...
        'OutputView', ra, 'FillValues', double(min(im55(:))));
    
%     % plot results
%     subplot(2, 1, 1)
%     imagesc(im55)
%     subplot(2, 1, 2)
%     imagesc(im55to90)
    
    % save result
    imwrite(im55to90, [outdir filesep files55(I).name]);
    
end
