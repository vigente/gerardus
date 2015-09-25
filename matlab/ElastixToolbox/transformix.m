function sciout = transformix(t, scimat, opts)
% TRANSFORMIX  Matlab interface to the image warping program "transformix".
%
% TRANSFORMIX is a simple interface to the command line program
% "transformix"
%
%   http://elastix.isi.uu.nl/
%
% FILENAMEOUT = TRANSFORMIX(T, FILENAMEIN)
% IMOUT = TRANSFORMIX(T, IMIN)
% SCIOUT = TRANSFORMIX(T, SCIMAT)
%
%   T is a struct, or the path and name of a text file with the transform
%   parameters. Typically, this is the output of a call to elastix, the
%   registration program. See "help elastix" for details. If T is empty, no
%   transformation is applied to the input.
%
%   The input image can be provided either as a string with the path and
%   filename (FILENAMEIN), as an image array (IMIN) or as a scimat struct
%   (SCIMAT). The output will have the same format (i.e. a path to an
%   output file, an image array or a scimat struct). The input image can be
%   grayscale (1 channel) or colour RGB (3 channels). The input image can
%   have metadata (i.e. pixel size and offset) if provided in the image
%   file header (e.g. .mha/.mhd, .png, .tif formats) or in the scimat
%   struct (see "help scimat" for details).
%
%   FILENAMEOUT, IMOUT or SCIOUT is the output transformed image.
%
% ... = TRANSFORMIX(..., OPTS)
%
%   OPTS is a struct with options:
%
%     verbose: (def 0) 1 = show transformix output on screen.
%
%     outfile: path and filename to output image. If none is provided, a
%       random filename in the temp directory is created. This option is
%       ignored if the input/output images are given in array form.
%
%     AutoDefaultPixelValue: (def false) If true, it overrides
%       t.DefaultPixelValue and estimates the typical colour of the image
%       to fill in newly created background pixels. Note that different
%       values are estimated per channel, so it works for coloured
%       backgrounds.
%
% See also: elastix, elastix_read_file2param, elastix_write_param2file,
% elastix_read_reg_output.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014-2015 University of Oxford
% Version: 0.3.1
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

% defaults
if (nargin < 3 || isempty(opts))
    opts = struct;
end
if (~isfield(opts, 'verbose') || isempty(opts.verbose))
    % capture the elastix output so that it's not output on the screen
    opts.verbose = 0;
end
if (~isfield(opts, 'outfile'))
    opts.outfile = '';
end
if (~isfield(opts, 'AutoDefaultPixelValue'))
    opts.AutoDefaultPixelValue = false;
end

% override t.DefaultPixelValue?
if (opts.AutoDefaultPixelValue)
    t.DefaultPixelValue = 0;
end

% preprocess transform
if (isempty(t))
    
    % if the transform is empty, we return the input image
    sciout = scimat;
    return
    
elseif (ischar(t))
    
    % if the transform is provided as a file, read it into a struct so that
    % we can access its fields, and possibly modify them
    t = elastix_read_file2param(t);

end

% extension of output file
if (isfield(t, 'ResultImageFormat'))

    % if user has provided an explicit extension for the output file,
    % that's what elastix is going to produce for the output
    ext = ['.' t.ResultImageFormat];
    
else % output format has to be deduced
    
    if (ischar(scimat)) % input image given as path and filename
    
        % if the user has not provided an extension for the output file, it
        % will have the same as the input
        [~, ~, ext] = fileparts(scimat);
        
    else % input image given as an array
        
        % opts.outfile is ignored
        % default output file format
        ext = '.mha';
        
    end
    
end

% input image format
if (ischar(scimat))
    imType = 'char';
elseif (isstruct(scimat))
    imType = 'scimat';
else
    imType = 'array';
end

% input image preprocessing depending on its type
switch (imType) 
    
    case 'char' % if image is provided as a path and filename

        % if user has not provided an output filename, we generate a random
        % one. If he has provided one, we use it
        if (isempty(opts.outfile))
            
            opts.outfile = [tempname ext];
            
        end
        
        % read image
        scimat = scimat_load(scimat);
    
    case 'array' % if image is provided as plain array, without metadata
    
        % convert to scimat format
        scimat = scimat_im2scimat(scimat);
    
end

% number of channels
nchannel = size(scimat.data, 5);

% name for a temp filename to save and transform each channel
tmpfilename = [tempname ext];

% auxiliary variable to keep one channel at a time
scich = scimat;

% loop channels
for I = 1:nchannel
    
    % extract channel from full image
    scich.data = scimat.data(:, :, :, :, I);
    
    % smart estimation of background colour for newly added pixels
    if (opts.AutoDefaultPixelValue)
        
        % we assume as background the typical intensity voxel for this
        % channel
        t.DefaultPixelValue = median(scich.data(:));
        
    end
    
    % save channel to temp file
    scimat_save(tmpfilename, scich);
    
    % apply transform to channel
    sciaux = warp_image(t, tmpfilename, ext, opts);
    
    % aggregate channel
    if (I == 1)
        sciout = sciaux;
    else
        sciout.data(:, :, :, :, I) = sciaux.data;
    end
end

% format output
switch (imType) 
    
    case 'char' % input image given as path and filename
        
        % write image to output file
        scimat_save(opts.outfile, sciout);
        
        % return to user the path and name of result file
        sciout = opts.outfile;
        
    case 'array' % input image given as plain array
        
        sciout = sciout.data;
end

% clean up temp file
delete(tmpfilename)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% nested functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% warp_image: auxiliary function to warp one channel of the image
function sciout = warp_image(t, imfile, ext, opts)

% create temp directory for the output
outdir = tempname;
mkdir(outdir)

% process the transform parameters
if (isstruct(t))
    
    % if the transformation parameters are provided as a struct, create a
    % temporary text file with them, so that we can pass it to transformix
    tfile = elastix_write_param2file([], t);
    delete_tfile = true;
    
elseif (ischar(t))
    
    % if the transform is provided as a file, we will pass it to
    % transformix, so we don't need to do anything else here
    tfile = t;
    delete_tfile = false;
    
else
    
    error('T must be a struct, or a path and filename')
    
end

% apply transformation to image
if (opts.verbose)
    status = system(['transformix'...
        ' -tp ' tfile ...
        ' -in ' imfile ...
        ' -out ' outdir ...
        ]);
else
    [status, ~] = system(['transformix'...
        ' -tp ' tfile ...
        ' -in ' imfile ...
        ' -out ' outdir ...
        ]);
end
if (status ~= 0)
    error('Transformation failed')
end

% read output image
sciout = scimat_load([outdir filesep 'result' ext]);

% clean up temp directory
rmdir(outdir, 's')

% delete temp files
if (delete_tfile)
    elastix_delete_param_file(tfile)
end

end
