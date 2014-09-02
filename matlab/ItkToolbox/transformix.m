function imout = transformix(t, im, opts)
% transformix  Matlab interface to the image warping program "transformix".
%
% transformix is a simple interface to the command line program
% "transformix"
%
%   http://elastix.isi.uu.nl/
%
% FILENAMEOUT = transformix(T, FILENAMEIN)
% IMOUT = transformix(T, IMIN)
%
%   T is a struct, or the path and name of a text file with the transform
%   parameters. Typically, this is the output of a call to elastix, the
%   registration program. See help to elastix for details.
%
%   The input image can be provided either as a string with the path and
%   filename (FILENAMEIN) or as an image array (IMIN). The output will have
%   the same format (i.e. a path to an output file, or an image array). The
%   input image can be grayscale (1 channel) or colour RGB (3 channels).
%
%   FILENAMEOUT or IMOUT is the output transformed image.
%
% ... = transformix(..., OPTS)
%
%   OPTS is a struct with options:
%
%     verbose: (def 0) 1 = show transformix output on screen.
%
%     outfile: path and filename to output image. If none is provided, a
%       random filename in the temp directory is created. This option is
%       ignored if the input/output images are given in array form.
%
% See also: elastix, elastix_read_file2param, elastix_write_param2file,
% elastix_read_reg_output.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.2.1
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
narginchk(2, 3);
nargoutchk(0, 1);

% defaults
if (nargin < 3 || isempty(opts))
    % capture the elastix output so that it's not output on the screen
    opts.verbose = 0;
end
if (~isfield(opts, 'verbose') || isempty(opts.verbose))
    opts.verbose = 0;
end
if (~isfield(opts, 'outfile'))
    opts.outfile = '';
end

% if the transformation parameters are provided as a struct, create a text
% file with them, so that we can pass it to transformix
if (isstruct(t))
    
    tfile = elastix_write_param2file([], t);
    delete_tfile = true;
    
elseif (ischar(t))
    
    tfile = t;
    delete_tfile = false;
    
else
    
    error('T must be a struct, or a path and filename')
    
end

% extension of output file
if (isfield(t, 'ResultImageFormat'))

    % if user has provided an explicit extension for the output file,
    % that's what elastix is going to produce for the output
    ext = ['.' t.ResultImageFormat];
    
else % output format has to be deduced
    
    if (ischar(im)) % input image given as path and filename
    
        % if the user has not provided an extension for the output file, it
        % will have the same as the input
        [~, ~, ext] = fileparts(im);
        
    else % input image given as an array
        
        % opts.outfile is ignored
        % for temp file with the image, we'll use the PNG format
        ext = '.png';
        
    end
    
end

% if image is provided as a path and filename
if (ischar(im))

    % clarify nomeclature
    imfile = im;
    
    % if user has not provided an output filename, we generate a random
    % one. If he has provided one, we use it
    if (isempty(opts.outfile))

        opts.outfile = [tempname ext];
        
    end
    
    % read image header
    info = imfinfo(imfile);
    
    switch (info.ColorType)
        case 'truecolor'
            
            % load image
            imaux = imread(imfile);
            
            % Elastix only processes the first channel, and ignores the
            % rest. Thus, we need to create 3 temp files, one per channel.            
            % We overwrite the name of the original colour file, as we are
            % not going to use it. We'll need to delete the temp files when
            % we finish with them
            delete_imfile = true;
            clear imfile
            imfile{1} = [tempname ext];
            imfile{2} = [tempname ext];
            imfile{3} = [tempname ext];
            imwrite(imaux(:, :, 1), imfile{1});
            imwrite(imaux(:, :, 2), imfile{2});
            imwrite(imaux(:, :, 3), imfile{3});
            
        case 'grayscale'
            
            % input is a filename already, and image is grayscale. Thus we
            % don't need to create a temp file
            delete_imfile = false;
            
        otherwise
            
            error('Image in input file is neither RGB nor grayscale')
    end
    
else % if the image is provided as an image array, we write it to a temp 
     % so that it can be read by transformix
    
    % if image is colour
    if (size(im, 3) == 3)

        % Elastix only processes the first channel, and ignores the
        % rest. Thus, we need to create 3 temp files, one per channel.
        % We'll need to delete the temp files when we finish with them
        delete_imfile = true;
        imfile{1} = [tempname ext];
        imfile{2} = [tempname ext];
        imfile{3} = [tempname ext];
        imwrite(im(:, :, 1), imfile{1});
        imwrite(im(:, :, 2), imfile{2});
        imwrite(im(:, :, 3), imfile{3});
        
    elseif (size(im, 3) == 1) % image is grayscale

        % only one channel, so we simply need to save the image array to a
        % temp file that can be deleted after transformix processes it
        delete_imfile = true;
        imfile = [tempname ext];
        imwrite(im, imfile);
        
    else
        
        error('IM must be grayscale (1 channel) or RGB (3 channels)')
        
    end
    
end

% init output image
switch (t.ResultImagePixelType)
    case 'unsigned char'
        imout = zeros([t.Size([2 1]) length(imfile)], 'uint8');
    otherwise
        error(['ResultImagePixelType = ' t.ResultImagePixelType ...
            ' not implemented'])
end

% apply transformation to each image channel separatedly
if (iscell(imfile))
    for CH = 1:length(imfile)
        imout(:, :, CH) = warp_image(tfile, imfile{CH}, ext, opts);
    end
else
    imout = warp_image(tfile, imfile, ext, opts);
end

% format output
if (ischar(im)) % input image given as path and filename

    % write image to output file
    imwrite(imout, opts.outfile);
    
    % return to user the path and name of result file
    imout = opts.outfile;
    
end

% clean up temp files
if (delete_tfile)
    delete(tfile)
end
if (delete_imfile)
    if (iscell(imfile))
        for CH = 1:length(imfile)
            delete(imfile{CH})
        end
    else
        delete(imfile)
    end
end

end

% warp_image: auxiliary function to warp one channel of the image
function imout = warp_image(tfile, imfile, ext, opts)

% create temp directory for the output
outdir = tempname;
mkdir(outdir)

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

% check that the output transformed image was created
resultfile = dir([outdir filesep 'result' ext]);
if (isempty(resultfile))
    error(['transformix did not produce an output image in directory ' outdir])
end
resultfile = [outdir filesep resultfile.name];

% read result image into memory
imout = imread(resultfile);

% clean up temp directory
rmdir(outdir, 's')

end
