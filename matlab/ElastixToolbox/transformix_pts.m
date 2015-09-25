function xout = transformix_pts(t, x, opts)
% TRANSFORMIX_PTS  Matlab interface to the image warping program
% "transformix", but for point transformation.
%
% TRANSFORMIX_PTS is a simple interface to the command line program
% "transformix", but for point coordinates inputs, instead of images
%
%   http://elastix.isi.uu.nl/
%
% FILENAMEOUT = TRANSFORMIX_PTS(T, FILENAMEIN)
% XOUT = TRANSFORMIX_PTS(T, X)
%
%   T is a struct, or the path and name of a text file with the transform
%   parameters. Typically, this is the output of a call to elastix, the
%   registration program. See "help elastix" for details. If T is empty, no
%   transformation is applied to the input.
%
%   The input points can be provided either as a string with the path and
%   filename (FILENAMEIN), or as an array with point coordinates (X). The
%   output will have the same format as the input.
%
%   Note 1: If points are provided as an array X, we expect real world
%   coordinates. If they are provided as a file FILENAMEIN, then they can
%   be provided in real world coordinates or as voxel indices.
%
%   Note 2: Counterintuitively, the input points are expected to be in the
%   *fixed* image space, and get mapped onto the moving image. This is
%   because of the way that transforms are defined in elastix.
%
%   If the points are provided in file FILENAMEIN, the file format expected
%   by transformix is, according to section 4.2 of the Elastix Manual v4.7
%   (http://elastix.isi.uu.nl/download/elastix_manual_v4.7.pdf)
%
%     <index, point>
%     <number of points>
%     point1_x point1_y [point1_z]
%     point2_x point2_y [point2_z]
%     . . .
%
%   FILENAMEOUT, or XOUT is the output transformed points.
%
% TODO: Transformix output provides 6 decimal digits for the coordinates.
% Maybe it'd be better to read the indices, and the convert to real-world
% coordinates ourselves.
%
% ... = TRANSFORMIX_PTS(..., OPTS)
%
%   OPTS is a struct with options:
%
%     verbose: (def 0) 1 = show transformix output on screen.
%
%     outfile: path and filename to output image. If none is provided, a
%       random filename in the temp directory is created. This option is
%       ignored if the input points are given in array form.
%
% See also: transformix, elastix, elastix_read_file2param,
% elastix_write_param2file, elastix_read_reg_output.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2015 University of Oxford
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
    
    error('T must be a struct, or a filename')
    
end

% process the input points
if (isnumeric(x))
    
    % if points are provided as an array, we need to create a temp file for
    % transformix
    xfile = [tempname '.txt'];
    fid = fopen(xfile, 'w');
    if (fid == -1)
        error(['Cannot open temp file to write X: ' xfile]);
    end
    
    % write header according to transformix format
    fprintf(fid, 'point\n');
    fprintf(fid, '%d\n', size(x, 1));
    fclose(fid);
    
    % write points, using Matlab's long format precission
    dlmwrite(xfile, x, 'delimiter', ' ', '-append', 'precision', '%.15f');
    
    % we'll delete the temp file at the end of the function
    delete_xfile = true;
    
    % number and dimension of input points
    N = size(x, 1);
    D = size(x, 2);
    
elseif (ischar(x))
    
    % if user provides the points in a file, we won't delete it at the end
    % of the function
    xfile = x;
    delete_xfile = false;
    
    % number of input points
    fid = fopen(xfile, 'r');
    if (fid == -1)
        error(['Cannot read file with X points: ' xfile]);
    end
    tline = fgetl(fid); % 'point' or 'index'
    if (~strcmp(tline, 'point') && ~strcmp(tline, 'index'))
        error('First line of points file must be either "point" or "index"')
    end
    tline = fgetl(fid); % number of points
    N = str2double(tline);
    if (isnan(N))
        error('Second line of points file must be a scalar with the number of input points')
    end
    tline = fgetl(fid); % first point
    
    % dimension of input points
    D = size(str2num(tline), 2);
    fclose(fid);
    
else
    
    error('X must be an array or a filename');
    
end

%% apply transformation to image
if (opts.verbose)
    status = system(['transformix'...
        ' -tp ' tfile ...
        ' -def ' xfile ...
        ' -out ' outdir ...
        ]);
else
    [status, ~] = system(['transformix'...
        ' -tp ' tfile ...
        ' -def ' xfile ...
        ' -out ' outdir ...
        ]);
end
if (status ~= 0)
    error('Transformation failed')
end

%% read warped points as an array

% open file with the results
fid = fopen([outdir filesep 'outputpoints.txt'], 'r');
if (fid == -1)
    error(['Cannot open file to read transformed points: ' outdir ...
        filesep 'outputpoints.txt']);
end

% init array for the warped points
xout = zeros(N, D);

% read line by line (one point per line)
for I = 1:N
    
    % each line from output file looks a bit like this
    % Point   0       ; InputIndex = [ 1582 0 ]       ; InputPoint = [ 0.014450 0.000000 ]    ; OutputIndexFixed = [ 1423 -363 ]      ; OutputPoint = [ 0.012998 -0.003333 ]  ; Deformation = [ -0.001452 -0.003333 ]
    tline = fgetl(fid);
    
    % keep only the OutputPoint coordinates
    tline = regexp(tline, 'OutputPoint = \[.*?\]', 'match');
    tline = regexp(tline{1}, '\d.*\d', 'match');
    tline = tline{1};
    
    % convert to numeric format
    xout(I, :) = str2num(tline);
    
end
fclose(fid);

%% format output

if (isnumeric(x))
    
    % the output is already the numeric array, nothing to do here
    
elseif (ischar(x))
    
    % we need to create a file with the points
    
    % if the user has not provided a name for the output file
    if (isempty(opts.outfile))
        
        % we create one in the temp directory
        opts.outfile = [tempname '.txt'];
        
    end
    fid = fopen(opts.outfile, 'w');
    if (fid == -1)
        error(['Cannot open temp file to write XOUT: ' opts.outfile]);
    end
    
    % write header according to transformix format
    fprintf(fid, 'point\n');
    fprintf(fid, '%d\n', N);
    fclose(fid);
    
    % write points
    dlmwrite(opts.outfile, xout, 'delimiter', ' ', '-append', ...
        'precision', '%.15f');
    
    % we provide at the output the name of the file, not the array itself
    xout = opts.outfile;
    
end

%% clean up

% clean up temp directory
rmdir(outdir, 's')

% delete temp files
if (delete_tfile)
    elastix_delete_param_file(tfile)
end
if (delete_xfile)
    delete(xfile)
end
