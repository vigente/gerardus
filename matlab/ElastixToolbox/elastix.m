function [t, movingReg, iterInfo] = elastix(regParam, fixed, moving, opts)
% ELASTIX  Matlab interface to the image registration program "elastix".
%
% ELASTIX is a simple interface to the command line program "elastix"
%
%   http://elastix.isi.uu.nl/
%
% Command-line "elastix" is a powerful image registration program, but it
% requires that input images are provided as files, and returns the results
% (transform parameters, registered image, etc) as files in a directory.
% This is not so convenient when using Matlab, so this interface function
% allows to pass images either as filenames or as image arrays, and
% transparently takes care of creating temporary files and directories,
% reading the result to Matlab variables, and cleaning up afterwards.
%
% [T, MOVINGREG, ITERINFO] = ELASTIX(REGPARAM, FIXED, MOVING, OPTS)
%
%   REGPARARM are the registration parameters for elastix, given either as
%   a struct or a string with the path and name of a text file, e.g.
%   '/path/to/ParametersTranslation2D.txt'.
%
%   FIXED, MOVING are the images to register. They can be given as:
%
%     * file names (e.g. 'im1.mha', 'im2.tif')
%     * scimat structs (see "help scimat" for details)
%     * image arrays (e.g. im = checkerboard(10))
%
%   Images can have one channel (grayscale) or more (e.g. RGB colour
%   images). In the case of multi-channel images, we use the "multi-image"
%   feature in elastix (see section 6.1.1 "Image registration with multiple
%   metrics and/or images" in the elastix manual), with all channels
%   contributing equally to the combined cost function. FIXED and MOVING
%   must have the same number of channels.
%
%   OPTS is a struct where the fields contain options for the algorithm:
%
%     verbose: (def 0) If 0, it hides all the elastix command output to the
%              screen.
%
%     outfile: (def '') Path and name of output image file to save
%              registered image to. This option is ignored when MOVING is
%              an image array.
%
%     t0:      (def '') Struct or filename with an initial transform (see T
%              below for format). t0 is applied to the image before the
%              registration is run. Field
%              t0.InitialTransformParametersFileName allows to provide
%              another transform that will be applied before t0, and so on
%              iteratively. Note: It seems that t0 does not work well in
%              conjuction with masks, and it causes "Too many samples map
%              outside moving image buffer".
%
%     fMask:   (def '') Mask for the fixed image. Only voxels == 1 are
%              considered for the registration.
%
%     mMask:   (def '') Mask for the moving image. Only voxels == 1 are
%              considered for the registration.
%
%   T is a struct with the contents of the parameter transform file
%   (OUTDIR/TransformParameters.0.txt), e.g.
%
%                             Transform: 'TranslationTransform'
%                    NumberOfParameters: 2
%                   TransformParameters: [-2.4571 0.3162]
%    InitialTransformParametersFileName: 'NoInitialTransform'
%                HowToCombineTransforms: 'Compose'
%                   FixedImageDimension: 2
%                  MovingImageDimension: 2
%           FixedInternalImagePixelType: 'float'
%          MovingInternalImagePixelType: 'float'
%                                  Size: [2592 1944]
%                                 Index: [0 0]
%                               Spacing: [1 1]
%                                Origin: [0 0]
%                             Direction: [1 0 0 1]
%                   UseDirectionCosines: 'true'
%                  ResampleInterpolator: 'FinalBSplineInterpolator'
%        FinalBSplineInterpolationOrder: 3
%                             Resampler: 'DefaultResampler'
%                     DefaultPixelValue: 0
%                     ResultImageFormat: 'png'
%                  ResultImagePixelType: 'unsigned char'
%                   CompressResultImage: 'false'
%
%   The format of the transforms is:
%
%     'TranslationTransform': [tx ty].
%     'EulerTransform':       [rotation(rad) tx ty].
%     'BSplineTransform':     [dx1 dx2 ... dxN dy1 dy2 ... dyN]. These
%                             values represent the displacement of the
%                             control grid points, not their coordinates.
%
%   Note: In Elastix x -> rows, y -> columns, the opposite of Matlab's
%   convention. If you use elastix_bspline_grid or
%   elastix_bspline_grid2param, though, B-spline parameters are transposed
%   to follow Matlab's convention.
%
%   If T is a sequence of transforms, e.g.
%
%     tb.Transform = 'BSplineTransform';
%     tb.InitialTransformParametersFileName = ta;
%
%     ta.Transform = 'EulerTransform';
%     ta.InitialTransformParametersFileName = 'NoInitialTransform';
%
%   The EulerTransform (ta) is applied first to the image, followed by the
%   BSplineTransform (tb).
%
%   MOVINGREG is the result of registering MOVING onto FIXED. MOVINGREG is
%   the same type as MOVING (i.e. image array, scimat struct or path and
%   filename). In the path and filename case, an image file will we created
%   with path and name MOVINGREG=PARAM.outfile. If PARAM.outfile is not
%   provided or empty, then the registered image is deleted and
%   MOVINGREG=''. For multi-channel images (e.g. RGB colour images),
%   MOVINGREG will correspond to the first channel only, as this is the
%   output provided by elastix. To obtain the multi-channel result, it is
%   necessary to apply the transform to MOVING with transformix().
%
%   ITERINFO is a struct with the details of the elastix optimization
%   (OUTDIR/IterationInfo.0.R0.txt), e.g.
%
%        ItNr: [10x1 double]
%      Metric: [10x1 double]
%    stepSize: [10x1 double]
%    Gradient: [10x1 double]
%        Time: [10x1 double]
%
%
% See also: transformix, elastix_read_file2param, elastix_write_param2file,
% elastix_read_reg_output.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014-2015 University of Oxford
% Version: 0.5.4
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
narginchk(3, 4);
nargoutchk(0, 3);

if (isempty(regParam))
    error('REGPARAM must be a .txt file or a struct with the registration parameters for elastix')
end

% defaults
if (nargin < 4 || isempty(opts) || ~isfield(opts, 'verbose'))
    % capture the elastix output so that it's not output on the screen
    opts.verbose = 0;
end
if (nargin < 4 || isempty(opts) || ~isfield(opts, 'outfile'))
    % no name provided for the output file
    opts.outfile = '';
end
if (nargin < 4 || isempty(opts) || ~isfield(opts, 't0'))
    % no initial transform provided
    opts.t0 = '';
end
if (nargin < 4 || isempty(opts) || ~isfield(opts, 'fMask'))
    % no fixed image mask provided
    opts.fMask = '';
end
if (nargin < 4 || isempty(opts) || ~isfield(opts, 'mMask'))
    % no moving image mask provided
    opts.mMask = '';
end

% input parameters and initial transform can be either filename or struct.
% We want to have both, as we need the file for elastix, and the struct for
% some checks in this function
[regParam, paramfile, delete_paramfile] = filename_struct_param(regParam);
[opts.t0, t0file, delete_t0file] = filename_struct_param(opts.t0);

% if images are given as arrays instead of filenames, we need to create 
% temporary files with the images so that elastix can work with them
[fixedfile, delete_fixedfile] = create_temp_file_if_array_image(fixed, '.mha');
[movingfile, delete_movingfile] = create_temp_file_if_array_image(moving, '.mha');

% check that fixed and moving images have the same number of channels
if (size(fixedfile, 1) ~= size(movingfile, 1))
    error('Fixed and moving image must have the same number of channels')
end

% check that registration parameters are correct for multichannel images
if (size(fixedfile, 1) > 1)
    if (~strcmp(regParam.Registration, 'MultiMetricMultiResolutionRegistration'))
        error('Multi-channel image registration requires regParam.Registration = ''MultiMetricMultiResolutionRegistration''.')
    end
end

% if masks are given as arrays instead of filenames, we need to create 
% temporary files with the images so that elastix can work with them
[fMaskfile, delete_fMaskfile] = create_temp_file_if_array_image(opts.fMask, '.mha');
[mMaskfile, delete_mMaskfile] = create_temp_file_if_array_image(opts.mMask, '.mha');

if (ischar(moving))
    
    % if moving image is given as a filename, then we'll save the
    % registered image to the path provided by the user. If the user has
    % not provided an output name, then the result image will not be
    % returned. This is useful when we only need the transformation
    if (~isfield(opts, 'outfile'))
        opts.outfile = '';
    end
    
else
    
    % if moving image was provided as an image, the output will be an
    % image array in memory, and we can ignore the output directory, as it
    % will not be saved to a file
    opts.outfile = '';
    
end

% create temp directory for output
tempoutdir = tempname;
success = mkdir(tempoutdir);
if (~success)
    error(['Cannot create temp directory for registration output: ' tempoutdir])
end

% create text string with the command to run elastix
comm = 'elastix ';
for I = 1:size(fixedfile, 1)
    comm = [...
        comm ...
        ' -f' num2str(I-1) ' "' fixedfile(I, :) '"' ...
        ' -m' num2str(I-1) ' "' movingfile(I, :) '"' ...
        ];
end
comm = [...
    comm ...
    ' -out ' tempoutdir ...
    ' -p ' paramfile
    ];
if (~isempty(opts.t0))
    comm = [...
        comm ...
        ' -t0 ' t0file ...
        ];
end
if (~isempty(opts.fMask))
    comm = [...
        comm ...
        ' -fMask "' fMaskfile '"' ...
        ];
end
if (~isempty(opts.mMask))
    comm = [...
        comm ...
        ' -mMask "' mMaskfile '"' ...
        ];
end

% run elastix
if (opts.verbose)
    status = system(comm);
else
    % hide command output from elastix
    [status, ~] = system(comm);
end
if (status ~= 0)
    error('Registration failed')
end

if (ischar(moving))
    
    % check that the registration produced the expected output image
    regfile = dir([tempoutdir filesep 'result.0.' regParam.ResultImageFormat]);
    if (isempty(regfile))
        error('Elastix did not produce an output registration')
    end
    
    % read elastix result
    [t, ~, iterInfo] = elastix_read_reg_output(tempoutdir);
    
    % check that the output file extension matches the extension of the
    % output filename, and give a warning if they don't match
    if (~isempty(opts.outfile))
        
        [~, ~, regfile_ext] = fileparts(regfile.name);
        [~, ~, outfile_ext] = fileparts(opts.outfile);
        if (~strcmpi(regfile_ext, outfile_ext))
            warning('Gerardus:BadFileExtension', ...
                [ 'Elastix produced a ' regfile_ext ...
                ' image, but user asked to save to format ' outfile_ext])
        end
        
        % the first output argument will be the path to the output file,
        % rather than the image itself
        movingReg = opts.outfile;
        
        % move the result image to the directory requested by the user
        movefile([tempoutdir filesep regfile.name], movingReg);
        
    else
        
        % the first output argument will be the path to the output file, rather
        % than the image itself
        movingReg = '';
        
    end
    
    
else
    
    % read elastix result
    [t, movingReg, iterInfo] = elastix_read_reg_output(tempoutdir);
    
end

% If multi-channel registration, use transformix to get other channel
% images (elastix only returns a single channel image after multi-channel registration)
% For now, only when elastix is given a structure.
if ndims(movingfile)>1 && (~ischar(moving))
    for i=2:size(movingfile,1)
        % Use tranformix to register other channel
        [status, ~] = system(['transformix -in ' movingfile(i,:) ' -out ' tempoutdir  ' -tp ' tempoutdir filesep 'TransformParameters.0.txt']);

        % the result file can be in many different formats
        resultfile = dir([tempoutdir filesep 'result.' t.ResultImageFormat]);
        if (isempty(resultfile))
            error(['No image result file for additional channels: ' outdir filesep 'result.' t.ResultImageFormat])
        end
        % read the image
        movingRegTemp = scimat_load([tempoutdir filesep resultfile.name]);
        movingReg.data= cat(5, movingReg.data, movingRegTemp.data);
        delete_image([tempoutdir filesep resultfile.name]);
    end
end
    
% delete temp files and directories
if (delete_fixedfile)
    delete_image(fixedfile)
end
if (delete_movingfile)
    delete_image(movingfile)
end
if (delete_paramfile)
    elastix_delete_param_file(paramfile)
end
if (delete_t0file)
    elastix_delete_param_file(t0file)
    if (ischar(t.InitialTransformParametersFileName))
        % if the transform is returned as a filename, the only one that is
        % going to be preserved is the top transform. The nested ones will
        % be deleted, so we cannot return filenames for them
        t.InitialTransformParametersFileName = 'NoInitialTransform';
    end
end
if (delete_fMaskfile)
    delete(fMaskfile)
end
if (delete_mMaskfile)
    delete(mMaskfile)
end
rmdir(tempoutdir, 's')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% nested functions

%% create_temp_file_if_array_image()
%
% check whether an input image is provided as an array or as a path to a
% file. If the former, then save the image to a temp file so that it can be
% processed by elastix
%
% file: filename or array with the input image
% extOut: extension of the temp files we use internally
function [filename, delete_tempfile] = create_temp_file_if_array_image(file, extTemp)

if (isempty(file)) % no filename or array provided by the user
    
    % this option is useful for the fixed and moving masks, when they are
    % not provided by the user, and thus file==''
    filename = '';
    delete_tempfile = false;
    
    return
    
end

% find number of channels
if (ischar(file)) % user has provided a filename

    [~, ~, extIn] = fileparts(file);
    
    % read file metadata block (we need to know whether the image is
    % grayscale or colour    
    switch (extIn) % type of input file
        
        case {'.mha', '.mhd'} % MetaImage format from ITK
            
            % open file to read
            fid = fopen(file, 'r');
            if (fid <= 0)
                error(['Cannot open file: ' file])
            end
            
            % initialize variable to save number of channels
            nchannel = [];
            
            % process text header, and stop if we get to the raw data
            while (true)
                
                % read text header line
                tline = fgetl(fid);
                
                % if end of text header stop reading
                % Warning: this only works if the first voxel=0. Otherwise, we
                % would read some voxels as part of a header line. To avoid it,
                % we now break after finding ElementDataFile
                if (tline(1) == 0), break, end
                
                % find location of "=" sign
                idx = strfind(tline, '=');
                if isempty(idx), break, end
                
                switch getlabel(tline, idx)
                    case 'ElementNumberOfChannels'
                        nchannel = getnumval(tline, idx);
                        break
                end
                
            end
            
            % close file
            fclose(fid);
            
        otherwise % '.png', '.tif'
    
            % read image header
            info = imfinfo(file);
            
            % number of image channels
            if (isfield(info, 'SamplesPerPixel'))
                nchannel = info.SamplesPerPixel;
            elseif (isfield(info, 'BitDepth'))
                nchannel = info.BitDepth / 8;
            else
                warning('File does not provide information about number of channels. Assuming 1 channel')
                nchannel = 1;
            end
            
    end
    
elseif (isstruct(file)) % image is a scimat struct
    
        nchannel = size(file.data, 5);
        
else % image is an array
    
        nchannel = size(file, 3);
    
end % find number of channels

% deal with every type of input image (array, filename, 1 or more
% channels...). If it's a grayscale file, then we don't need to do
% anything, because we can just use that file. Otherwise, we read a file
% into a scimat struct, so that we can later write the channels to
% separate files
if (nchannel == 1 && ischar(file)) % grayscale image provided in a file
    
    % input is a filename already, and image is grayscale. Thus we don't
    % need to create a temp file
    delete_tempfile = false;
    filename = file;
    
    return

elseif (ischar(file)) % filename, more than one channel
    
    % load image
    file = scimat_load(file);
    
elseif (~isstruct(file)) % array, but no metainformation (no scimat)
    
    % convert the array to a scimat struct
    file = scimat_im2scimat(reshape(file, ...
        size(file, 1), size(file, 2), 1, 1, nchannel));
    
end

% if we reach this point, the image is in a scimat struct. We are going to
% save each channel into a separate file (whether its 1 channel or more)
    
% create temp files to save the channels of the image
delete_tempfile = true;
filename = char(zeros(nchannel, length(tempname) + length(extTemp), 'uint8'));
scichan = file;
if size(scichan.data,3)==1  % Truncate time dimension but avoid discarding z-axis if present 
    scichan.axis = scichan.axis(1:2);
else
    scichan.axis = scichan.axis(1:3);
end
for I = 1:nchannel
    
    % create a temp filename for the image
    [pathstr, name] = fileparts(tempname);
    filename(I, :) = [pathstr filesep name extTemp];
    
    % extract one channel
    scichan.data = file.data(:, :, :, :, I);
    
    % write the channel to temp file
    scimat_save(filename(I, :), scichan);
    
end

end % function ... = create_temp_file_if_array_image(...)

%% filename_struct_param()
%
% parameters can be provided as a struct or filename. This function deals
% with both cases. If it's a struct, it writes it to a temp file. If it's a
% filename, it loads the file as a struct.
%
% Input:
%
%   param: filename or scimat struct
%
% Output:
%
%   param: always a scimat struct
%
%   filename: always a filename
%
%   delete_tempfile: whether the param file has to be deleted when we are
%          done with it. That's the case when we had to create a temp file.
function [param, filename, delete_tempfile] = filename_struct_param(param)

if (isempty(param))

    % nothing provided. This is useful for T0, as the user will not always
    % provide an initial transform
    filename = '';
    delete_tempfile = false;
    
elseif (isstruct(param))

    % create a temp file for the parameter file
    filename = elastix_write_param2file('', param);
    delete_tempfile = true;
    
elseif (ischar(param))
    
    % we are going to use the file provided by the user, and will not
    % delete it afterwards
    filename = ['"' param '"'];
    delete_tempfile = false;

    % we read the file
    param = elastix_read_file2param(param);
    
else
    
    error('REGPARAM must be a struct or file name')
    
end

end

%% delete_image()
%
% delete one or more files containing the channels of an image
function delete_image(filename)

for I = 1:size(filename, 1)
    delete(filename(I, :))
end

end

% helper function to parse lines of .mha header
function s = getlabel(s, idx)
s = strtrim(s(1:idx-1));
end
function n = getnumval(s, idx)
n = str2num(s(idx+1:end));
end
