function [ttot, info, infoTransfdiff, imout] = transfdiffreg(transform, im, optReg, optDiff)
% TRANSFDIFFREG Transform diffusion registration of a sequence of images.
%
% TRANSFDIFFREG implements a registration algorithm for a sequence of
% images I = 1, 2, ..., N, such that each image I is aligned to its two
% neighbours I-1 and I+1. This algorithm is based on diffusion of the
% transforms between pairs of images.
%
%     |         |         |         |         |         |         |
%     |  <--->  |  <--->  |  <--->  |  <--->  |  <--->  |  <--->  |
%     |         |         |         |         |         |         |
%
% Note: Currently, only implemented for matched filter rigid registration
% (regmatchedfilt) and elastix B-spline transforms (BSplineTransform).
%
% [TTOT, INFO, INFO_TRANSFDIFF, IMOUT] = TRANSFDIFFREG(TRANSFORM, IM, OPTREG, OPTDIFF)
%
%   TRANSFORM is a string with the type of transform that will be computed
%   by the registration algorithm and diffused (more details below):
%
%     'regmatchedfilt': matched filter (i.e. cross-correlation) rigid
%        transforms.
%
%     'BSplineTransform': B-spline transforms.
%
%   IM is a sequence of 2D images, grayscale or colour. The images can be
%   in one of three formats:
%
%     * A plain array IM(R, C, I, CH), R=row, C=column, I=image index,
%       CH=colour channel.
%
%     * A vector of scimat structs.  See "help scimat" for details. IM(I)
%       has the I-th image.
%
%     * A cell vector. IM{I} has the I-th image.
%
%   OPTREG and OPTDIFF are structs with parameters for the registration and
%   diffusion parts of the algorithm, respectively (see below for details).
%
%   TTOT is a vector of N structs with the transforms (in elastix format)
%   that applied to each image forms an aligned volume.
%
%   INFO is a struct (see below for details).
%
%   INFO_TRANSFDIFF is a vector of structs. INFO_TRANSFDIFF(I) contains the
%   output information from the diffusion step in registration iteration I.
%
%   IMOUT is result of transforming each IM with the corresponding TTOT.
%   Its format (array or scimat struct) is the same as IM's.
%
% =========================================================================
%
% [TTOT, IMOUT, INFO] = TRANSFDIFFREG('regmatchedfilt', IM, OPTREG, OPTDIFF)
%
%   Diffusion of rigid transforms computed by matched filter or
%   cross-correlation registration. Global optimum.
%
%   OPTREG:
%
%     'verbose': (def false) Verbose output.
%
%     'Angle': (def (-45:45)/180*pi) Vector of angle values that the
%         matched filter will try to align adjacent slices.
%
%     'tp': (def {}) User provided intra-slice registrations. tp(I) is the
%           result of registering image I to image I+1. This is by default
%           computed internally, but with this option the user can skip the
%           slow registration process and go directly into diffusion. You
%           can obtain tp from a previous run of the algorithm in output
%           INFO.tp.
%
%     'outdir': (def tempname) If IM are provided as filenames, path to the
%        output directory. For other IM formats, it is ignored.
%
%   OPTDIFF:
%
%     'verbose': (def false) Verbose output.
%
%     'MaxIter': (def 100) Stopping criterion. The algorithm stops after
%                MaxIter diffusion iterations.
%
%     'Alpha':   (def 0.45) Diffusion coefficient. Values >0.45 may not
%                smooth high frequencies very well. Values >=0.5 cause
%                oscillations.
%
%   INFO is a struct with the output information from the top-level
%   registration algorithm.
%
%     'StopCondition': Reason why the algorithm stopped.
%
%     'tp', 'tm': Are vectors with the I -> I+1 and I -> I-1 registration
%                 results from the last registration iteration.
%
% =========================================================================
%
% [TTOT, IMOUT, INFO] = TRANSFDIFFREG('BSplineTransform', IM, OPTREG, OPTDIFF)
%
%   Diffusion of B-spline transforms. Because of B-spline non-linearity,
%   the algorithm needs to run several iterations of registration followed
%   by diffusion.
%
%   This method uses Choi and Lee (2000)'s injectivity conditions to
%   guarantee the diffused B-splines are injective and have no fold-overs.
%   Note that these conditions are conservative, so it may be possible to
%   diffuse a bit more without fold-overs.
%
%   OPTREG:
%
%     'verbose': (def false) Verbose output.
%
%     'SpatialSamplesRatio': (def 1.0) Scalar value in [0.0, 1.0], with the
%                            percentage of pixels sampled by the
%                            registration algorithm.
%
%     'MaxIter': (def 5) Stopping criterion. Maximum number of registration
%                iterations.
%
%     'MaxVal': (def 0) Stopping criterion. After a registration sweep and
%               diffusion, if the update of all B-spline coefficients is <=
%               MaxVal, then no more registration sweeps are computed and
%               the algorithm stops.
%
%     'RegParam': (def {}) A string with the path to a file, or a struct in
%                 elastix format with the registration parameters.
%
%     'mask': (def {}) An array of binary masks with the same dimensions as
%                      the images. Anything outside the binary masks will
%                      be ignored for registration. This can be used to
%                      avoid artifacts and to speed up the algorithm.
%                      Furthermore, note that masks don't need to be very
%                      precise, so performance can be improved a lot if
%                      they are provided as low resolution images (with the
%                      correct metainformation so that they overlap in real
%                      world coordinates).
%
%     't0': (def {}) A list of structs with initial transforms in elastix
%           format. Each t0(I) is applied to image I before registration.
%           This can be used to provide a pre-alignment of the images,
%           typically rigid.
%
%     'tp': See 'regmatchedfilt' above.
%
%     'tm': Like tp, but tm(I) is the result of registering image I to
%           image I-1. Note that for B-splines we cannot compute tm(I) from
%           tp(I).
%
%     'CacheFile': (def '') Name of a .mat file to save/read  internal
%                  variables to and from after each full registration
%                  sweep. This allows to kill the algorithm, and restart it
%                  without having to repeat previous registration sweeps.
%
%   OPTDIFF:
%
%     'verbose': (def false) Verbose output.
%
%     'MaxIter': (def 100) Stopping criterion. Maximum number of diffusion
%                iterations.
%
%     'Epsilon': (def 0.0) Stopping criterion. It stops when all
%                coefficient displacement components are small:
%                for all t, abs(t)<=Epsilon.
%
%     'Alpha':   (def 0.45) Diffusion coefficient. Values >0.45 may not
%                smooth high frequencies very well. Values >=0.5 cause
%                oscillations.
%
%   INFO is a struct with the output information from the top-level
%   registration algorithm.
%
%     'StopCondition': Reason why the algorithm stopped.
%
%     'tp', 'tm': Are vectors with the I -> I+1 and I -> I-1 registration
%                 results from the last registration iteration.
%
%     'MaxVal': MaxVal(I, J) is the maximum absolute value of the
%               coefficients of the I-th transform at registration
%               iteration J.
%
%     'MedVal': MedVal(I, J) is the median absolute value of the
%               coefficients of the I-th transform at registration
%               iteration J.
%
% =========================================================================

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2015-2016 University of Oxford
% Version: 0.7.1
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

%% Process input arguments

% check number of arguments
switch (transform)
    case {'regmatchedfilt', 'BSplineTransform'}
        narginchk(2, 4);
        nargoutchk(0, 4);
    otherwise
        error('Not implemented for this transform')
end

% defaults
if (nargin < 3 || isempty(optReg))
    optReg = struct;
end
if (nargin < 4 || isempty(optDiff))
    optDiff = struct;
end

% common diffusion parameters
if (~isfield(optReg, 'verbose'))
    optReg.verbose = false;
end

% list of transforms, the options they admit and their default values
optRegList = {
    'regmatchedfilt',   'verbose',             false
    'regmatchedfilt',   'Angle',               (-45:45)/180*pi
    'regmatchedfilt',   'tp',                  {}
    'regmatchedfilt',   'outdir',              tempname
    'BSplineTransform', 'verbose',             false
    'BSplineTransform', 'SpatialSamplesRatio', 1.0
    'BSplineTransform', 'MaxIter',             5
    'BSplineTransform', 'MaxVal',              0
    'BSplineTransform', 'RegParam',            {}
    'BSplineTransform', 'mask',                {}
    'BSplineTransform', 't0',                  {}
    'BSplineTransform', 'tp',                  {}
    'BSplineTransform', 'tm',                  {}
    'BSplineTransform', 'outdir',              tempname
    'BSplineTransform', 'CacheFile',           ''
    };
optDiffList = {
    'regmatchedfilt',   'Alpha',               0.45
    'regmatchedfilt',   'MaxIter',             100
    'BSplineTransform', 'Alpha',               0.45
    'BSplineTransform', 'MaxIter',             100
    'BSplineTransform', 'Epsilon',             0.0
    };
optReg = check_user_options(optRegList, transform, optReg);
optDiff = check_user_options(optDiffList, transform, optDiff);

% number of input images
if (iscell(im) || isstruct(im)) % filenames or scimats
    N = length(im);
elseif (ischar(im))
    error('IM cannot be a single filename')
else % array
    N = size(im, 3);
end

% transform-specific checks
switch (transform)
    
    case 'regmatchedfilt'
        
        % with affine transforms, registrations between images need to be
        % computed only once
        optReg.MaxIter = 1;
        
        % transform any input image into scimat format
        [im, imOutType] = imToScimat(im);
        
    case 'BSplineTransform'

        if (optReg.SpatialSamplesRatio < 0.0 ...
                || optReg.SpatialSamplesRatio > 1.0)
            error('optReg.SpatialSamplesRatio must be in [0, 1]')
        end
        
        % if registration parameters provided as file name, read the file
        % (we may need to change some of the parameters, e.g.
        % SpatialSamplesRatio, and it also makes it easier to check the
        % values of the parameters)
        if (isempty(optReg.RegParam))
            error('No registration parameters provided in optReg.RegParam')
        elseif (ischar(optReg.RegParam))
            optReg.RegParam = elastix_read_file2param(optReg.RegParam);
        end
        
        % hard-coded options for registration between two images
        optElastix.verbose = false;

        % if images are in memory instead of as files, we have to
        % create temp files for them
        [im, deleteImFiles, imOutType] = imToFiles(im);
        N = length(im);
        
        % if the user has provided masks, and they are in memory, they need to be
        % saved as temp files
        [optReg.mask, deleteMaskFiles, maskType, nVoxMask] = imToFiles(optReg.mask);
        Nmask = length(optReg.mask);
        if (Nmask > 0 && N ~= Nmask)
            error('There must be either no masks or the same number of masks (OPTREG.mask) as images (IM)')
        end
        
end

%% Main algorithm

% init number of registration iterations
RegIter = 0; % number of registration iterations

switch (transform)

    case 'BSplineTransform'
        
    % initialize total transforms applied to each image as empty
    if (isempty(optReg.t0))
        ttot = cell(1, N);
    else
        ttot = optReg.t0;
    end
    
    % if the user is running in cached mode, we set the internal state of
    % the algorithm to the latest full registration sweep
    if (exist(optReg.CacheFile, 'file') == 2)
        
        if (optReg.verbose)
            disp('***************************************************')
            disp('** READING CACHED INFO ****************************')
            disp('***************************************************')
        end
    
        % load intermediate results, except for optReg, im because the temp
        % files will have different random names
        load(optReg.CacheFile, ...
            'ttot', 'info', 'infoTransfdiff', 'transform', 'optDiff', ...
            'optRegList', 'optDiffList', ...
            'N', 'imOutType', 'optElastix', 'deleteImFiles', ...
            'deleteMaskFiles', 'nVoxMask', ...
            'maskType', 'Nmask', 'RegIter', 'J', 'imoutaux', ...
            'iterinfo', 'tp', 'tpMat',  't', 'tref')
        
    end

end

% loop registration steps
while (true)
    
    if (optReg.verbose)
        disp('***************************************************')
        disp('** REGISTRATION BLOCK *****************************')
        disp('***************************************************')
    end
    
    % increase iteration counter
    RegIter = RegIter + 1;
    if (optReg.verbose)
        disp(['** IterReg = ' num2str(RegIter) ...
            '/' num2str(optReg.MaxIter) ...
            '                                 **'])
    end
    
    %% registration block
    if (optReg.verbose)
        disp('** Register image i to image i+1                 **')
    end
    
    switch (transform)
        
        case 'regmatchedfilt'
            
            % register image to next neighbour
            if (isempty(optReg.tp))
                
                % register each image i to image i+1
                % images loop
                parfor J = 1:N-1
                    
                    if (optReg.verbose)
                        fprintf('** ... %d/%d', J, N-1)
                        tic
                    end
                    
                    % register each image J to its next neighbour J+1
                    aux(J) = regmatchedfilt(im(J+1), im(J), ...
                        optReg.Angle);
                    
                    if (optReg.verbose)
                        fprintf(' (%0.2f sec)\n', toc)
                    end
                    
                end % images loop
                info.tp = aux;
                
            % if the user has provided tp already, we are going to skip the
            % registration step. This allows for fast testing of diffusion
            % without having to recompute the registration again and again
            else % if (isempty(optReg.tp))
                
                if (optReg.verbose)
                    disp('** TP provided by user, skipping registration block')
                    info.tp = optReg.tp;
                end
                
            end % if (isempty(optReg.tp))
            
        case 'BSplineTransform'
            
            % register each image i to image i+1
            % images loop
            for J = 1:N-1
                
                if (optReg.verbose)
                    fprintf('** ... %d/%d', J, N-1)
                end
                
                % in case we only want to use a subset of the voxels within
                % the mask. We use the fixed image to set the number of
                % samped voxels for both images
                if (strcmp(optReg.RegParam.ImageSampler, 'RandomCoordinate'))
                    
                    % number of samples we are going to take for
                    % registration
                    numSamples = round(nVoxMask(J) ...
                        * optReg.SpatialSamplesRatio);
                    
                    % make sure that the number of pixels is not too low,
                    % or that we are taking more than what's in the mask
                    numSamples = max([numSamples 2000]);
                    numSamples = min([numSamples nVoxMask(J)]);
                    
                    optReg.RegParam.NumberOfSpatialSamples = numSamples;
                end
                
                % if the user has provided tp already, we are going to skip
                % the registration step. This allows for fast testing of
                % diffusion without having to recompute the registration
                % again and again
                if (isempty(optReg.tp))
                    
                    % register each image J to its next neighbour J+1
                    tic
                    [info.tp(J), imoutaux, iterinfo] = register_2_images(...
                        optReg.RegParam, ...         % registration param
                        im{J+1}, ...                 % fixed
                        im{J}, ...                   % moving
                        ttot(J+1), ...               % fixedT0
                        ttot(J), ...                 % movingT0
                        optReg.mask{J+1}, ...        % fixedMask
                        optReg.mask{J}, ...          % movingMask
                        optElastix);                 % elastix options
                    
                    % DEBUG
                    if (optReg.verbose)
                        fprintf(' (tp %0.2f sec)', toc)
                    end
                    
                else
                    
                    info.tp(J) = optReg.tp(J);
                    fprintf(' (tp provided)')
                    
                end
                
                % same as above but for tm
                if (isempty(optReg.tm))
                    
                    % register each image J to its previous neighbour J-1
                    tic
                    [info.tm(J), imoutaux, iterinfo] = register_2_images(...
                        optReg.RegParam, ...         % registration param
                        im{J}, ...                   % fixed
                        im{J+1}, ...                 % moving
                        ttot(J), ...                 % fixedT0
                        ttot(J+1), ...               % movingT0
                        optReg.mask{J}, ...          % fixedMask
                        optReg.mask{J+1}, ...        % movingMask
                        optElastix);                 % elastix options
                    
                    % DEBUG
                    if (optReg.verbose)
                        fprintf(' (tm %0.2f sec)\n', toc)
                    end
                    
                else
                    
                    info.tm(J) = optReg.tm(J);
                    fprintf(' (tm provided)\n')
                    
                end
                
            end % images loop
            
            % average transforms so that the transform J->J+1 is the
            % inverse of J+1->J. This is important because otherwise
            % diffusion doesn't converge
            %
            % we only need to correct tp, because tm will be easily
            % computed as -tp inside of transfdiff()
            for J = 1:N-1
                
                tp(J) = elastix_colon(info.tp(J), 1);
                tm(J) = elastix_colon(info.tm(J), 1);
                
                tp(J).TransformParameters ...
                    = (tp(J).TransformParameters ...
                    - tm(J).TransformParameters) * .5;
                
            end
            clear tm
            
            % if the user provided tp and tm at the input, they were valid
            % only for the first registration iteration, and need to be
            % removed for subsequent iterations
            if (~isempty(optReg.tp))
                optReg.tp = [];
            end
            if (~isempty(optReg.tm))
                optReg.tm = [];
            end
            
    end % switch transform
    
    %% diffusion block
    
    if (optReg.verbose)
        disp('***************************************************')
        disp('** DIFFUSION BLOCK ********************************')
        disp('***************************************************')
    end
    
    switch (transform)
        
        case 'regmatchedfilt'
            
            % convert elastix parameter vectors to affine matrices
            tpMat = elastix_affine_struct2matrix(info.tp);
            
            % transform diffusion of affine parameters
            optDiff.Transform = 'EulerTransform';
            [t, infoTransfdiff] = transfdiff(optDiff, tpMat);
            
            % transfer diffused parameters to elastix structs
            ttot = elastix_affine_matrix2struct(t, info.tp(1));
            
        case 'BSplineTransform'
            
            % extract transform vector from each TP and concatenate into a
            % matrix
            tpMat = cat(1, tp(:).TransformParameters);
            
            % transform diffusion of B-spline parameters
            optDiff.Transform = transform; % transform = 'BSplineTransform'
            optDiff.tbsp = tp(1); % used to validate B-spline injectivity
            [t, infoTransfdiff(RegIter)] = transfdiff(optDiff, tpMat);
            
            % template for the current B-spline transforms to apply to each
            % image
            tref(1:N) = tp(1);
            
            % transfer diffused parameters to elastix structs
            for J = 1:N
                tref(J).TransformParameters = t(J, :);
            end
            
            % add the new level of B-splines to the series of transforms to
            % apply to each image
            ttot = elastix_cat(tref, ttot);
            
            % stopping criterion: all control point translations are
            % very small
            info.MaxVal(:, RegIter) = max(abs(t), [], 2);
            info.MedVal(:, RegIter) = median(abs(t), 2);

            % save intermediate results to cache file
            if (~isempty(optReg.CacheFile))
                save(optReg.CacheFile)
            end

            if (optReg.verbose)
                disp(['** max(info.MaxVal) = ' num2str(max(info.MaxVal(:, RegIter)))])
            end
            if (max(info.MaxVal(:, RegIter)) <= optReg.MaxVal)
                info.StopCondition = 'MaxVal';
                break
            end

    end
    
    % stop criterion
    if (RegIter == optReg.MaxIter)
        info.StopCondition = 'MaxIter';
        break
    end
    
end % end registration step loop

% apply transform results to the input images
if (nargout > 3)
    
    if (optReg.verbose)
        disp('***************************************************')
        disp('** OUTPUT IMAGE TRANSFORM BLOCK *******************')
        disp('***************************************************')
    end
    
    % create output directory if it doesn't exist already and the output is
    % going to be filenames
    if (iscell(im) && ~exist(optReg.outdir))
        mkdir(optReg.outdir);
    end

    % if the output is not filenames, then we will delete the files, as the
    % output will be provided in memory
    deleteImOutFiles = ~strcmp(imOutType, 'char') ...
        && ~strcmp(imOutType, 'filename');

    for I = 1:N
        
        if (optReg.verbose)
            fprintf('** ... %d/%d\n', I, N)
            tic
        end
        
        % if the output is going to be given as filenames, then we want
        % control over where those files are saved to
        if (strcmp(imOutType, 'filename') || strcmp(imOutType, 'char'))

            % filename we are going to provide to the user
            [~, filename, ext] = fileparts(im{I});
            optTransformix.outfile = [optReg.outdir filesep filename ext];
            
        else
            
            optTransformix.outfile = '';

        end
        
        % at this point, "im" can be:
        % * vector of scimat (regmatchedfilt)
        % * cell vector of filenames (BSplineTransform)
        if (isstruct(im))

            % transform input image
            aux = transformix(ttot(I), im(I));
            
        elseif (iscell(im))
            
            % transform input image
            aux = transformix(ttot(I), im{I}, optTransformix);
        
        else
            
            error('Internal check fail: IM should not have this type')
            
        end
        
        % make the output image match the format of the input image
        if (isstruct(im) && strcmp(imOutType, 'char'))
            
            error('TODO: Not implemented')
            
        elseif (iscell(im) && strcmp(imOutType, 'char'))
            
            imout = aux;
            
        elseif (isstruct(im) && strcmp(imOutType, 'filename'))
            
            imout{I} = optTransformix.outfile;
            scimat_save(imout{I}, aux);
            
        elseif (iscell(im) && strcmp(imOutType, 'filename'))
            
            imout{I} = aux;
            
        elseif (isstruct(im) && strcmp(imOutType, 'scimat'))
            
            imout(I) = aux;
            
        elseif (iscell(im) && strcmp(imOutType, 'scimat'))
            
            imout(I) = scimat_load(aux);
            delete(aux);
            
        elseif (isstruct(im) && strcmp(imOutType, 'array'))
            
            imout(:, :, I, :) = aux.data;
            
        elseif (iscell(im) && strcmp(imOutType, 'array'))
            
            error('TODO: Not implemented')

        end
            
    end
    
    % delete output files if the output is provided in memory, not files
    if (exist(optReg.outdir) && deleteImOutFiles)
        rmdir(optReg.outdir);
    end
    
end

% delete temp files
switch (transform)
        
    case 'BSplineTransform'

        % if temp files had to be created to pass images to elastix, we can
        % delete them now
        if (deleteImFiles)
            for I = 1:N
                delete(im{I})
            end
        end
        if (deleteMaskFiles && ~isempty(optReg.mask))
            for I = 1:N
                delete(optReg.mask{I})
            end
        end
        
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% nested functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check_user_options
%
% check invalid options and set defaults
function opt = check_user_options(optList, transform, opt)

% list of options allowed for this transform
idx = cellfun(@(s) strcmp(s, transform), optList(:, 1));
optList = optList(idx, :);

% list of options in the user input
names = fieldnames(opt);

% check whether any option provided by the user is not used for this
% transform
for I = 1:length(names)
    
    if (~ismember(names{I}, optList(:, 2)))
        warning([transform ' does not use option ' inputname(3) ...
            '.' names{I}])
    end
    
end

% set options not provided by the user to their defaults
isProvided = ismember(optList(:, 2), names);
for I = find(~isProvided)'
    opt = setfield(opt, optList{I, 2}, optList{I, 3});
end

end

%% imToScimat()
%
% convert any input image format (filenames, arrays or scimat) to scimat
function [im, imType] = imToScimat(im)

% check the type and number of the input images
if (isempty(im)) % empty input
    
    % nothing to do
    return
    
elseif (ischar(im)) % single filename
    
    imType = 'char';
    im = scimat_load(im);
    
elseif (iscell(im)) % list of filenames
    
    imType = 'filename';
    for I = 1:length(im)
        aux(I) = scimat_load(im{I});
    end
    im = aux;
    
elseif (isstruct(im))
    
    imType = 'scimat';
    
else % plain array

    imType = 'array';

    for I = 1:size(im, 3)
        aux(I) = scimat_im2scimat(im(:, :, I, :));
    end
    im = aux;
    
end

end

%% imToFiles()
%
% convert any input format to filenames, creating temp files if necessary
%
% * Inputs
%
%    im: input images (filenames, scimats or arrays)
%
% * Outputs:
%
%    im: cell vector with paths and filenames to images
%    deleteFiles: whether we need to delete the image files when we finish
%                 registering them (i.e. whether they are temp files)
%    nVox: number of voxels ~= 0 (only useful for binary masks)
function [im, deleteFiles, imType, nVox] = imToFiles(im)

% check the type and number of the input images
if (isempty(im)) % empty input (typically, no masks provided)
    
    imType = '';
    deleteFiles = false;
    nVox = 0;
    return
    
elseif (ischar(im)) % single filename

    imType = 'char';
    deleteFiles = false;
    im = {im};
    
    % get metainformation
    if (nargout > 3) % user wants to know number of voxels ~= 0
        scimat0 = scimat_load(im{1});
        nVox = nnz(scimat0.data);
    end
    
elseif (iscell(im)) % list of filenames

    imType = 'filename';
    deleteFiles = false;
    N = length(im);
    
    % get metainformation (we assume that all slices have the same)
    if (nargout > 3) % user wants to know number of voxels ~= 0
        for I = 1:N
            scimat0 = scimat_load(im{I});
            nVox(I) = nnz(scimat0.data);
        end
    end
    
elseif (isstruct(im)) % list of scimat structs

    imType = 'scimat';

    % number of images
    N = length(im);
    
    % if user wants to know number of voxels ~= 0
    if (nargout > 3)
        for I = 1:N
            nVox(I) = nnz(im(I).data);
        end
    end
    
    % create a temp file for each input image
    deleteFiles = true;
    filename = cell(1, N);
    for I = 1:N
        filename{I} = [tempname '.mha'];
        scimat_save(filename{I}, im(I));
    end
    im = filename;
    
else % plain array
    
    imType = 'array';

    % number of images
    N = size(im, 3);
    
    % if user wants to know number of voxels ~= 0
    if (nargout > 3)
        if (size(im, 4) ~= 1)
            warning('User wants to count voxels~=0, so we assume binary mask, but it has more than one channel')
        end
        for I = 1:N
            nVox(I) = nnz(im(:, :, I, :));
        end
    end
    
    % create a temp file for each input image
    deleteFiles = true;
    filename = cell(1, N);
    for I = 1:N
        filename{I} = [tempname '.mha'];
        scimat_save(filename{I}, scimat_im2scimat(im(:, :, I, :)));
    end
    im = filename;
    
end

end

%% register_2_images: local function to register two images
function [t, imout, iterinfo] = register_2_images(regParam, fixed, moving, fixedT0, movingT0, fixedMask, movingMask, optElastix)
tic
if (iscell(fixedT0)) % no initial transform
    
    % the only way to pass empty initial transorms is as a cell with an
    % empty, but if we pass them to elastix, we need to remove the cell
    fixedT0 = fixedT0{1};
    optElastix.fMask = fixedMask;
    
else % initial transform
    
    % apply initial transform
    fixed = transformix(fixedT0, fixed);
    optElastix.fMask = transformix(fixedT0, fixedMask);
    
end
if (iscell(movingT0)) % no initial transform
    
    % the only way to pass empty initial transorms is as a cell with an
    % empty, but if we pass them to elastix, we need to remove the cell
    optElastix.mMask = movingMask;
    
else % initial transform

    % we apply the initial transform explictly to the moving image
    % and mask, because it's substantially faster than providing a list of
    % transforms in optElastix.t0
    optElastix.t0 = '';
    optElastix.mMask = transformix(movingT0, movingMask);
    moving = transformix(movingT0, moving);

end

% register both images
[t, imout, iterinfo] ...
    = elastix(regParam, fixed, moving, optElastix);

% concatenate new level of B-splines and previous initial transforms.
% Note: because of elastix internals, if we apply t0 explicitly before
% registration, then t has to be concatenated BEFORE t0, not after
t = elastix_cat(t, movingT0);

end

%% dbg_plot_reg_results
% debug function to plot the registration result
function dbg_plot_reg_results(imf, imm, fT0, mT0, mT)

% plot input images
subplot(2, 1, 1)
imagesc(imfuse(imf, imm))

% transform fixed image with its initial transform
imf = transformix(fT0, imf);

% apply the registration output transform to the moving image
immreg = transformix(mT, imm);

% transform moving image with its initial transform
imm = transformix(mT0, imm);

% plot registration results
subplot(2, 1, 2)
imagesc(imfuse(imf, immreg))

end
