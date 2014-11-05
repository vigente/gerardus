function [tParam, regParam] = histology_blockface_reg(pathstrh, fileh, pathstrbf, filebf, regParam, polymask, ellipmask, opts)
% HISTOLOGY_BLOCKFACE_REG  Register deformed histology volume to blockface.
%
% Correct a histology volume obtained by intraslice registration with
% histology_blockface_reg(). The steps in this function are:
%
%   1) Each "valid" histology slice is registered to the corresponding
%   blockface slice (typically with a B-spline, but the user can provide
%   any registration parameters). "Valid" slices are pre-selected by the
%   user, and are those where the registration will give a good result,
%   because there are no artifacts in the wax, the tissue has no tears,
%   etc.
%
%   2) B-spline transforms for "invalid" slices are interpolated and
%   extrapolated from the valid ones.
%
%   3) Each B-spline coefficient is smoothed so that it changes slowly
%   between each slice and its neighbours. This preserves the structural
%   integrity of the histology reconstruction.
%
% [TPARAM, REGPARAM] = histology_blockface_reg(PATHSTRH, FILEH, PATHSTRBF, FILEBF, REGPARAM, POLYMASK, ELLIPMASK)
%
%   PATHSTRH and PATHSTRBF are strings with the paths to the histology and
%   blockface images, respectively.
%
%   FILEH, FILEBF are the lists of histology and blockface files, obtained
%   with dir(). Each file contains a slice. FILEH(i) corresponds to
%   FILEBF(i).
%
%   REGPARAM is a struct (or a string with the full path to a file) with
%   the registration parameters for elastix. Typically, you want this
%   registration to use a multi-level B-spline.
%
%   POLYMASK, ELLIPMASK are two binary masks with the polygon highlighting
%   the wax blockface (POLYMASK) and an ellipse around the heart at it's
%   biggest (ELLIPMASK). The masks can be produced with
%   blockface_create_masks().
%
% ... = histology_blockface_reg(..., OPTS)
%
%   OPTS is a struct with parameters for the algorithm.
%
%     'verbose':   (def 0) Verbose output of the registration algorithm.
%
%     't0':        (def []) Initial transform to apply to histology slices
%                  before registration.
%
%     'IdxIgnore': (def []) Indices of slices where instead from a
%                  registration, the transform is obtained by
%                  interpolation/extrapolation from the other slices.
%                  This parameter can be provided as a vector of indices or
%                  as a boolean vector.

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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

DEBUG = 0;

% check arguments
narginchk(7, 8);
nargoutchk(0, 2);

% defaults
if (nargin < 8 || isempty(opts) || ~isfield(opts, 'verbose'))
    opts.verbose = 0;
end
if (isempty(opts) || ~isfield(opts, 't0'))
    error('OPTS.t0 expected');
end
if (~isfield(opts, 'IdxIgnore'))
    opts.IdxValid = true(1, length(fileh));
elseif (~islogical(opts.IdxValid))
    % if user provides indices as integers, convert to logical vector so
    % that it's easier to get a list of invalid slices
    aux = false(1, length(fileh));
    aux(opts.IdxValid) = true;
    opts.IdxValid = aux;
end
if (isempty(pathstrh))
    pathstrh = '.';
end
if (isempty(pathstrbf))
    pathstrbf = '.';
end

if (~isempty(opts.t0))
    if (ischar(opts.t0))
        opts.t0 = elastix_read_file2param(opts.t0);
    end
    if ((length(opts.t0) == 1))
        opts.t0 = opts.t0(ones(1, length(fileh)));
    end
end

szMask = size(ellipmask);
if (any(szMask ~= size(polymask)))
    error('ELLIPMASK and POLYMASK must be the same size')
end

% if registration parameters are provided as a filename, convert to struct
if (ischar(regParam))
    regParam = elastix_read_file2param(regParam);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial alignment of histology to blockface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we use 10% of the voxels in the slice
regParam.NumberOfSpatialSamples = round(numel(ellipmask) / 10);

% iterate blockface images
for I = find(opts.IdxValid)

    if (opts.verbose)
        disp(['I = ' num2str(I)])
    end
    
    % load blockface and histology
    imbf = imread([pathstrbf filesep filebf(I).name]);
    imh = imread([pathstrh filesep fileh(I).name]);
    
    if (any(szMask ~= size(imbf)))
        error(['Blockface image is not the same size as masks: ' ...
            pathstrbf filesep filebf(I).name])
    end
    
    % plot slices
    if (DEBUG)
        subplot(2, 1, 1)
        imagesc(imbf)
        subplot(2, 1, 2)
        imagesc(imh)
    end
    
    % preprocessing of blockface and histology to prepare them for
    % registration
    [imh, imbf] = histology_blockface_preprocessing(imh, imbf, ...
        polymask, ellipmask);
    
    % plot slices
    if (DEBUG)
        subplot(2, 1, 1)
        imagesc(imbf)
        axis([1 size(imbf, 2) 1 size(imbf, 1)])
        subplot(2, 1, 2)
        imagesc(imh)
    end
    
    % we are going to use the user-provided initial transform for the
    % registration
    optsI.t0 = opts.t0(I);
    
    % plot initial alignment
    if (DEBUG)
        imh0 = transformix(optsI.t0, imh);
        subplot(2, 1, 2)
        imagesc(imh0)
        axis([1 size(imbf, 2) 1 size(imbf, 1)])
        imshowpair(imbf, imh0)
    end
    
    % registration of histology on blockface: Similarity
    optsI.verbose = 0;
    optsI.fMask = ellipmask;
    tic
    [tParam(I), imhreg, iterinfo] = elastix(...
        regParam, imbf, imh, optsI);
    toc
    
    % plot result of initial registration
    if (DEBUG)
        subplot(2, 1, 1)
        imagesc(imhreg)
        axis([1 size(imbf, 2) 1 size(imbf, 1)])
        subplot(1, 1, 1)
        imshowpair(imbf, imhreg)
    end
    
    
end

% DEBUG: check result
if (DEBUG)
    for I = 1:length(filebf)
        
        disp(['I = ' num2str(I)])
        
        % load blockface and histology
        imbf = imread([pathstrbf filesep filebf(I).name]);
        imh = imread([pathstrh filesep fileh(I).name]);
        
        imh0 = transformix(tParam(I), imh);
        imshowpair(imbf, imh0)
        title(['I = ' num2str(I)])
        pause
    end
end

% extract transform parameters to matrix
tParamMat = cat(1, tParam(:).TransformParameters);

% smooth out parameters so that we can fit the histology to the
% blockface without breaking the slice correspondence we already have
for I = 1:tParam(1).NumberOfParameters
    
    % interpolate/extrapolate values of the sections where the user has
    % signaled that the registration fails
    tParamMat(:, I) = interp1(find(opts.IdxValid), tParamMat(opts.IdxValid, I), ...
        1:length(fileh), 'linear', 'extrap');
    
    % smooth the values of the B-spline coefficients
    tParamMat(:, I) = smooth(tParamMat(:, I), 'rlowess');
    
end

% overwrite transform parameters with the smoothed values
for I = 1:length(filebf)
    tParam(I).TransformParameters = tParamMat(I, :);
end

end
