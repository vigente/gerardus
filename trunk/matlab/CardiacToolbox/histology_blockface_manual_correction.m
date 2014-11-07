function [tParam, ptsh, ptsbf] = histology_blockface_manual_correction(pathstrh, fileh, pathstrbf, filebf, idx, opts)
% HISTOLOGY_BLOCKFACE_MANUAL_CORRECTION  Manual landmark correction of
% 3D histology volume to match blockface volume.
%
% HISTOLOGY_BLOCKFACE_MANUAL_CORRECTION is used after
% histology_intraframe_reg() to provide a rough correction of the histology
% volume. The user selects a few slices, and the GUI allows to click
% corresponding landmarks in the blockface and histology (when finished
% with a slice, click on "Close Control Point Selection Tool"). Similarity
% transforms are computed for those slices, and interpolated/extrapolated
% for the rest.
%
% [TPARAM, PTSH, PTSBF] = HISTOLOGY_BLOCKFACE_MANUAL_CORRECTION(PATHSTRH, FILEH, PATHSTRBF, FILEBF, IDX)
%
%   PATHSTRH and PATHSTRBF are strings with the paths to the histology and
%   blockface images, respectively.
%
%   FILEH, FILEBF are the lists of histology and blockface files, obtained
%   with dir(). Each file contains a slice. FILEH(i) corresponds to
%   FILEBF(i).
%
%   IDX is a vector with the indices of the frames that are going to be
%   used for the correction. These are the frames where the user will click
%   corresponding landmarks between the histology and blockface.
%
%   TPARAM is a struct array with the transform parameters for each
%   histology file, in elastix format.
%
%   PTSH, PTBF are two cell arrays with the coordinates of the landmarks
%   clicked by the user. PTSH{i}, PTSBF{i} contain the landmarks for the
%   i-th slice.
%
% ... = HISTOLOGY_BLOCKFACE_MANUAL_CORRECTION(..., OPTS)
%
%  OPTS is a struct with optional parameters:
%
%    'outdir': (def []) Output directory. If this parameter is provided,
%              then the tParam transforms are applied to the histology
%              images and the results saved to that directory.
%
% See also: histology_intraframe_reg, histology_blockface_reg.

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

% check arguments
narginchk(4, 6);
nargoutchk(0, 3);

if (length(fileh) ~= length(filebf))
    error('There must be the same number of files in FILEH and FILEBF')
end

% defaults
if (nargin < 5 || isempty(idx))
    % if the user doesn't pick up any frames, use all frames for manual
    % correction
    idx = 1:length(fileh);
end
if (isempty(pathstrh))
    pathstrh = '.';
end
if (isempty(pathstrbf))
    pathstrbf = '.';
end
if (nargin < 6 || ~isfield(opts, 'outdir'))
    opts.outdir = [];
end

% init cell array to store landmarks
ptsbf = cell(1, length(fileh));
ptsh = cell(1, length(fileh));

% ask user to enter corresponding landmarks in selected frames
for I = idx(:)'
    
    % load images
    imbf = imread([pathstrbf filesep filebf(I).name]);
    imh = imread([pathstrh filesep fileh(I).name]);

    % extend histogram to dynamic range for better visualization
    minbf = double(min(imbf(:)));
    maxbf = double(max(imbf(:)));
    imbf = uint8((double(imbf) - minbf) / (maxbf - minbf) ...
        * double(intmax(class(imbf))));
    minh = double(min(imh(:)));
    maxh = double(max(imh(:)));
    imh = uint8((double(imh) - minh) / (maxh - minh) ...
        * double(intmax(class(imbf))));
    
    % find matched landmarks in both images
    [ptsbf{I}, ptsh{I}] = cpselect(imbf, imh, 'Wait', true);
    
end

% compute similarity transform between previous frames. Note that we have
% to switch blockface and histology because of how elastix defines the
% image transforms vs point transforms
for I = idx(:)'
    
    tform(I) = fitgeotrans(ptsbf{I}, ptsh{I}, 'similarity');
    
end

% convert transform matrices to parameter vectors
param = elastix_affine_tform2param(tform, 'similarity');

% create matrix with parameters for all slices. We use linear interpolation
% to avoid oscillation
paramH = zeros(length(fileh), 4);
paramH(:, 1) = interp1(idx, param(idx, 1), 1:length(fileh), 'linear');
paramH(:, 2) = interp1(idx, param(idx, 2), 1:length(fileh), 'linear');
paramH(:, 3) = interp1(idx, param(idx, 3), 1:length(fileh), 'linear');
paramH(:, 4) = interp1(idx, param(idx, 4), 1:length(fileh), 'linear');

% compute transform structs of transforms in elastix format
for I = 1:length(fileh)
    
    % load histology image
    imh = imread([pathstrh filesep fileh(I).name]);
    
    % load info of blockface image
    infobf = imfinfo([pathstrbf filesep filebf(I).name]);
    
    % create the transform parameters struct
    tParam(I) = tParamSimilarity(imh, infobf, paramH(I, :));
    
end

% if user provides output directory, we apply the transforms to the
% histology images
if (~isempty(opts.outdir))
    
    % apply similarity correction to histology files
    for I = 1:length(fileh)
        
        % apply transform to image
        opts.outfile = [opts.outdir filesep fileh(I).name];
        transformix(tParam(I), [pathstrh filesep fileh(I).name], opts);
        
    end
    
end

end

% auxiliary function to generate the similarity transform in elastix format
% for one image
function tParam = tParamSimilarity(imh, infobf, params)

tParam.Transform = 'SimilarityTransform';
tParam.NumberOfParameters = 4;
tParam.TransformParameters = params;
tParam.InitialTransformParametersFileName = 'NoInitialTransform';
tParam.HowToCombineTransforms = 'Compose';
tParam.FixedImageDimension = 2;
tParam.MovingImageDimension = 2;
tParam.FixedInternalImagePixelType = 'float';
tParam.MovingInternalImagePixelType = 'float';
tParam.Size = [infobf.Width infobf.Height];
tParam.Index = [0 0];
tParam.Spacing = [1 1];
tParam.Origin = [0 0];
tParam.Direction = [1 0 0 1];
tParam.UseDirectionCosines = 'true';
tParam.CenterOfRotationPoint = [0 0];
tParam.ResampleInterpolator = 'FinalBSplineInterpolator';
tParam.FinalBSplineInterpolationOrder = 1;
tParam.Resampler = 'DefaultResampler';
tParam.DefaultPixelValue = mode(imh(:));
tParam.ResultImageFormat = 'png';
tParam.ResultImagePixelType = 'unsigned char';
tParam.CompressResultImage = 'false';

end

