function [t, tParam, iterInfo, regParam] = blockface_intraframe_reg(pathstr, files)
% blockface_intraframe_reg  Register consecutive frames in a list of
% blockface images.
%
% blockface_intraframe_reg takes a list of files with blockface images, and
% registers each file(I) to file(I-1), using a similarity transformation to
% account for small translation and zoom-in/zoom-out movements.
%
% [T, TPARAM, ITERINFO, REGPARAM] = blockface_intraframe_reg(PATHSTR, FILES)
%
%   PATHSTR is the full path to the files. If empty, PATHTOFILES='.'.
%
%   FILES is the result of a dir() command, e.g. dir('*_55_*.png'). The
%   function expects a list of blockface images.
%
%   T is a matrix where row I contains the transform parameters of the
%   registration from FILES(I) to FILES(I-1).
%
%   TPARAM is a struct array. TPARAM(I) has all the details of the
%   transform from FILES(I) to FILES(I-1).
%
%   ITERINFO is a struct array. ITERINFO(I) has the optimization measures
%   corresponding to TPARAM(I). See elastix for details.
%
%   REGPARAM is a struct with the registration parameters generated
%   internally by the function.
%
% % See also: elastix, elastix_read_reg_output.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.4.1
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

DEBUG = 0;

% check arguments
narginchk(2, 2);
nargoutchk(0, 4);

if (length(files) < 2)
    error('At least two files are needed for registration')
end

% defaults
if (isempty(pathstr))
    pathstr = '.';
end

% create a struct with the registration parameters
regParam = generate_registration_parameters([pathstr filesep files(1).name]);

% expand list of files names to clarify notation
fixed = files(1:end-1);
moving = files(2:end);

% memory allocation for output
switch (regParam.Transform)
    case 'SimilarityTransform'
        t = zeros(length(files), 4);
    otherwise
        error('Transform not implemented')
end

% register 2nd frame onto 1st frame to obtain a param struct that we can
% use as reference
[tParam, iterInfo] = frame_registration(...
    [pathstr filesep fixed(1).name], ...
    [pathstr filesep moving(1).name], ...
    regParam, DEBUG);
t(2, :) = tParam.TransformParameters;

% Allocate memory for the rest
tParam(1:length(files)) = tParam;
iterInfo(1:length(files)) = iterInfo;

% frame 1 is the reference and doesn't get registered to anything
tParam(1).TransformParameters(:) = 0;
iterInfo(1).ItNr(:) = 0;
iterInfo(1).Metric(:) = 0;
iterInfo(1).stepSize(:) = 0;
iterInfo(1).Gradient(:) = 0;
iterInfo(1).Time(:) = 0;

% iterate images that we have to register
% note: even though elastix uses parallel processing, it is slightly faster
% (factor of ~0.93) using parfor to run 4 registrations on 1 processor at
% the same time than 1 registration on 4 processors
for I = 2:length(fixed)
    
    % display frame number
    if (DEBUG)
        disp(['Registering frame ' moving(I).name ' -> ' fixed(I).name])
    end
    
    % register I+1 frame onto I frame
    [tParam(I+1), iterInfo(I+1)] = frame_registration(...
        [pathstr filesep fixed(I).name], ...
        [pathstr filesep moving(I).name], ...
        regParam, DEBUG);
    
    % transformation parameters
    t(I+1, :) = tParam(I+1).TransformParameters;
    
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frame_registration()
%
%   Translate registration of two frames.
function [t, iterInfo] = frame_registration(fixed, moving, paramfile, DEBUG)

% load two consecutive frames
im1 = imread(fixed);
im2 = imread(moving);

% convert to grayscale
im1 = rgb2gray(im1);
im2 = rgb2gray(im2);

% match 2nd frame's histogram to 1st frame's histogram
im2 = imhistmatch(im2, im1);

% plot image 1 and difference with image 2
if (DEBUG)
    subplot(2, 1, 1)
    hold off
    imagesc(im1)
    subplot(2, 1, 2)
    hold off
    imshowpair(im1, im2, 'Scaling', 'joint');
end

% register corrected 2nd frame onto 1st frame
param.verbose = DEBUG;
[t, im2_reg, iterInfo] = elastix(paramfile, im1, im2, param);

% plot result
if (DEBUG)
    subplot(2, 1, 1)
    hold off
    imshowpair(im1, im2, 'Scaling', 'joint');
    subplot(2, 1, 2)
    hold off
    imshowpair(im1, im2_reg, 'Scaling', 'joint');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate_registration_parameters()
%
%   Create a struct with the parameters necessary to find 
function regParam = generate_registration_parameters(file)

info = imfinfo(file);

regParam.FixedInternalImagePixelType = 'float';
regParam.MovingInternalImagePixelType = 'float';
regParam.FixedImageDimension = 2;
regParam.MovingImageDimension = 2;
regParam.Registration = 'MultiResolutionRegistration';
regParam.FixedImagePyramid = 'FixedRecursiveImagePyramid';
regParam.MovingImagePyramid = 'MovingRecursiveImagePyramid';
regParam.Interpolator = 'BSplineInterpolator';
regParam.ResampleInterpolator = 'FinalBSplineInterpolator';
regParam.Resampler = 'DefaultResampler';
regParam.Optimizer = 'RegularStepGradientDescent';
regParam.Transform = 'SimilarityTransform';
regParam.Metric = 'AdvancedMattesMutualInformation';
regParam.NumberOfHistogramBins = 32;
% Scales the rotations compared to the translations, to make
% sure they are in the same range. The higher this parameter,
% the smaller the changes in rotation angle in each iteration.
% If you have the feeling that rotations are not found by elastix,
% decrease it; if elastix crashes after a few iterations, with
% the message that all samples map outside the moving image 
% buffer, you may have to increase this parameter.
regParam.Scales = 50000.0;
regParam.AutomaticTransformInitialization = 'true';
regParam.UseDirectionCosines = 'true';
%param.CenterOfRotation = blockc;
regParam.NumberOfResolutions = 1;
regParam.ResultImagePixelType = 'unsigned char';
regParam.HowToCombineTransforms = 'Compose';
regParam.NumberOfSpatialSamples = round(info.Width * info.Height * .3);
regParam.ImageSampler = 'Random';
regParam.BSplineInterpolationOrder = 3;
regParam.FinalBSplineInterpolationOrder = 3;
regParam.DefaultPixelValue = 0;
regParam.MaximumNumberOfIterations = 500;
% Maximum step size of the RSGD optimizer for each resolution level.
% The higher you set this, the more aggressive steps are taken.
regParam.MaximumStepLength = 1.0;
% Minimum step size of the RSGD optimizer for each resolution level.
% The lower you set this, the more accurate the final result
regParam.MinimumStepLength = 0.005;
% Minimum magnitude of the gradient (stopping criterion) for the RSGD optimizer:
% The lower you set this, the more accurate the final result may be
regParam.MinimumGradientMagnitude = 0.000001;
regParam.ResultImageFormat = 'png';
regParam.CompressResultImage = 'true';

end
