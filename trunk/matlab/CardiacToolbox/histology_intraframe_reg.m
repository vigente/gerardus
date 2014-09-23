function [tParam, regParam] = histology_intraframe_reg(pathstr, files, opts)
% histology_intraframe_reg  Registration diffusion of histology slices.
%
% histology_intraframe_reg computes the registration of each slices to the
% two adjacent slices. The two resulting transforms are combined to
% produce a compromise transform. The process is applied iteratively.
%
% [TPARAM, REGPARAM] = histology_intraframe_reg(PATHSTR, FILES)
%
%   PATHSTR is the full path to the files. If empty, PATHTOFILES='.'.
%
%   FILES is the result of a dir() command, e.g. dir('*.png'). The
%   function expects a list of histology images.
%
% ... = histology_intraframe_reg(..., OPTS)
%
%   OPTS is a struct with algorithm parameters:
%
%     'MaxIter': (def 10) Maximum number of diffusion iterations.
%
%     'DiffusionCoefficient': (def 0.4) Value of diffusion coefficient.
%
%     'verbose': (default 0) Show elastix output on the screen.
%
% See also: histology_preprocessing.

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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

DEBUG = 0;

% check arguments
narginchk(2, 3);
nargoutchk(0, 4);

if (length(files) < 3)
    error('At least three files are needed for registration')
end

% defaults
if (nargin < 3 || isempty(opts) || ~isfield(opts, 'MaxIter'))
    opts.MaxIter = 10;
end
if (isempty(opts) || ~isfield(opts, 'DiffusionCoefficient'))
    opts.DiffusionCoefficient = 0.4;
end
if (isempty(opts) || ~isfield(opts, 'verbose'))
    opts.verbose = 0;
end

% defaults
if (isempty(pathstr))
    pathstr = '.';
end

% create a struct with the registration parameters
regParam = generate_registration_parameters();

% create a duplicate of the images that we will deform repeteadly in the
% diffusion process
pathtemp = '/tmp';
for I = 1:length(files)
    copyfile([pathstr filesep files(I).name], ...
        [pathtemp filesep files(I).name]); 
end

% diffusion computation loop
for S = 1:opts.MaxIter
    
    % process extreme slices (these only have one other adjacent slice)
    im2 = imread([pathstr filesep files(1).name]);
    im3 = imread([pathstr filesep files(2).name]);
    [im2, im3] = histology_preprocessing(im2, im3);
    if (S > 1)
        im2 = transformix(tParam(1), im2);
        im3 = transformix(tParam(2), im3);
    end
    [t3, im2to3reg, iterinfo3] = elastix(regParam, im3, im2, opts);
    t3.TransformParameters ...
        = 2 * t3.TransformParameters * opts.DiffusionCoefficient;
    if (S == 1)
        tParam(1) = t3;
    else
        t3.InitialTransformParametersFileName = tParam(1);
        tParam(1) = t3;
    end
    
    N = length(files);
    im1 = imread([pathstr filesep files(end-1).name]);
    im2 = imread([pathstr filesep files(end).name]);
    [im2, im1] = histology_preprocessing(im2, im1);
    if (S > 1)
        im1 = transformix(tParam(N-1), im1);
        im2 = transformix(tParam(N), im2);
    end
    [t1, im2to1reg, iterinfo1] = elastix(regParam, im1, im2, opts);
    t1.TransformParameters ...
        = 2 * t1.TransformParameters * opts.DiffusionCoefficient;
    if (S == 1)
        tParam(N) = t1;
    else
        t1.InitialTransformParametersFileName = tParam(N);
        tParam(N) = t1;
    end
    
    % loop for all the other slices
    for I = 2:length(files)-1
        disp(['S = ' num2str(S) ', I = ' num2str(I)])
       tic
        % load 3 consecutive images
        im1 = imread([pathstr filesep files(I-1).name]);
        im2 = imread([pathstr filesep files(I).name]);
        im3 = imread([pathstr filesep files(I+1).name]);
        
        % image preprocessing
        [~, im1] = histology_preprocessing(im2, im1);
        [im2, im3] = histology_preprocessing(im2, im3);
        
        if (DEBUG)
            % plot slices
            subplot(3, 1, 1)
            imagesc(im1)
            title(files(I-1).name, 'interpreter', 'none')
            axis([885 1776 475 1208])
            subplot(3, 1, 2)
            imagesc(im2)
            title(files(I).name, 'interpreter', 'none')
            axis([885 1776 475 1208])
            subplot(3, 1, 3)
            imagesc(im3)
            title(files(I+1).name, 'interpreter', 'none')
            axis([885 1776 475 1208])
        end
        
        % initial transform
        if (S > 1)
            % we cannot pass an initial transform to the adjacent slices,
            % so we have to transform them. Note that because we apply the
            % cumulative transforms to the original images, there's not a
            % problem with resampling images at each diffusion iteration
            im1 = transformix(tParam(I-1), im1);
            im2 = transformix(tParam(I), im2);
            im3 = transformix(tParam(I+1), im3);
        end
        
        % intra-histology registration  with B-spline
        opts.verbose = 0;
        [t1, im2to1reg, iterinfo1] = elastix(regParam, im1, im2, opts);
        [t3, im2to3reg, iterinfo3] = elastix(regParam, im3, im2, opts);
        
        if (DEBUG)
            % plot slices
            subplot(3, 1, 1)
            imagesc(double(im1) - double(im2to1reg))
            title(files(I-1).name, 'interpreter', 'none')
            axis([885 1776 475 1208])
            subplot(3, 1, 3)
            imagesc(double(im3) - double(im2to3reg))
            title(files(I+1).name, 'interpreter', 'none')
            axis([885 1776 475 1208])
        end
        
        % weighted average of the two transforms
        t1.TransformParameters ...
            = (t1.TransformParameters + t3.TransformParameters) ...
            * opts.DiffusionCoefficient;

        % update cumulative transform
        if (S == 1)
            tParam(I) = t1;
        else
            t1.InitialTransformParametersFileName = tParam(I);
            tParam(I) = t1;
        end
    toc    
    end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate_registration_parameters()
%
%   Create a struct with the parameters necessary for registration
function regParam = generate_registration_parameters()

regParam.FixedInternalImagePixelType = 'float';
regParam.FixedImageDimension = 2;
regParam.MovingInternalImagePixelType = 'float';
regParam.MovingImageDimension = 2;
regParam.Registration = 'MultiResolutionRegistration';
regParam.FixedImagePyramid = 'FixedSmoothingImagePyramid';
regParam.MovingImagePyramid = 'MovingSmoothingImagePyramid';
regParam.Interpolator = 'BSplineInterpolator';
regParam.Metric = 'AdvancedMeanSquares';
regParam.Optimizer = 'RegularStepGradientDescent';
regParam.ResampleInterpolator = 'FinalBSplineInterpolator';
regParam.Resampler = 'DefaultResampler';
regParam.Transform = 'BSplineTransform';
regParam.AutomaticTransformInitialization = 'false';
regParam.AutomaticTransformInitializationMethod = 'GeometricalCenter';
regParam.ErodeMask = 'false';
regParam.NumberOfResolutions = 4;
regParam.FinalGridSpacingInVoxels = 16;
regParam.GridSpacingSchedule = [32 16 8 4];
regParam.HowToCombineTransforms = 'Compose';
regParam.MaximumStepLength = 50;
% regParam.MinimumStepLength = 0.0500;
regParam.MinimumStepLength = 0.500;
regParam.MinimumGradientMagnitude = 1.0000e-08;
regParam.NumberOfHistogramBins = 32;
regParam.FixedLimitRangeRatio = 0;
regParam.MovingLimitRangeRatio = 0;
regParam.ImageSampler = 'Full';
regParam.FixedImageBSplineInterpolationOrder = 1;
regParam.UseRandomSampleRegion = 'false';
regParam.BSplineInterpolationOrder = 1;
regParam.FinalBSplineInterpolationOrder = 3;
regParam.DefaultPixelValue = 0;
regParam.WriteResultImage = 'true';
regParam.WriteResultImageAfterEachResolution = 'false';
regParam.WriteTransformParametersEachIteration = 'false';
regParam.WriteTransformParametersEachResolution = 'true';
regParam.ResultImagePixelType = 'unsigned char';
regParam.ResultImageFormat = 'png';
regParam.CompressResultImage = 'true';

end
