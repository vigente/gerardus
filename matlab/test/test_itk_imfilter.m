% test_itk_imfilter.m
%
% Script to test the filters provided by itk_imfilter.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
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

%% itk::AnisotropicDiffusionVesselEnhancementImageFilter (advess) filter

% load test data
nrrd = scinrrd_load('../../cpp/src/third-party/IJ-Vessel_Enhancement_Diffusion.1/CroppedWholeLungCTScan.mhd');
nrrd.data = single(nrrd.data);

% plot data
close all
imagesc(nrrd.data(:, :, 3))

% user-provided parameters
sigmaMin = nrrd.axis(1).spacing * 5;
sigmaMax = nrrd.axis(1).spacing * 10;
sigmaSteps = 15;
isSigmaStepLog = false;
iterations = 30;
wStrength = 24;
sensitivity = 4.0;
timeStep = 1e-3;
epsilon = 1e-2;

% smooth along vessels (11.4 s)
tic
im = itk_imfilter('advess', nrrd, sigmaMin, sigmaMax, sigmaSteps, ...
    isSigmaStepLog, iterations, wStrength, sensitivity, timeStep, epsilon);
toc

figure
imagesc(im(:, :, 3))

figure
imagesc(nrrd.data(:, 2:end, 3) - im(:, 2:end, 3))

%% itk::MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter (hesves)

% load test data
nrrd = scinrrd_load('../../cpp/src/third-party/IJ-Vessel_Enhancement_Diffusion.1/CroppedWholeLungCTScan.mhd');

% plot data
close all
imagesc(nrrd.data(:, :, 5))

% user-provided parameters
sigmaMin = 1;
sigmaMax = 8;
sigmaSteps = 4;
isSigmaStepLog = false;

% compute vesselness measure
tic
im = itk_imfilter('hesves', nrrd.data, sigmaMin, sigmaMax, sigmaSteps, ...
    isSigmaStepLog);
toc

figure
imagesc(im(:, :, 5))

figure
imagesc(im(:, :, 5) .* double(nrrd.data(:, :, 5)));
