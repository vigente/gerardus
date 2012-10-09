% test_itk_imfilter.m
%
% Script to test the filters provided by itk_imfilter.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.3.0
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

%% itk::BinaryThinningImageFilter3D (skel) filter

% load test data
nrrd = scinrrd_load('../../cpp/src/third-party/IJ-Vessel_Enhancement_Diffusion.1/CroppedWholeLungCTScan.mhd');
nrrd.data = nrrd.data > 150;

% plot data
close all
imagesc(nrrd.data(:, :, 3))

% skeletonise (0.04 s)
tic
im = itk_imfilter('skel', nrrd);
toc

figure
imagesc(im(:, :, 3))


%% itk::DanielssonDistanceMapImageFilter (dandist) filter

% load test data
nrrd = scinrrd_load('../../cpp/src/third-party/IJ-Vessel_Enhancement_Diffusion.1/CroppedWholeLungCTScan.mhd');
nrrd2 = nrrd;
nrrd2.data(:) = 0;
nrrd2.data(:, 20, :) = 1;
nrrd2.data(:, 45, :) = 2;

% plot data
close all
imagesc(nrrd2.data(:, :, 3))

% compute distance (0.11 s)
tic
[im, v, w] = itk_imfilter('dandist', nrrd2);
toc

% distance map
figure
imagesc(im(:, :, 3))

% voronoi diagram
figure
imagesc(v(:, :, 3))


%% itk::SignedMaurerDistanceMapImageFilter (maudist) filter

% load test data
nrrd = scinrrd_load('../../cpp/src/third-party/IJ-Vessel_Enhancement_Diffusion.1/CroppedWholeLungCTScan.mhd');
nrrd2 = nrrd;
nrrd2.data(:) = 0;
nrrd2.data(:, 20:23, :) = 1;
nrrd2.data = uint16(nrrd2.data);

% plot data
close all
imagesc(nrrd2.data(:, :, 3))

% compute distance (0.11 s)
tic
im = itk_imfilter('maudist', nrrd2.data);
toc

% distance map
figure
imagesc(im(:, :, 3))


%% itk:::BinaryDilateImageFilter (bwdilate) filter

% load test data
nrrd = scinrrd_load('../../cpp/src/third-party/IJ-Vessel_Enhancement_Diffusion.1/CroppedWholeLungCTScan.mhd');
nrrd.data = nrrd.data > 150;

% plot data
close all
imagesc(nrrd.data(:, :, 3))

% dilate (0.02 s)
tic
im = itk_imfilter('bwdilate', nrrd, 4);
toc

figure
imagesc(im(:, :, 3))


%% itk:::BinaryErodeImageFilter (bwerode) filter

% load test data
nrrd = scinrrd_load('../../cpp/src/third-party/IJ-Vessel_Enhancement_Diffusion.1/CroppedWholeLungCTScan.mhd');
nrrd.data = nrrd.data > 150;

% plot data
close all
imagesc(nrrd.data(:, :, 3))

% erode (0.02 s)
tic
im = itk_imfilter('bwerode', nrrd, 1);
toc

figure
imagesc(im(:, :, 3))


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
nrrd.data = single(nrrd.data);

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


%% itk::MedianImageFilter (median)

% load test data
nrrd = scinrrd_load('../../cpp/src/third-party/IJ-Vessel_Enhancement_Diffusion.1/CroppedWholeLungCTScan.mhd');

% default median filtering (no filtering)
im2 = itk_imfilter('median', nrrd);

err = im2 - nrrd.data;
% expected result = 0
any(err(:) ~= 0)

% median filtering only along rows, line of length 5+1+5=11
im2 = itk_imfilter('median', nrrd, [5 0 0]);

% plot difference
subplot(2, 1, 1)
imagesc(nrrd.data(:, :, 4))
subplot(2, 1, 2)
imagesc(im2(:, :, 4))

% median filtering only along columns, line of length 20+1+20=41
im2 = itk_imfilter('median', nrrd, [0, 20, 0]);

% plot difference
subplot(2, 1, 1)
imagesc(nrrd.data(:, :, 4))
subplot(2, 1, 2)
imagesc(im2(:, :, 4))

% compare speed and results with Matlab's implementation
tic
im2 = itk_imfilter('median', nrrd, [10, 10, 10]); % 5.9 sec
toc
tic
im3 = medfilt3(nrrd.data, [21, 21, 21]); % 36.5 sec
toc

% difference: expected result = 0
any(double(im2(:)) - double(im3(:)))

% plot both results
subplot(2, 1, 1)
imagesc(im2(:, :, 4))
subplot(2, 1, 2)
imagesc(im3(:, :, 4))

%% itk::MRFImageFilter (mrf)

% load test data
scimat = scinrrd_load('../../cpp/src/third-party/IJ-Vessel_Enhancement_Diffusion.1/CroppedWholeLungCTScan.mhd');

% plot image
hold off
subplot(2, 2, 1)
imagesc(scimat.data(:,:,4))

% initial segmentation
[thr, q, obj, seg] = gmthr_seg(double(scimat.data), 2);

% quality of the separation between classes (expected: q = 0.9020)
q

% plot segmented image using the Gaussian mixture model
subplot(2, 2, 2)
imagesc(seg(:,:,4))

% compute neighbourhood weights that decrease with Euclidean distance
[gr, gc, gs] = ndgrid(-3:3, -3:3, -2:2);
weights = 1./sqrt(gr.^2 + gc.^2 + gs.^2);
weights(isinf(weights)) = 0.0;

% segment the image using the Markov Random Filter algorithm
seg = itk_imfilter('mrf', scimat.data, obj.mu', weights);

% plot segmented image using the Gaussian mixture model
subplot(2, 2, 3)
imagesc(seg(:,:,4))

% repeat segmentation, now with significant smoothing
seg = itk_imfilter('mrf', scimat.data, obj.mu', weights, 2);

% plot segmented image using the Gaussian mixture model
subplot(2, 2, 4)
imagesc(seg(:,:,4))

