% toy_examples.m
%
% Obsolete, not used in final paper.
%
% Version: 0.1.6
% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2016 University of Oxford
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

BASEDIR = '/data2';
% BASEDIR = '/media/rcasero/mvldata/data';

DOCDIR = '/home/orie1416/Documents/20150112_3d_histology_diffusion_registration';
% DOCDIR = '/home/rcasero/Documents/20150112_3d_histology_diffusion_registration';

FIGDIR = [DOCDIR filesep 'figures'];

MOUSE = 'Q53';

IMPROCDIR = [BASEDIR '/Mouse/' MOUSE '/ImageProcessing'];
BFDATADIR = [BASEDIR '/Mouse/' MOUSE '/Blockface'];
BFSDATADIR = [BASEDIR '/Mouse/' MOUSE '/BlockfaceStabilised'];
BFCDATADIR = [BASEDIR '/Mouse/' MOUSE '/BlockfaceCorrected'];
% HISTO = 'HistologyTrichrome';
HISTO = 'HistologySiriusRed';
HISTOLO = [HISTO 'LoRes'];
HISTLODIR = [BASEDIR '/Mouse/' MOUSE filesep HISTOLO];
% HISTDIR = [BASEDIR '/Mouse/' MOUSE filesep HISTO];
% HISTLOREGDIR = [BASEDIR '/Mouse/' MOUSE filesep HISTOLO 'ToBlockfaceCorrectedRegistration'];
% HISTLO2BFREGSIMILARITYDIR = [IMPROCDIR filesep HISTOLO 'GrayToBlockfaceRegistrationSimilarity'];
% HISTLO2BFREGBSPLINEDIR = [IMPROCDIR filesep HISTOLO 'GrayToBlockfaceRegistrationBspline'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% vertical translation using transfdiff()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of slices
N = 10;

% random initialization of slice positions away from their perfect
% alignment
rng(5)
t = 3 * (rand(N, 1)-0.5) + 5;

% plot slice positions
hold off
plot(t, 'o', 'LineWidth', 4)
set(gca, 'FontSize', 18)
hold on
plot(t, 'LineWidth', 4)

% apply registration diffusion
opt.transform = 'TranslationTransform';
opt.MaxIter = Inf;
opt.Alpha = 0.49;
opt.Epsilon = 1e-4;
[ftot, info] = transfdiff(opt, t(2:end)-t(1:end-1));

% plot result
hold on
plot(t + ftot, 'r', 'LineWidth', 4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% registration diffusion with affine transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coefficient
alpha = 0.5;

% number of slices
N = 10;

% random initialization of slice positions away from their perfect
% alignment
rng(5)
for I = 1:N

    while(1)
        
        % generate a random affine transform
        aux = rand(3, 2);
        aux(1:2, 3, :) = 0;
        aux(3, 3, :) = 1;
        
        % check for invertibility
        if (abs(det(aux)) < 1e-3)
            % ill-conditioned matrix, so we generate another random
            % transform
            continue
        end
        
        % check for positive eigenvectors
        eigd = eig(aux);
        if (any(eig(aux) <= 0))
            % negative eigenvalues, for which logm is not defined, so
            % generate another random transform
            continue
        end
        
        % the transform is valid, so we can add it to the data set
        t(:, :, I) = aux;
        break
    
    end
    
end

% % plot slice positions
% hold off
% plot(t, 'o', 'LineWidth', 4)
% set(gca, 'FontSize', 18)
% hold on
% plot(t, 'LineWidth', 4)

% compute registration between slices
for I = 2:N
    fprev(:, :, I) = t(:, :, I) \ t(:, :, I-1);
end
fprev(:, :, 1) = [];
for I = 1:N-1
    fpos(:, :, I) = t(:, :, I) \ t(:, :, I+1);
end

% init accumulated transform to apply to each slice
ftot = [1 0 0; 0 1 0; 0 0 1];
ftot = repmat(ftot, 1, 1, N);

% apply registration diffusion
for iter = 1:400

    % compute the transform to apply to each slice in this iteration
    f = real(fpos(:, :, 1)^alpha);
    for I = 1:N-2
        f(:, :, I+1) = ...
            real(expm(alpha * (logm(fprev(:, :, I)) + logm(fpos(:, :, I+1)))));
    end
    f(:, :, N) = real(fprev(:, :, N-1)^alpha);
    
    % compose with the previous transforms to obtain a total transform
    for I = 1:N
        ftot(:, :, I) = ftot(:, :, I) * f(:, :, I);
    end
    
    % update the transforms between each slice and its neighbours
    for I = 2:N
        % recompute the registration of each slice to its neighbours
        fprev(:, :, I-1) = f(:, :, I) \ fprev(:, :, I-1) * f(:, :, I-1);
    end
    for I = 1:N-1
        % recompute the registration of each slice to its neighbours
        fpos(:, :, I) = f(:, :, I) \ fpos(:, :, I) * f(:, :, I+1);
    end
    
end

% check whether the slices are aligned after the correction
tcorr = t;
for I = 1:N
    tcorr(:, :, I) = real(t(:, :, I) * ftot(:, :, I));
end
err = squeeze(sqrt(sum(sum(diff(tcorr, 1, 3).^2, 1), 2)));
max(err)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% affine transform using transfdiff()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of slices
N = 10;

% random initialization of slice positions away from their perfect
% alignment
rng(5)
t = zeros(3, 3, N);
for I = 1:N

    while(1)
        
        % generate a random affine transform
        aux = rand(3, 2);
        aux(1:2, 3, :) = 0;
        aux(3, 3, :) = 1;
        
        % check for invertibility
        if (abs(det(aux)) < 1e-3)
            % ill-conditioned matrix, so we generate another random
            % transform
            continue
        end
        
        % check for positive eigenvectors
        eigd = eig(aux);
        if (any(eig(aux) <= 0))
            % negative eigenvalues, for which logm is not defined, so
            % generate another random transform
            continue
        end
        
        % the transform is valid, so we can add it to the data set
        t(:, :, I) = aux;
        break
    
    end
    
end

% compute registration between slices
fpos = zeros(3, 3, N-1);
for I = 1:N-1
    fpos(:, :, I) = t(:, :, I) \ t(:, :, I+1);
end

% apply registration diffusion
opt.transform = 'AffineTransform';
opt.MaxIter = Inf;
opt.Alpha = 0.49;
opt.Epsilon = 1e-4;
[ftot, info] = transfdiff(opt, fpos);

% check whether the slices are aligned after the correction
tcorr = t;
for I = 1:N
    tcorr(:, :, I) = real(t(:, :, I) * ftot(:, :, I));
end
err = squeeze(sqrt(sum(sum(diff(tcorr, 1, 3).^2, 1), 2)));
max(err)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% format of B-spline transform
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % enable access to elastix shared object
% libpath = getenv('LD_LIBRARY_PATH');
% setenv('LD_LIBRARY_PATH', [ '/usr/local/bin/elastix/lib:' libpath]);
% 
% BFCDATADIR = ['/data2/Mouse/' MOUSE '/BlockfaceCorrected'];
% 
% % load image
% im0 = imread([BFCDATADIR filesep 'Q53_55_0150.png']);
% 
% % register image to itself with a B-spline
% opts.verbose = 1;
% opts.t0 = '';
% [t, imout, iterinfo] = elastix(...
%     'basic_bspline_params.txt', ...
%     im0, im0, opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% common code to load one example histology slice. We are then going to 
%% repeat and deform this slice to create toy examples to illustrate
%% transform diffusion registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% enable access to elastix shared object
libpath = getenv('LD_LIBRARY_PATH');
setenv('LD_LIBRARY_PATH', [ '/usr/local/bin/elastix/lib:' libpath]);

HISTO = 'HistologyTrichrome';
HISTOLO = [HISTO 'LoRes'];
HISTLODIR = [BASEDIR '/Mouse/' MOUSE filesep HISTOLO];

% load image
im0 = imread([HISTLODIR filesep 'Q53_151 - 2014-07-01 17.16.41.png']);

% reduce size of image to make things faster
im0 = imresize(im0, 0.5);

% make grayscale copy
im0g = rgb2gray(im0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1 slice repeated 5 times, each with a random B-spline deformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optTransformix.AutoDefaultPixelValue = true;

% load registration parameters
regParamBsp = elastix_read_file2param([DOCDIR '/matlab/regParam-bspline-gray-intrahisto.txt']);

% image size has been reduced to make things faster
regParamBsp.FinalGridSpacingInVoxels = 200 * 0.5;

% make default pixel value some value that is similar to the background
regParamBsp.DefaultPixelValue = median(im0(:));

% register image to itself with a B-spline
opts.verbose = 0;
opts.t0 = '';
[t0, ~, iterinfo] = elastix(regParamBsp, ...
    im0g, im0g, opts);

% put control points in grid form
[gx0, gy0] = elastix_bspline_grid(t0);

% plot grid
hold off
plot(gx0(:), gy0(:), 'o')
box = [
    t0.Origin
    t0.Origin(1)+(t0.Size(1)-1)*t0.Spacing(1), t0.Origin(2)
    t0.Origin+(t0.Size-1).*t0.Spacing
    t0.Origin(1), t0.Origin(2)+(t0.Size(2)-1)*t0.Spacing(2)
    t0.Origin
    ];
hold on
plot(box(:, 1), box(:, 2), 'r', 'LineWidth', 2)

% now that we have a template for the B-spline transform, we are going to
% generate a few randomly deformed transforms of the base image
N = 5;
im = zeros([size(im0) N], 'like', im0);
clear t
rng(5)
for I = 1:N

    % we only move the control points within the image domain (I think this
    % is not exactly right, because a cubic B-spline has 2 extra control
    % points on one side, and 1 on the other)
    gx = gx0;
    gx(3:end-2, 3:end-2) = gx(3:end-2, 3:end-2) ...
        + (rand(size(gx)-4) - 0.5) * t0.GridSpacing(1)/2;
    gy = gy0;
    gy(3:end-2, 3:end-2) = gy(3:end-2, 3:end-2) ...
        + (rand(size(gy)-4) - 0.5) * t0.GridSpacing(2)/2;
    
%     % plot grid with perturbation
%     hold on
%     plot(gx(:), gy(:), 'x')
    
    % convert grid back to vector format
    t(I) = t0;
    t(I).TransformParameters = elastix_bspline_grid2param(gx-gx0, gy-gy0);
    
    % transform image
    im(:, :, :, I) = transformix(t(I), im0, optTransformix);

end

% plot images
for I = 1:N
   
    % plot image
    hold off
    imagesc(im(:, :, :, I))
    axis ij
    pause

end

% preprocess all images for registration
im2 = zeros(size(im), 'like', im);
mask2 = uint8(im(:, :, 1, 1));
for I = 1:N
    [~, im2(:, :, :, I), mask2(:, :, I)] ...
        = histology_preprocessing(im0, im(:, :, :, I));
end
im2g = zeros(size(im, 1), size(im, 2), size(im, 4), 'like', im);
for I = 1:N
    im2g(:, :, I) = rgb2gray(im2(:, :, :, I));
end

% plot images
for I = 1:N
   
    % plot image
    subplot(2, 1, 1)
    hold off
    imagesc(im2(:, :, :, I))
    axis ij

    % plot mask
    subplot(2, 1, 2)
    hold off
    imagesc(mask2(:, :, I))
    axis ij
    pause

end


% % read registration parameters
% regParamBsp = elastix_read_file2param([ DOCDIR '/matlab/regParam-bspline-intrahisto.txt']);
% regParamBsp.NumberOfResolutions = 3;
% regParamBsp.ImagePyramidSchedule = [16 16 4 4 1 1];
% regParamBsp.FinalGridSpacingInVoxels = 25;
% regParamBsp.GridSpacingSchedule = [16 4 1];
% regParamBsp.MaximumStepLength = regParamBsp.MaximumStepLength * 0.5;

% apply registration diffusion algorithm
clear optReg optDiff
optReg.verbose = true;
optReg.SpatialSamplesRatio = 0.1;
optReg.MaxIter = 16;
optReg.MaxVal = 1;
optReg.mask = mask2;
optReg.RegParam = [DOCDIR '/matlab/regParam-bspline-gray-intrahisto.txt'];
optDiff.MaxIter = 5;

[tout, ~, ~, imout] = transfdiffreg(...
    'BSplineTransform', im2g, optReg, optDiff);

clear tMax tMed
for L = 1:elastix_length(tout)
    taux = elastix_colon(tout, L);
    tMat = cat(1, taux(:).TransformParameters);
    tMed(:, L) = median(abs(tMat), 2);
    tMax(:, L) = max(abs(tMat), [], 2);
end

% apply transform to images
taux = elastix_colon(tout, 1:elastix_length(tout));
imout = im;
for I = 1:N
    imout(:, :, :, I) = transformix(tout(I), im(:, :, :, I), optTransformix);
end

% plot images
for I = 1:N
   
    % plot image
    subplot(2, 1, 1)
    hold off
    imagesc(imout(:, :, :, I))
    axis ij
    title(['image ' num2str(I)])
    pause

end

% plot difference between images and reference
for I = 1:N
   
    % plot image
    subplot(2, 1, 1)
    hold off
    aux = imfuse(im(:, :, :, I), im0);
    imagesc(aux)
    axis ij

    % plot image
    subplot(2, 1, 2)
    hold off
    aux = imfuse(imout(:, :, :, I), im0);
    imagesc(aux)
    axis ij
    pause

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Deterministic rigid deformation correction toy example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% reference image with the rigid deformation, but without B-spline noise.
%% This should only be run when we want to generate a ground truth image.
%% Otherwise, skip this section, and run the 2 next ones to generate the
%% B-spline noise and then the deterministic rigid deformation

% load rigid registration parameters
regParamRig = elastix_read_file2param([DOCDIR filesep 'matlab' ...
    filesep 'regParam-rigid.txt']);

% make default pixel value some value that is similar to the background
regParamRig.DefaultPixelValue = median(im0(:));

% register image to itself to get a rigid transform template
opts.verbose = 0;
opts.t0 = '';
[trig, ~, iterinfo] = elastix(regParamRig, ...
    im0, im0, opts);

% change the rotation center to the middle of the left ventricle
trig.CenterOfRotationPoint = [364 284];

% TransformParameters = [rotation, tx, ty]

% introduce a "Pisa Tower" effect, plus a rotation
tx = linspace(-30, 100, N);
rot = linspace(0, pi/2, N);
for I = 1:N
   
    % add a slowly varying rigid deformation
    trig.TransformParameters = [rot(I) -50 tx(I)];
    
    im(:, :, :, I) = transformix(trig, im0);
    
end

% plot cross horizontal section
subplot(2, 1, 1)
imagesc(permute(squeeze(im(300, :, :, :)), [3 1 2]))

% plot cross vertical section
subplot(2, 1, 2)
imagesc(permute(squeeze(im(:, 350, :, :)), [3 1 2]))

% save images
saveas(gca, ...
    [DOCDIR '/figures/results-bspline-correction-data-deterministic-rigid.png'])


%% B-spline noisy deformations

% load registration parameters
regParamBsp = elastix_read_file2param([DOCDIR '/matlab/regParam-bspline-intrahisto.txt']);

% image size has been reduced to make things faster
regParamBsp.FinalGridSpacingInVoxels = regParamBsp.FinalGridSpacingInVoxels * 0.5;

% make default pixel value some value that is similar to the background
regParamBsp.DefaultPixelValue = median(im0(:));

% register image to itself with a B-spline
opts.verbose = 0;
opts.t0 = '';
[tbsp, ~, iterinfo] = elastix(regParamBsp, ...
    im0, im0, opts);

% put control points in grid form
[gx0, gy0] = elastix_bspline_grid(tbsp);

% now that we have a template for the B-spline transform, we are going to
% generate a few randomly deformed transforms of the base image
N = 250;
im = zeros([size(im0) N], 'like', im0);
clear t
rng(5)
for I = 1:N

    % we only move the control points within the image domain (I think this
    % is not exactly right, because a cubic B-spline has 2 extra control
    % points on one side, and 1 on the other)
    gx = gx0;
    gx(3:end-2, 3:end-2) = gx(3:end-2, 3:end-2) ...
        + (rand(size(gx)-4) - 0.5) * tbsp.GridSpacing(1)/4;
    gy = gy0;
    gy(3:end-2, 3:end-2) = gy(3:end-2, 3:end-2) ...
        + (rand(size(gy)-4) - 0.5) * tbsp.GridSpacing(2)/4;
    
%     % plot grid with perturbation
%     hold on
%     plot(gx(:), gy(:), 'x')
    
    % convert grid back to vector format
    t(I) = tbsp;
    t(I).TransformParameters = elastix_bspline_grid2param(gx-gx0, gy-gy0);
    
    % transform image
    im(:, :, :, I) = transformix(t(I), im0);

end

%% deterministic rigid deformations

% load rigid registration parameters
regParamRig = elastix_read_file2param([DOCDIR filesep 'matlab' ...
    filesep 'regParam-rigid.txt']);

% make default pixel value some value that is similar to the background
regParamRig.DefaultPixelValue = median(im0(:));

% register image to itself to get a rigid transform template
opts.verbose = 0;
opts.t0 = '';
[trig, ~, iterinfo] = elastix(regParamRig, ...
    im0, im0, opts);

% change the rotation center to the middle of the left ventricle
trig.CenterOfRotationPoint = [364 284];

% TransformParameters = [rotation, tx, ty]

% introduce a "Pisa Tower" effect, plus a rotation
tx = linspace(-30, 100, N);
rot = linspace(0, pi/2, N);
for I = 1:N
   
    % add a slowly varying rigid deformation
    trig.TransformParameters = [rot(I) -50 tx(I)];
    
    im(:, :, :, I) = transformix(trig, im(:, :, :, I));
    
end

% plot images
for I = 1:N
   
    % plot image
    hold off
    imagesc(im(:, :, :, I))
    axis ij
    pause

end

% plot cross horizontal section
subplot(2, 1, 1)
imagesc(permute(squeeze(im(300, :, :, :)), [3 1 2]))

% plot cross vertical section
subplot(2, 1, 2)
imagesc(permute(squeeze(im(:, 350, :, :)), [3 1 2]))

% save images
saveas(gca, ...
    [DOCDIR '/figures/results-bspline-correction-data-noisy-bspline-deterministic-rigid.png'])

% preprocess all images for registration
im2 = im;
mask2 = uint8(im(:, :, 1, 1));
for I = 1:N
    [~, im2(:, :, :, I), mask2(:, :, I)] = histology_preprocessing(im0, im(:, :, :, I));
end

% plot cross horizontal section
subplot(2, 1, 1)
imagesc(permute(squeeze(im2(300, :, :, :)), [3 1 2]))

% plot cross vertical section
subplot(2, 1, 2)
imagesc(permute(squeeze(im2(:, 350, :, :)), [3 1 2]))

%% correction with only B-spline deformations, and not letting much
%% diffusion, so that we don't lose the deterministic rigid deformation

% read registration parameters
regParamBsp = elastix_read_file2param([DOCDIR '/matlab/regParam-bspline-intrahisto.txt']);
regParamBsp.NumberOfResolutions = 3;
regParamBsp.ImagePyramidSchedule = [16 16 4 4 1 1];
regParamBsp.FinalGridSpacingInVoxels = 25;
regParamBsp.GridSpacingSchedule = [16 4 1];
regParamBsp.MaximumStepLength = regParamBsp.MaximumStepLength * 0.5;

% apply registration diffusion algorithm
clear opt optTransfdiff optElastix
optElastix.verbose = 0; 
opt.MaxIter = 4;
opt.MaxVal = 15;
opt.SpatialSamplesRatio = 0.1;
optTransfdiff = [];
optTransfdiff.MaxIter = 25;
[tout, info, infoTransfdiff, imout] = transfdiffreg(...
    regParamBsp, ...
    im2, mask2, opt, optTransfdiff, optElastix);

% apply transform to images
for I = 1:N
    tout(I).DefaultPixelValue = t(I).DefaultPixelValue;
    imout(:, :, :, I) = transformix(tout(I), im(:, :, :, I));
end

% plot cross horizontal section
subplot(2, 1, 1)
hold off
imagesc(permute(squeeze(imout(300, :, :, :)), [3 1 2]))

% plot cross vertical section
subplot(2, 1, 2)
hold off
imagesc(permute(squeeze(imout(:, 350, :, :)), [3 1 2]))

% plot images
for I = 1:N
   
    % plot original image
    subplot(2, 1, 1)
    hold off
    imagesc(im(:, :, :, I))
    axis ij
    
    % plot corrected image
    subplot(2, 1, 2)
    hold off
    imagesc(imout(:, :, :, I))
    axis ij
    pause

end

% save images
saveas(gca, ...
    [DOCDIR '/figures/results-bspline-correction-' ...
    num2str(optTransfdiff.MaxIter) '-diffusion-' ...
    num2str(opt.MaxIter) '-bspline-iter.png'])

% save results
save([DOCDIR '/matlab/results-bspline-correction-25-diffusion-4-bspline-iter.mat'], ...
    'tout', 'info', 'infoTransfdiff', 'regParamBsp', 'opt', ...
    'optTransfdiff', 'regParamRig', 'tx', 'rot')

%% correction with only B-spline deformations, full
%% diffusion, to lose the deterministic rigid deformation. We only need 1
%% B-spline iteration

% read registration parameters
regParamBsp = elastix_read_file2param([DOCDIR '/matlab/regParam-bspline-intrahisto.txt']);
regParamBsp.NumberOfResolutions = 3;
regParamBsp.ImagePyramidSchedule = [16 16 4 4 1 1];
regParamBsp.FinalGridSpacingInVoxels = 25;
regParamBsp.GridSpacingSchedule = [16 4 1];
regParamBsp.MaximumStepLength = regParamBsp.MaximumStepLength * 0.5;

% apply registration diffusion algorithm
clear opt optTransfdiff optElastix
optElastix.verbose = 0; 
opt.MaxIter = 1;
opt.MaxVal = 15;
opt.SpatialSamplesRatio = 0.1;
optTransfdiff = [];
[tout, info, infoTransfdiff, imout] = transfdiffreg(...
    regParamBsp, ...
    im2, mask2, opt, optTransfdiff, optElastix);

% apply transform to images
for I = 1:N
    tout(I).DefaultPixelValue = t(I).DefaultPixelValue;
    imout(:, :, :, I) = transformix(tout(I), im(:, :, :, I));
end

% plot cross horizontal section
subplot(2, 1, 1)
hold off
imagesc(permute(squeeze(imout(300, :, :, :)), [3 1 2]))

% plot cross vertical section
subplot(2, 1, 2)
hold off
imagesc(permute(squeeze(imout(:, 350, :, :)), [3 1 2]))

% save images
saveas(gca, ...
    [DOCDIR '/figures/results-bspline-correction-' ...
    '-max-diffusion-' ...
    num2str(opt.MaxIter) '-bspline-iter.png'])

% save results
save([DOCDIR '/matlab/results-bspline-correction-' ...
    'max-diffusion-' num2str(opt.MaxIter) '-bspline-iter.mat'], ...
    'tout', 'info', 'infoTransfdiff', 'regParamBsp', 'opt', ...
    'optTransfdiff', 'regParamRig', 'tx', 'rot')

