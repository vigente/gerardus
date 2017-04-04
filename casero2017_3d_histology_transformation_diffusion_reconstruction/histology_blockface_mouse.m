% histology_blockface_mouse.m
%
% Script to reconstruct Mouse Q53 heart's low-res and hi-res histology,
% using blockface as external reference

% version: 0.6.5
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

% % enable access to elastix shared object
% libpath = getenv('LD_LIBRARY_PATH');
% setenv('LD_LIBRARY_PATH', [ '/usr/local/bin/elastix/lib:' libpath]);

% mouse
MOUSE = 'Q53';

% laptop at home
TOPDIR = '/home/rcasero/Software/casero2017_3d_histology_transformation_diffusion_reconstruction';
DATADIR = ['/media/rcasero/mvldata/data/Mouse/' MOUSE];

% workstation at the office
TOPDIR = '/home/orie1416/Software/casero2017_3d_histology_transformation_diffusion_reconstruction';
DATADIR = ['/data2/Mouse/' MOUSE];
% common directories for documentation
DOCDIR = [TOPDIR '/doc'];
SRCDIR = [TOPDIR '/src'];
FIGDIR = [TOPDIR '/doc/figures'];
RESULT_DIR = [SRCDIR '/results'];

% directory to save processed images to
IMPROC_DIR = [DATADIR '/Image_Processing'];

% binary masks for blockface stacks
BF_MASKS = [IMPROC_DIR filesep 'Blockface_Masks.mat'];

% blockface external reference (jumps corrected, scratches corrected,
% perspective corrected, illumination corrected)
BFP_DIR = [IMPROC_DIR '/Blockface_Preprocessed'];

% histology directories for data
% HISTO = 'HistologyTrichrome';
HISTO = 'HistologySiriusRed';
HISTOLO = [HISTO 'LoRes'];
HISTO_DIR = [DATADIR filesep HISTO];
HISTOLO_DIR = [DATADIR filesep HISTOLO];

% histology pre-aligned to blockface
HISTOLO_BF_DIR = [IMPROC_DIR filesep HISTOLO '_To_Blockface'];

% pre-processing of the pre-aligned histology so that it can be internally
% registered
HISTOLO_BF_P_DIR = [HISTOLO_BF_DIR '_Preprocessed'];
HISTOLO_BF_P2_DIR = [HISTOLO_BF_DIR '_Bsp_Preprocessed'];

% pre-processing of full resolution histology
HISTO_INIT_DIR = [IMPROC_DIR filesep HISTO '_Initial'];
HISTO_INIT_P_DIR = [HISTO_INIT_DIR '_Preprocessed'];

% high resolution histology refined with B-splines
HISTO_BSP_DIR = [IMPROC_DIR filesep HISTO '_Bsp_Diff'];

% transformations
T_HISTOLO_BF = [RESULT_DIR filesep 't_' HISTOLO '_to_Blockface.mat'];
T_HISTO_BF = [RESULT_DIR filesep 't_' HISTO '_to_Blockface.mat'];

% landmark errors from script study_stopping_criterion_rigid_refinement.m
ERR = [RESULT_DIR filesep 'err_' MOUSE '_' HISTOLO '.mat'];
ERRHI = [RESULT_DIR filesep 'err_' MOUSE '_' HISTO '.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visual assessment of histology slices that are too bad to work with
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Histology files that where acquired with big white bands, so that not all
% tissue is visible
filehBadWhiteBands = {
    'Q53-344.tif'
    'Q53-372.tif'
    'Q53-440.tif'
    'Q53-464.tif'
    'Q53-476.tif'
    'Q53-494.tif'
    'Q53-508.tif'
    'Q53-522.tif'
    'Q53-530.tif'
    'Q53-534.tif'
    'Q53-536.tif'
    'Q53-538.tif'
    'Q53-542.tif'
    'Q53-548.tif'
    'Q53-554.tif'
    'Q53-558.tif'
    'Q53-560.tif'
    'Q53-564.tif'
    'Q53-574.tif'
    'Q53-590.tif'
    'Q53-594.tif'
    'Q53-598.tif'
    };

filehBadWhiteBandsIdx = [
    344 372 440 464 476 494 508 522 530 534 536 538 542 548 554 558 560 ...
    564 574 590 594 598];

% histology slices that have too much broken tissue
filehBadTissue = {
    'Q53_120 - 2014-06-24 18.34.33.tif'
    'Q53_144 - 2014-06-24 19.02.10.tif'
    'Q53-316.tif'
    'Q53-360.tif'
    'Q53-396.tif'
    'Q53-398.tif'
    'Q53-400.tif'
    'Q53-402.tif'
    'Q53-404.tif'
    'Q53-406.tif'
    'Q53-412.tif'
    'Q53-414.tif'
    'Q53-416.tif'
    'Q53-418.tif'
    'Q53-422.tif'
    'Q53-424.tif'
    'Q53-426.tif'
    'Q53-428.tif'
    'Q53-430.tif'
    'Q53-432.tif'
    'Q53-434.tif'
    'Q53-456.tif'
    'Q53-460.tif'
    'Q53-462.tif'
    'Q53-466.tif'
    'Q53-474.tif'
    'Q53-480.tif'
    'Q53-482.tif'
    'Q53-484.tif'
    'Q53-486.tif'
    'Q53-488.tif'
    'Q53-502.tif'
    };

filehBadTissueIdx = [
    120 144 316 360 396 398 400 402 404 406 412 414 416 418 422 424 426 ...
    428 430 432 434 456 460 462 466 474 480 482 484 486 488 502];

filehBadCombined = union(filehBadWhiteBands, filehBadTissue);
filehBadCombinedIdx = union(filehBadWhiteBandsIdx, filehBadTissueIdx);

save(T_HISTOLO_BF, 'filehBadCombined', 'filehBadCombinedIdx', '-append')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rigid registration of histology to blockface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% list files

% list of matching blockface and histology files
[filebf, fileh] = match_file_lists(...
    dir([BFP_DIR filesep '*.mha']), ...
    dir([HISTOLO_DIR filesep '*.tif']), ...
    9:11, ...
    5:7 ...
    );

% indices of all existing pairs of blockface/histology slices
fileAllIdx = ...
    cell2mat(cellfun(@(x) str2num(x(9:11)), {filebf.name}, 'UniformOutput', false));

% remove bad slices
idx = ismember(fileAllIdx, filehBadCombinedIdx);
filebf(idx) = [];
fileh(idx) = [];

% save the list of files we are going to align
save(T_HISTOLO_BF, 'filebf', 'fileh', '-append')

% make an explicit list of good blockface/histology slice pairs
fileGoodIdx = setdiff(fileAllIdx, filehBadCombinedIdx);
save(T_HISTOLO_BF, 'fileAllIdx', 'fileGoodIdx', '-append')

N = length(filebf);

%% register slices

% load blockface masks
load(BF_MASKS, 'ellipmask90Crop', 'polymask90Crop')

% iterate histology images registering them to blockface
delete(gcp)
mypool = parpool('local', 3);
parfor I = 1:N

    % compute the transforms but don't apply them to the images. We have to
    % correct the output size, pixel size and offset
    disp(['I = ' num2str(I) ', Histology = ' fileh(I).name])
    th2bf(I) = aux_histology2blockface(...
        [HISTOLO_DIR filesep fileh(I).name], ...
        [BFP_DIR filesep filebf(I).name], ...
        polymask90Crop, ellipmask90Crop);

end

% compute an average output size, origin and pixel size for all the transforms
mSize = round(mean(cat(1, th2bf.Size)));
mSpacing = mean(cat(1, th2bf.Spacing));
mOrigin = mean(cat(1, th2bf.Origin));

% use auto estimation of background colour in transformix
optTransformix.AutoDefaultPixelValue = true;

% apply the mean output size, origin and pixel size to all the transforms,
% so that they produce images that can be directly stacked and compared
for I = 1:N
    th2bf(I).Size = mSize;
    th2bf(I).Spacing = mSpacing;
    th2bf(I).Origin = mOrigin;
    
    % load original histology slice
    imh0 = scimat_load([HISTOLO_DIR filesep fileh(I).name]);
    
    % apply rigid transform to original histology
    imh2bf = transformix(th2bf(I), imh0, optTransformix);
    
    % save registered image
    [~, name] = fileparts(fileh(I).name);
    scimat_save([HISTOLO_BF_DIR filesep name '.mha'], imh2bf);
end

% save transforms
save(T_HISTOLO_BF, 'th2bf', '-append')

%% visually check histology to blockface alignment

% list of blockface and histology files
[filebf, fileh] = match_file_lists(...
    dir([BFP_DIR filesep '*.mha']), ...
    dir([HISTOLO_BF_DIR filesep '*.mha']), ...
    9:11, ...
    5:7 ...
    );

% remove bad slices from the stack
idx_filehBadTissue = ismember({fileh.name}, filehBadCombined);
filebf(idx_filehBadTissue) = [];
fileh(idx_filehBadTissue) = [];

N = length(filebf);

aux = scimat_load([HISTOLO_BF_DIR filesep fileh(1).name]);
imh = aux.data;
imh = zeros([size(imh) N], 'like', imh);

for I = 1:N
    
    % load blockface and registered histology
    imbf = scimat_load([BFP_DIR filesep filebf(I).name]);
    aux = scimat_load([HISTOLO_BF_DIR filesep fileh(I).name]);
    imh(:, :, :, :, :, I) = aux.data;
    
    % save copy of the image in .tif format
    [~, filename] = fileparts([HISTOLO_BF_DIR filesep fileh(I).name]);
    scimat_save([HISTOLO_BF_DIR filesep filename '.tif'], aux);
    
    % plot alignment
%     subplot(2, 2, 1)
%     imagesc(imbf.data)
%     subplot(2, 2, 2)
%     imagesc(squeeze(imh(:, :, :, :, :, I)))
%     subplot(2, 1, 2)
%     imagesc(imfuse(imbf.data, squeeze(imh(:, :, :, :, :, I))))
%     pause
    
end

% plot cross horizontal section
subplot(2, 1, 1)
hold off 
imagesc(permute(squeeze(imh(485, :, :, :, :, :)), [3 1 2]))
set(gca, 'FontSize', 16)

% plot cross vertical section
subplot(2, 1, 2)
hold off
imagesc(permute(squeeze(imh(:, 550, :, :, :, :)), [3 1 2]))
set(gca, 'FontSize', 16)

saveas(gca, [FIGDIR filesep HISTOLO '_to_Blockface_cross_sections.png'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Affine matching of the second half of the stack to the first half
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load last slice from the 1st half, and first slice from the 2nd half
imh1 = scimat_load([HISTOLO_BF_DIR filesep fileh(148).name]); % moving
imh2 = scimat_load([HISTOLO_BF_DIR filesep fileh(149).name]); % fixed
imbf1 = scimat_load([BFP_DIR filesep filebf(148).name]);
imbf2 = scimat_load([BFP_DIR filesep filebf(149).name]);

% see how well each histology overlaps its blockface
subplot(2, 1, 1)
imagesc(imfuse(squeeze(imh1.data), imbf1.data))
subplot(2, 1, 2)
imagesc(imfuse(squeeze(imh2.data), imbf2.data))

% open interface so that we can click corresponding points on slice(I)
% and slice(I-1)
[movingPoints, fixedPoints] = cpselect(...
    squeeze(imh1.data), ... % moving
    squeeze(imh2.data), ... % fixed
    'Wait', true);

% DEBUG: Transformation in Matlab format
% % create reference frame
% r1 = imref2d(size(imh1.data));
% 
% % compute a projective transformation to correct the affine deformation
% % between the first and second halves of the stack
% tformCorrect2ndHalf = fitgeotrans(movingPoints, fixedPoints, 'affine');
% 
% % apply the transformation
% moving_reg = imh1;
% [moving_reg.data, rfixed] = imwarp(imh1.data, r1, tformCorrect2ndHalf, ...
%     'OutputView', imref2d(size(imh2.data)));
%
% % check registration
% subplot(2, 1, 1)
% imagesc(imfuse(squeeze(imh1.data), squeeze(imh2.data)))
% subplot(2, 1, 2)
% imagesc(imfuse(squeeze(moving_reg.data), squeeze(imh2.data)))

% convert point indices to real world coordinates (note that movingPoints,
% fixedPoints are given in the order X, Y, but as pixel indices)
movingPointsXY = scimat_index2world(movingPoints(:, [2 1]), imh1);
fixedPointsXY = scimat_index2world(fixedPoints(:, [2 1]), imh2);

% deform first half of the stack halfway towards second half
for I = 1:148
    
    % compute image transformation to apply to the slices
    tCorrectHalvesMismatch(I) = elastix_fitgeotrans(movingPointsXY, fixedPointsXY, ...
        'affine', th2bf(I).Size, th2bf(I).Spacing, th2bf(I).Origin, 0.5);
    
    
end

% deform second half of the stack halfway towards first half
for I = 149:N
    
    % compute image transformation to apply to the slices
    tCorrectHalvesMismatch(I) = elastix_fitgeotrans(fixedPointsXY, movingPointsXY, ...
        'affine', th2bf(I).Size, th2bf(I).Spacing, th2bf(I).Origin, 0.5);
    
    
end

% save transform
save(T_HISTOLO_BF, 'tCorrectHalvesMismatch', '-append')

% DEBUG: apply the transformation
clear opts
opts.AutoDefaultPixelValue = true;
foo = transformix(tCorrectHalvesMismatch(148), imh1, opts);
bar = transformix(tCorrectHalvesMismatch(149), imh2, opts);

% DEBUG: check registration
subplot(2, 1, 1)
imagesc(imfuse(squeeze(imh1.data), squeeze(imh2.data)))
subplot(2, 1, 2)
imagesc(imfuse(squeeze(foo.data), squeeze(bar.data)))

%% Visually check the mismatch refinement to the histology

% list of blockface and histology files
[filebf, fileh] = match_file_lists(...
    dir([BFP_DIR filesep '*.mha']), ...
    dir([HISTOLO_DIR filesep '*.tif']), ...
    9:11, ...
    5:7 ...
    );

% remove bad slices from the stack
idx_filehBadTissue = ismember({fileh.name}, filehBadCombined);
filebf(idx_filehBadTissue) = [];
fileh(idx_filehBadTissue) = [];

N = length(filebf);

% compose transforms
load(T_HISTOLO_BF, 'th2bf', 'tCorrectHalvesMismatch')
ttot = elastix_cat(tCorrectHalvesMismatch, th2bf);

% use auto estimation of background colour in transformix
optTransformix.AutoDefaultPixelValue = true;

% apply the mean output size, origin and pixel size to all the transforms,
% so that they produce images that can be directly stacked and compared
imh = zeros(971, 1099, 1, 1, 3, N, 'uint8');
for I = 1:N

    disp(['I = ' num2str(I) ', Histology = ' fileh(I).name])
    
    % load original histology slice
    imh0 = scimat_load([HISTOLO_DIR filesep fileh(I).name]);
    
    % apply transform to original histology
    aux = transformix(ttot(I), imh0, optTransformix);
    imh(:, :, :, :, :, I) = aux.data;

end

% plot cross horizontal section
subplot(2, 1, 1)
hold off 
imagesc(permute(squeeze(imh(485, :, :, :, :, :)), [3 1 2]))
set(gca, 'FontSize', 16)

% plot cross vertical section
subplot(2, 1, 2)
hold off
imagesc(permute(squeeze(imh(:, 550, :, :, :, :)), [3 1 2]))
set(gca, 'FontSize', 16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Intra-histology refinement: rigid diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% pre-processing of histology
%% It's faster to pre-process than transform 

% list of blockface and histology files
load(T_HISTOLO_BF, 'filebf', 'fileh')

N = length(fileh);

% compose transforms
load(T_HISTOLO_BF, 'th2bf', 'tCorrectHalvesMismatch')
ttot = elastix_cat(tCorrectHalvesMismatch, th2bf);

% use auto estimation of background colour in transformix
optTransformix.AutoDefaultPixelValue = true;

% load original histo
for I = 1:N

    disp(['I = ' num2str(I) ', Histology = ' fileh(I).name])
    
    % load original histology slice
    imh0 = scimat_load([HISTOLO_DIR filesep fileh(I).name]);
    
    % apply transform to original histology
    imh(I) = transformix(ttot(I), imh0, optTransformix);

end

% average all slices to create a reference histogram to match the rest of
% the slices to
imhref.data = zeros(size(imh(1).data));
for I = 1:N
    imhref.data = imhref.data + double(imh(I).data);
end
imhref.data = uint8(imhref.data / N);

% preprocess all images for registration
opts.Grayscale = true;
[~, imh, mask] = histology_preprocessing(imhref, imh, opts);

% save to directory
for I = 1:N
    
    disp(['I = ' num2str(I) ', Histology = ' fileh(I).name])
    scimat_save([HISTOLO_BF_P_DIR filesep fileh(I).name], imh(I));
    scimat_save([HISTOLO_BF_P_DIR filesep 'mask-' fileh(I).name], mask(I));
    
end

%% intra-histology refinement (here we only compute the registrations, 
%% refinements for different number of diffusions will be computed in 
%% study_stopping_criterion_rigid_refinement.m)

% load all preprocessed histology and mask
for I = 1:N
    
    imh(I) = scimat_load([HISTOLO_BF_P_DIR filesep fileh(I).name]);
    mask(I) = scimat_load([HISTOLO_BF_P_DIR filesep 'mask-' fileh(I).name]);
    
end

% rigid transform diffusion registration of the histology volume
clear optReg optDiff
optReg_rig.Angle = (-10:.25:10) / 180 * pi;
optReg_rig.verbose = true;
optDiff_rig.MaxIter = 1;
tic
[th2h_rig, infoReg_rig, infoRigDiff_rig] = transfdiffreg(...
    'regmatchedfilt', imh, optReg_rig, optDiff_rig);
toc

save(T_HISTOLO_BF, 'th2h_rig', 'infoReg_rig', 'infoRigDiff_rig', ...
    'optReg_rig', 'optDiff_rig', '-append')

%% Go to script study_stopping_criterion_rigid_refinement.m to test
%% different number of stack sweeps

% then we come back here

%% intra-histology refinement (we keep only the optimal number of sweeps
%% from all the refinements computed in 
%% study_stopping_criterion_rigid_refinement.m)

load(ERR, 'maxIter', 't_rig', 'infoTransfdiff_rig', 'optDiff_rig')
load(T_HISTOLO_BF, 'infoReg_rig', 'infoRigDiff_rig', 'optReg_rig', ...
    'optDiff_rig')

I = 30;
disp(['maxIter = ' num2str(maxIter(I))])

% transfer diffused parameters to elastix structs
th2h_rig = elastix_affine_matrix2struct(t_rig{I}, infoReg_rig.tp(1));

% save transformation and associated variables to where we would save them
% if we computed them in this script, instead of loading them from
% study_stopping_criterion_rigid_refinement.m
save(T_HISTOLO_BF, 'th2h_rig', 'infoReg_rig', 'infoRigDiff_rig', ...
    'optReg_rig', 'optDiff_rig', '-append')

%% visually check histology to histology alignment

clear imh
opts.AutoDefaultPixelValue = true;

% compose transforms
load(T_HISTOLO_BF, 'th2bf', 'tCorrectHalvesMismatch', 'th2h_rig')
ttot = elastix_cat(th2h_rig, tCorrectHalvesMismatch, th2bf);

% init memory to store transformed histology
imh = zeros([ttot(1).Size([2 1]) 1 1 3 N], 'uint8');

for I = 1:N
    
    % load histology
    [~, nameh] = fileparts(fileh(I).name);
    aux = scimat_load([HISTOLO_DIR filesep nameh '.tif']);
    
    % apply transforms to histology
    aux = transformix(ttot(I), aux, opts);
    imh(:, :, :, :, :, I) = aux.data;
    
end

% coordinates of image limits
aux0 = scimat_index2world([1 1], aux);
auxend = scimat_index2world([size(aux.data, 1) size(aux.data, 2)], aux);
bbx = [aux0(1), auxend(1)]; % x-coordinates of the first and last corners of bounding box
bby = [aux0(2), auxend(2)]; % y-coordinates of the first and last corners of bounding box
bbz = 10e-6*[0 N-1]; % z-coordinates of eacn pixel

% plot cross horizontal section
subplot(2, 1, 1)
hold off 
imagesc(bbx * 1e3, bbz * 1e3, permute(squeeze(imh(485, :, :, :, :, :)), [3 1 2]))
set(gca, 'FontSize', 16)

% plot cross vertical section
subplot(2, 1, 2)
hold off
imagesc(bby * 1e3, bbz * 1e3, permute(squeeze(imh(:, 550, :, :, :, :)), [3 1 2]))
set(gca, 'FontSize', 16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Intra-histology registration refinement: B-spline diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear imh
opts.AutoDefaultPixelValue = true;

% compose transforms
load(T_HISTOLO_BF, 'th2bf', 'tCorrectHalvesMismatch', 'th2h_rig')
ttot = elastix_cat(th2h_rig, tCorrectHalvesMismatch, th2bf);

% init memory to store transformed histology
[~, nameh] = fileparts(fileh(1).name);
imh = scimat_load([HISTOLO_DIR filesep nameh '.tif']);
imh(1:N) = imh;

% apply previous transforms to original histology (this is the same we do
% above to visualize the whole histology)
for I = 1:N
    
    % load histology
    [~, nameh] = fileparts(fileh(I).name);
    imh(I) = scimat_load([HISTOLO_DIR filesep nameh '.tif']);
    
    % apply transforms to histology
    imh(I) = transformix(ttot(I), imh(I), opts);
    
end

% average all slices to create a reference histogram to match the rest of
% the slices to
imhref.data = zeros(size(imh(1).data));
for I = 1:N
    imhref.data = imhref.data + double(imh(I).data);
end
imhref.data = uint8(imhref.data / N);

% preprocess all images for registration
opts.Grayscale = true;
mask = imh(1);
mask.data = uint8(mask.data(:, :, 1, 1, 1));
mask(1:N) = mask;
for I = 1:N
    
    [~, imh(I), mask(I)] ...
        = histology_preprocessing(imhref, imh(I), opts);
    
end

% save to directory to save time when we are testing different parameters
% for the B-spline refinement
for I = 1:N
    
    disp(['I = ' num2str(I) ', Histology = ' fileh(I).name])
    scimat_save([HISTOLO_BF_P2_DIR filesep fileh(I).name], imh(I));
    scimat_save([HISTOLO_BF_P2_DIR filesep 'mask-' fileh(I).name], mask(I));
    
end

%% intra-histology refinement (here we only compute the registrations, 
%% refinements for different number of diffusions will be computed in 
%% study_stopping_criterion_bsp_refinement.m)

% load histology slices and masks
clear imh mask
for I = 1:N
    
    disp(['I = ' num2str(I) ', Histology = ' fileh(I).name])
    imh(I) = scimat_load([HISTOLO_BF_P2_DIR filesep fileh(I).name]);
    mask(I) = scimat_load([HISTOLO_BF_P2_DIR filesep 'mask-' fileh(I).name]);
    
end

% B-spline transform diffusion registration of the histology volume
clear optReg_bsp optDiff_bsp
optReg_bsp.verbose = true;
optReg_bsp.SpatialSamplesRatio = 0.1;
optReg_bsp.MaxIter = 1;
optReg_bsp.mask = mask;
optReg_bsp.t0 = [];
optReg_bsp.RegParam = [SRCDIR '/regParam-bspline-gray-intrahisto.txt'];
optReg_bsp.MaxVal = 0.0; % don't stop before reaching maximum number of iterations
optDiff_bsp.Alpha = 0.45;
optDiff_bsp.MaxIter = 1;
tic
[th2h_bsp, infoReg_bsp, infoDiff_bsp, imout] = transfdiffreg(...
    'BSplineTransform', imh, optReg_bsp, optDiff_bsp);
toc

% save result
save(T_HISTOLO_BF, 'th2h_bsp', 'infoReg_bsp', 'infoDiff_bsp', ...
    'optReg_bsp', 'optDiff_bsp', '-append')

%% visually check histology to histology alignment

clear imh
opts.AutoDefaultPixelValue = true;

% compose transforms
load(T_HISTOLO_BF, 'th2bf', 'tCorrectHalvesMismatch', 'th2h_rig', 'th2h_bsp')
ttot = elastix_cat(th2h_bsp, th2h_rig, tCorrectHalvesMismatch, th2bf);

% init memory to store transformed histology
imh = zeros([ttot(1).Size([2 1]) 1 1 3 N], 'uint8');

for I = 1:N
    
    % load histology
    [~, nameh] = fileparts(fileh(I).name);
    aux = scimat_load([HISTOLO_DIR filesep nameh '.tif']);
    
    % apply transforms to histology
    aux = transformix(ttot(I), aux, opts);
    imh(:, :, :, :, :, I) = aux.data;
    
end

% coordinates of image limits
aux0 = scimat_index2world([1 1], aux);
auxend = scimat_index2world([size(aux.data, 1) size(aux.data, 2)], aux);
bbx = [aux0(1), auxend(1)]; % x-coordinates of the first and last corners of bounding box
bby = [aux0(2), auxend(2)]; % y-coordinates of the first and last corners of bounding box
bbz = 10e-6*[0 N-1]; % z-coordinates of eacn pixel

% plot cross horizontal section
subplot(2, 1, 1)
hold off 
imagesc(bbx * 1e3, bbz * 1e3, permute(squeeze(imh(485, :, :, :, :, :)), [3 1 2]))
set(gca, 'FontSize', 16)

% plot cross vertical section
subplot(2, 1, 2)
hold off
imagesc(bby * 1e3, bbz * 1e3, permute(squeeze(imh(:, 550, :, :, :, :)), [3 1 2]))
set(gca, 'FontSize', 16)

saveas(gca, [FIGDIR filesep INTRAHISTLOBSP '-cross-sections.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply transforms to high resolution histology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Applying to full resolution, 7.3 min / slice.
% Total 239 slices -> 29.1 h
% 779M / slice (histology), 448M / slice (preprocessed), 448M / slice(mask)
% Total 239 slices -> 391 GB
%
% At full resolution, distance between cell centres ~ 34 pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load reconstruction for low resolution. Here the B-spline refinement is
% only 1 registration sweep. The diffusion was computed in script
% study_stopping_criterion_bsp_refinement.m
load(T_HISTOLO_BF, 'th2bf', 'tCorrectHalvesMismatch', 'th2h_rig', ...
    'th2h_bsp', ...
    'infoReg_rig', 'infoDiff_rig', 'optReg_rig', 'optDiff_rig', ...
    'infoReg_bsp', 'infoDiff_bsp', 'optReg_bsp', 'optDiff_bsp', ...
    'filehBadCombinedIdx', 'filebf', 'fileh')

N = length(fileh);

% load
load(ERR, 'maxIter_bsp', 't_bsp', 'infoTransfdiff_bsp')

% we select low resolution refinement for 80 diffusion iterations
I = 38;
maxIter_bsp(I)

% transfer diffused parameters to elastix structs
for J = 1:N
    th2h_bsp(J) = infoReg_bsp.tp(1);
    th2h_bsp(J).TransformParameters = t_bsp{I}(J, :);
end

% compose transforms (for the histogram, we don't need B-spline refinement,
% and this way it's faster)
ttot = elastix_cat(th2h_rig, tCorrectHalvesMismatch, th2bf);

% load original histo
% use auto estimation of background colour in transformix
optTransformix.AutoDefaultPixelValue = true;
for I = 1:N

    disp(['I = ' num2str(I) ', Histology = ' fileh(I).name])
    
    % load original histology slice
    imh0 = scimat_load([HISTOLO_DIR filesep fileh(I).name]);
    
    % apply transform to original histology
    imh(I) = transformix(ttot(I), imh0, optTransformix);

end

% average all slices to create a reference histogram to match the rest of
% the slices to
imhref.data = zeros(size(imh(1).data));
for I = 1:N
    imhref.data = imhref.data + double(imh(I).data);
end
imhref.data = uint8(imhref.data / N);

clear imh



% compose transforms
ttot_lo2hi = elastix_cat(th2h_bsp, th2h_rig, tCorrectHalvesMismatch, th2bf);

% cropping parameters
originCrop = [0.001547298366564   0.002231223753707];
sizeCrop = [230 298];

% load spacing information from first slice (1st and 2nd halves of the
% block have different spacing)
filein = [HISTO_DIR filesep fileh(1).name];
info = scimat_load(filein, 'HeaderOnly', true);
highResDownsample = 2;
outSpacing = [info.axis.spacing] * highResDownsample;
outSpacing = outSpacing([2 1]);
K = ttot_lo2hi(J).Spacing ./ outSpacing;
outSize = round((ttot_lo2hi(J).Size - sizeCrop) .* K);

% create a high-res alignment from the low-res alignment
for J = 1:N

    tic
    disp(['J = ' num2str(J)])
    
    % convert low res transform to high res
    ttot_lo2hi(J).Spacing = outSpacing;
    ttot_lo2hi(J).Size = outSize;
    ttot_lo2hi(J).CompressResultImage = 'false';
    ttot_lo2hi(J).Origin = originCrop;
    
    ttot_lo2hi(J).InitialTransformParametersFileName.Spacing = outSpacing;
    ttot_lo2hi(J).InitialTransformParametersFileName.Size = outSize;
    ttot_lo2hi(J).InitialTransformParametersFileName.CompressResultImage = 'false';
    ttot_lo2hi(J).InitialTransformParametersFileName.Origin = originCrop;
    
    ttot_lo2hi(J).InitialTransformParametersFileName.InitialTransformParametersFileName.Spacing = outSpacing;
    ttot_lo2hi(J).InitialTransformParametersFileName.InitialTransformParametersFileName.Size = outSize;
    ttot_lo2hi(J).InitialTransformParametersFileName.InitialTransformParametersFileName.CompressResultImage = 'false';
    ttot_lo2hi(J).InitialTransformParametersFileName.InitialTransformParametersFileName.Origin = originCrop;
    
    ttot_lo2hi(J).InitialTransformParametersFileName.InitialTransformParametersFileName.InitialTransformParametersFileName.Spacing = outSpacing;
    ttot_lo2hi(J).InitialTransformParametersFileName.InitialTransformParametersFileName.InitialTransformParametersFileName.Size = outSize;
    ttot_lo2hi(J).InitialTransformParametersFileName.InitialTransformParametersFileName.InitialTransformParametersFileName.CompressResultImage = 'false';
    ttot_lo2hi(J).InitialTransformParametersFileName.InitialTransformParametersFileName.InitialTransformParametersFileName.Origin = originCrop;
    
end

save(T_HISTO_BF, 'ttot_lo2hi', '-append')


% apply the low res alignment to high res histology
for J = 1:N

    tic
    disp(['J = ' num2str(J)])
    
    %% apply low res transformations to high res images
    
    % input file
    filein = [HISTO_DIR filesep fileh(J).name];
    
    % output filename   
    optTransformix.outfile = [HISTO_INIT_DIR filesep fileh(J).name];
    
    % apply low res reconstruction to high res histology
    transformix(ttot_lo2hi(J), filein, optTransformix);
    
    %% pre-process histology slices
    
    % load high res slice
    imh = scimat_load(optTransformix.outfile);
    
    % preprocess all images for registration
    opts.Grayscale = true;
    [~, imh, mask] = histology_preprocessing(imhref, imh, opts);
    
    % save pre-processed image and mask
    [~, file] = fileparts(fileh(J).name);
    disp(['I = ' num2str(J) ', Histology = ' fileh(J).name])
    scimat_save([HISTO_INIT_P_DIR filesep file '.mha'], imh);
    % we use the low resolution masks, we don't need high resolution ones.
    % Anyway, masks are rough outlines
%     scimat_save([HISTO_INIT_P_DIR filesep 'mask-' file '.mha'], mask);
    toc
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Virtual slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear file

% we have cropped the edges of the high resolution histology. Using real
% world coordinates, find out which row and column correspond to the low
% resolution virtual slices
imh = scimat_load([HISTOLO_BF_P2_DIR filesep fileh(1).name]); % low res
aux = scimat_index2world([485, 550], imh);
imh = scimat_load([HISTO_INIT_DIR filesep fileh(1).name]); % high res
round(scimat_world2index(aux, imh)) % [R=4289, C=6237]

% create list of histology files with paths
for I = 1:length(fileh)
    file{I} = [HISTO_INIT_DIR filesep fileh(I).name];
end

% high res image size
NR = size(imh.data, 1);
NC = size(imh.data, 2);

tic
im = scimat_crop(file, [4289 nan nan], [4289 nan nan], 20e-6); % 16.2 min
toc
tic
im2 = scimat_crop(file, [nan 6237 nan], [nan 6237 nan], 20e-6); % 16.7 min
toc

im.data = permute(im.data, [3 2 1 4 5]);
im.axis = im.axis([3 2]);
im.axis(1).spacing = im.axis(1).spacing * 1e3;
im.axis(2).spacing = im.axis(2).spacing * 1e3;
im.rotmat = eye(2);

scimat_imagesc(im)
set(gca, 'FontSize', 18)
xlabel('x (mm)')
ylabel('z (mm)')
axis equal tight xy
drawnow

fig = gcf;
set(fig, 'Color', 'white')

export_fig(gca, [FIGDIR filesep 'virtual_slice_' MOUSE '_histo_hires_bsp_iter_0_LA_R4289.png'])

im2.data = permute(im2.data, [3 1 2 4 5]);
im2.axis = im2.axis([3 1]);
im2.axis(1).spacing = im2.axis(1).spacing * 1e3;
im2.axis(2).spacing = im2.axis(2).spacing * 1e3;
im2.rotmat = eye(2);

scimat_imagesc(im2)
set(gca, 'FontSize', 18)
xlabel('y (mm)')
ylabel('z (mm)')
axis equal tight xy
drawnow

fig = gcf;
set(fig, 'Color', 'white')

export_fig(gca, [FIGDIR filesep 'virtual_slice_' MOUSE '_histo_hires_bsp_iter_0_SAX_C6237.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Refine reconstruction in high resolution  (registration only, as we try
%% different levels of diffusion in script study_stopping_criterion_bsp_refinement.m
%% based on this registration sweep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% intra-histology refinement

N = length(fileh);

% list of preprocessed histology slices and masks
for I = 1:N
    
    % list of preprocessed histology and mask files
    [~, name] = fileparts(fileh(I).name);
    disp(['I = ' num2str(I) ', Histology = ' fileh(I).name])
    imh{I} = [HISTO_INIT_P_DIR filesep name '.mha'];
    mask{I} = [HISTOLO_BF_P_DIR filesep 'mask-' name '.mha'];
%     mask{I} = [HISTO_INIT_P_DIR filesep 'mask-' name '.mha'];

end

% B-spline transform diffusion registration of the histology volume
clear optReg_bsp_hi optDiff_bsp_hi
optReg_bsp_hi.verbose = true;
optReg_bsp_hi.SpatialSamplesRatio = 0.1;
optReg_bsp_hi.MaxIter = 1;
optReg_bsp_hi.mask = mask;
optReg_bsp_hi.t0 = [];
optReg_bsp_hi.RegParam = [SRCDIR '/regParam-bspline-gray-intrahisto-hires.txt'];
optReg_bsp_hi.MaxVal = 0.0; % don't stop before reaching maximum number of iterations
optDiff_bsp_hi.Alpha = 0.45;
optDiff_bsp_hi.MaxIter = 1;
tic
[th2h_bsp_hi, infoReg_bsp_hi, infoDiff_bsp_hi, imout] = transfdiffreg(...
    'BSplineTransform', imh, optReg_bsp_hi, optDiff_bsp_hi);
toc

% save result
save(T_HISTO_BF, 'th2h_bsp_hi', 'infoReg_bsp_hi', 'infoDiff_bsp_hi', ...
    'optReg_bsp_hi', 'optDiff_bsp_hi')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% High res virtual slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load low-res alignment modified for high-res histology
load(T_HISTO_BF, 'ttot_lo2hi', 'infoReg_bsp_hi')
load(T_HISTOLO_BF, 'fileh')

% load high res B-spline diffusion transformations
load(ERRHI, 'maxIter_bsp_hi', 't_bsp_hi')

% number of images
N = length(fileh);

%% long axis and short axis

% select the B-spline refinement from all the refinements tried in
% study_stopping_criterion_bsp_refinement.m
for ITER = 1:length(maxIter_bsp_hi)
    
    disp(['Diff iter = ' num2str(maxIter_bsp_hi(ITER)) ...
        ', idx = ' num2str(ITER) ' / ' num2str(length(maxIter_bsp_hi))])
    
    for J = 1:N
        th2h_bsp_hi(J) = infoReg_bsp_hi.tp(1);
        th2h_bsp_hi(J).TransformParameters = t_bsp_hi{ITER}(J, :);
    end
    
    % total transform
    % ttot = elastix_cat(th2h_bsp_hi, ttot_lo2hi); % from scratch
    ttot = elastix_cat(th2h_bsp_hi); % only last B-spline hi-res refinement
    
    % apply high res refinement to high res histology
    optTransformix.AutoDefaultPixelValue = true;
    for J = 1:N
        
        tic
        disp(['J = ' num2str(J)])
        
        %% apply low res transformations to high res images
        
        % input file
        %     filein = [HISTO_DIR filesep fileh(J).name]; % from scratch
        filein = [HISTO_INIT_DIR filesep fileh(J).name]; % only last B-spline hi-res refinement
        
        % output filename
        optTransformix.outfile = [HISTO_BSP_DIR num2str(maxIter_bsp_hi(ITER)) filesep fileh(J).name];
        
        % apply low res reconstruction to high res histology
        transformix(ttot(J), filein, optTransformix);
        
        toc
        
    end
    
    clear file
    
    % create list of histology files with paths
    for I = 1:length(fileh)
        file{I} = [HISTO_BSP_DIR num2str(maxIter_bsp_hi(ITER)) filesep fileh(I).name];
    end
    
    tic
    aux = scimat_crop(file, [4289 nan nan; nan 6237 nan], ...
        [4289 nan nan; nan 6237 nan], 20e-6); % 14.8 min
    toc
    
    im = aux(1);
    im2 = aux(2);
    
    im.data = permute(im.data, [3 2 1 4 5]);
    im.axis = im.axis([3 2]);
    im.axis(1).spacing = im.axis(1).spacing * 1e3;
    im.axis(2).spacing = im.axis(2).spacing * 1e3;
    im.rotmat = eye(2);
    
    scimat_imagesc(im)
    set(gca, 'FontSize', 18)
    xlabel('x (mm)')
    ylabel('z (mm)')
    axis equal tight xy
    drawnow
    
    set(gcf, 'Color', 'white')
    
    export_fig(gca, [FIGDIR filesep 'virtual_slice_' MOUSE '_histo_hires_bsp_iter_' num2str(maxIter_bsp_hi(ITER)) '_LA_R4289.png'])
    
    im2.data = permute(im2.data, [3 1 2 4 5]);
    im2.axis = im2.axis([3 1]);
    im2.axis(1).spacing = im2.axis(1).spacing * 1e3;
    im2.axis(2).spacing = im2.axis(2).spacing * 1e3;
    im2.rotmat = eye(2);
    
    scimat_imagesc(im2)
    set(gca, 'FontSize', 18)
    xlabel('y (mm)')
    ylabel('z (mm)')
    axis equal tight xy
    drawnow
    
    set(gcf, 'Color', 'white')
    
    export_fig(gca, [FIGDIR filesep 'virtual_slice_' MOUSE '_histo_hires_bsp_iter_' num2str(maxIter_bsp_hi(ITER)) '_SAX_C6237.png'])

end

%% small details

ITER = 5; % 10 sweeps

tic
aux2 = scimat_crop(file, ...
    [4544 nan 33; 4577 nan 159; 7110 nan 165], ...
    [6136 nan 77; 7032 nan 221; 8669 nan 220], ...
    20e-6); % 14.8 min
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auxiliary code to fix resolution values in PNG header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for I = 1:N
% 
%     info = imfinfo([HISTDIR filesep fileh(I).name]);
%     infolo = imfinfo([HISTLODIR filesep fileh(I).name]);
%     
%     sWidth = info.Width / infolo.Width;
%     sHeight = info.Height / infolo.Height;
%     
%     str = [...
%         'mogrify -density ' ...
%         num2str(round(info.XResolution / sWidth / 100)) ...
%         'x' ...
%         num2str(round(info.YResolution / sHeight / 100)) ...
%         ' -units PixelsPerCentimeter ' ...
%         '''' HISTLODIR filesep fileh(I).name ''''];
%     
%     foo = system(str);
%     
% end
