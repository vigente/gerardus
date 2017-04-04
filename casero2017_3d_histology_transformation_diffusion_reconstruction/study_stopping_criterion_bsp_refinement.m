% study_stopping_criterion_bsp_refinement.m
%
% Script to check how landmark error changes with the number of stack
% sweeps for low-res and hi-res B-spline refinement.
%
% This script works with the output of histology_blockface_mouse.m, section
% "Intra-histology registration refinement: B-spline diffusion" (low-res)
% and section "Refine reconstruction in high resolution" (hi-res).

% version: 0.3.1
% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2016,2017 University of Oxford
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

%% variables

DEBUG = false;

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
HISTOLO_BF_P2_DIR = [HISTOLO_BF_DIR '_Bsp_Preprocessed'];

% transformations
T_HISTOLO_BF = [RESULT_DIR filesep 't_' HISTOLO '_to_Blockface.mat'];
T_HISTO_BF = [RESULT_DIR filesep 't_' HISTO '_to_Blockface.mat'];

% manual landmarks between blockface and histology
PTS_HIST2BF = [SRCDIR '/hand_tracing/pts_' HISTOLO '_to_Blockface.mat'];
PTS_HIST2HIST = [SRCDIR '/hand_tracing/pts_Intra_' HISTOLO '.mat'];

% landmark errors
ERR = [RESULT_DIR filesep 'err_' MOUSE '_' HISTOLO '.mat'];
ERRHI = [RESULT_DIR filesep 'err_' MOUSE '_' HISTO '.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% B-spline intra-histology refinement, 1 registration, several diffusions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we assume the registrations have already been run in script
% histology_blockface_mouse.m
load(T_HISTOLO_BF, 'th2bf', 'tCorrectHalvesMismatch', 'th2h_rig', ...
    'th2h_bsp', ...
    'infoReg_rig', 'infoDiff_rig', 'optReg_rig', 'optDiff_rig', ...
    'infoReg_bsp', 'infoDiff_bsp', 'optReg_bsp', 'optDiff_bsp', ...
    'filehBadCombinedIdx', 'filebf', 'fileh')

N = length(fileh);

%% compute B-spline refinement for different values of diffusion iterations
%% (the transformations computed here are used in the other sections of 
%% this script)

% concatenate B-spline parameters as plain matrix
tpMat = cat(1, infoReg_bsp.tp(:).TransformParameters);

% try different number of diffusion iterations
clear t_bsp infoTransfdiff_bsp
maxIter_bsp = [0:22 23:2:43 50:10:90 100:100:500];
for I = 2:length(maxIter_bsp)

    disp(['maxIter = ' num2str(maxIter_bsp(I))])
    
    % transform diffusion of B-spline parameters
    optDiff_bsp.Transform = 'BSplineTransform';
    optDiff_bsp.MaxIter = maxIter_bsp(I);
    optDiff_bsp.Epsilon = 0; % don't stop the diffusion because the update is too small
    optDiff_bsp.tbsp = infoReg_bsp.tp(1);
    [t_bsp{I}, infoTransfdiff_bsp{I}] = transfdiff(optDiff_bsp, tpMat);
    
end

save(ERR, 'maxIter_bsp', 't_bsp', 'infoTransfdiff_bsp', '-append')


%% blockface-histology landmark error

% load landmarks
load(PTS_HIST2BF, 'ptsbf', 'ptsh', 'filebfidx')

% load transformations
load(ERR, 'maxIter_bsp', 't_bsp')

% compute landmark error for each value of maximum number of diffusion
% iterations
errMed_bsp_bf2h = zeros(1, length(maxIter_bsp));
err95_bsp_bf2h = zeros(1, length(maxIter_bsp));
errMax_bsp_bf2h = zeros(1, length(maxIter_bsp));
for I = 1:length(maxIter_bsp)
    
    disp(['maxIter = ' num2str(maxIter_bsp(I))])

    if (I == 1)
        
        % compose transforms
        ttot = elastix_cat(th2h_rig, tCorrectHalvesMismatch, th2bf);
        
    else
    
        % transfer diffused parameters to elastix structs
        for J = 1:N
            th2h_bsp(J) = infoReg_bsp.tp(1);
            th2h_bsp(J).TransformParameters = t_bsp{I}(J, :);
        end
        
        % compose transforms
        ttot = elastix_cat(th2h_bsp, th2h_rig, tCorrectHalvesMismatch, th2bf);
        
    end
    
    % compute landmark errors in each slice
    err = [];
    for J = filebfidx % J in [1,239], as there are 239 slices in this stack
        
        % apply transform to landmarks
        ptsbf_h = transformix_pts(ttot(J), ptsbf{J});
        
        % compute landmark errors
        err = [err ; sqrt(sum((ptsbf_h - ptsh{J}).^2, 2))];
        
    end
    
    % median and maximum error values
    errMed_bsp_bf2h(I) = median(err);
    err95_bsp_bf2h(I) = prctile(err, 95);
    errMax_bsp_bf2h(I) = max(err);
    
end

save(ERR, 'errMed_bsp_bf2h', 'err95_bsp_bf2h', 'errMax_bsp_bf2h', '-append')

% plot errors
hold off
semilogx(maxIter_bsp, errMed_bsp_bf2h * 1e3, 'k', 'LineWidth', 2)
hold on
semilogx(maxIter_bsp, err95_bsp_bf2h * 1e3, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 18)
set(gca, 'XTick', [1 10 100 500]);
xlabel('Number of stack sweeps')
ylabel('Landmark distance error (mm)')
title('Blockface-histology')

saveas(gca, [FIGDIR filesep 'errtot_' MOUSE '_h2bf_bsp_' HISTO '.tif'])

%% intra-histology landmark error

% load landmarks
load(PTS_HIST2HIST, 'ptsh1', 'ptsh2', 'filehidx')

% load transformations
load(ERR, 'maxIter_bsp', 't_bsp')

% load transformed slice so that we have the reference grid of the
% transformed histology
imh = scimat_load([HISTOLO_BF_P2_DIR filesep fileh(1).name]);

% interpolation grid to invert the refinement
R = size(imh.data, 1);
C = size(imh.data, 2);
[gx_t, gy_t] = scimat_ndgrid(imh, 1:8:R, 1:8:C);

% compute landmark error for each value of maximum number of diffusion
% iterations
errMed_bsp_h2h = zeros(1, length(maxIter_bsp));
err95_bsp_h2h = zeros(1, length(maxIter_bsp));
errMax_bsp_h2h = zeros(1, length(maxIter_bsp));
for I = 1:length(maxIter_bsp)
    
    disp(['maxIter = ' num2str(maxIter_bsp(I))])

    %% compute total transformations
    
    if (maxIter_bsp(I) == 0)
        
        % compose transforms
        ttot = elastix_cat(th2h_rig, tCorrectHalvesMismatch, th2bf);
        
    else
    
        % transfer diffused parameters to elastix structs
        for J = 1:N
            th2h_bsp(J) = infoReg_bsp.tp(1);
            th2h_bsp(J).TransformParameters = t_bsp{I}(J, :);
        end
        
        % compose transforms
        ttot = elastix_cat(th2h_bsp, th2h_rig, tCorrectHalvesMismatch, th2bf);
        
    end
    
    % compute landmark errors between pairs of histology slices
    err = [];
    parfor K = 1:length(filehidx) 
        
        J = filehidx(K); % J in [1,239], as there are 239 slices in this stack
        
        fprintf('    J = %d... ', J)
        tic
        
        % check that we have landmarks
        if (isempty(ptsh1{J}))
            error(['No landmarks in ptsh1{' num2str(J) '}'])
        end
        % check that we have landmarks
        if (isempty(ptsh2{J}))
            error(['No landmarks in ptsh2{' num2str(J) '}'])
        end
        
        %% invert B-spline using interpolation
        
        % apply all steps of the reconstruction to the interpolation grid
        % Note: image transform original->reconstruction corresponds to
        % point transform reconstruction->original
        aux = transformix_pts(ttot(J), [gx_t(:) gy_t(:)]);
        gx1 = reshape(aux(:, 1), size(gx_t, 1), size(gx_t, 2));
        gy1 = reshape(aux(:, 2), size(gy_t, 1), size(gy_t, 2));
        
        aux = transformix_pts(ttot(J+1), [gx_t(:) gy_t(:)]);
        gx2 = reshape(aux(:, 1), size(gx_t, 1), size(gx_t, 2));
        gy2 = reshape(aux(:, 2), size(gy_t, 1), size(gy_t, 2));
        
        % map landmarks from the original slice to the transformed slice
        % using linear interpolation
        fx1 = scatteredInterpolant(gx1(:), gy1(:), gx_t(:), 'natural', 'linear');
        fy1 = scatteredInterpolant(gx1(:), gy1(:), gy_t(:), 'natural', 'linear');
        ptsh1_t = zeros(size(ptsh1{J}));
        ptsh1_t(:, 1) = fx1(ptsh1{J}(:, 1), ptsh1{J}(:, 2));
        ptsh1_t(:, 2) = fy1(ptsh1{J}(:, 1), ptsh1{J}(:, 2));

        fx2 = scatteredInterpolant(gx2(:), gy2(:), gx_t(:));
        fy2 = scatteredInterpolant(gx2(:), gy2(:), gy_t(:));
        ptsh2_t = zeros(size(ptsh2{J}));
        ptsh2_t(:, 1) = fx2(ptsh2{J}(:, 1), ptsh2{J}(:, 2));
        ptsh2_t(:, 2) = fy2(ptsh2{J}(:, 1), ptsh2{J}(:, 2));
       
        % compute landmark errors
        err = [err ; sqrt(sum((ptsh1_t - ptsh2_t).^2, 2))];

        fprintf('%.2f sec\n', toc)
        
    end
    
    % median and maximum error values
    errMed_bsp_h2h(I) = median(err);
    err95_bsp_h2h(I) = prctile(err, 95);
    errMax_bsp_h2h(I) = max(err);
    
end

save(ERR, 'errMed_bsp_h2h', 'err95_bsp_h2h', 'errMax_bsp_h2h', '-append')

% plot errors
hold off
semilogx(maxIter_bsp, errMed_bsp_h2h * 1e3, 'k', 'LineWidth', 2)
hold on
semilogx(maxIter_bsp, err95_bsp_h2h * 1e3, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 18)
set(gca, 'XTick', [1 10 100 500]);
xlabel('Number of stack sweeps')
ylabel('Landmark distance error (mm)')
title('Histology-histology')

% plot(maxIter_bsp, errMed_bsp_h2h * 1e3, 'r', 'LineWidth', 2)
% hold on
% plot(maxIter_bsp, err95_bsp_h2h * 1e3, 'r', 'LineWidth', 2)

saveas(gca, [FIGDIR filesep 'errtot_' MOUSE '_h2h_bsp_' HISTO '.tif'])

%% generate histology virtual slices

% load infoReg_bsp for transformix structs
load(T_HISTOLO_BF, 'infoReg_bsp')

% load transformations
load(ERR, 'maxIter_bsp', 't_bsp', 'optDiff_bsp')

% extrapolate background pixels when applying transformations to histology
% slices
opts.AutoDefaultPixelValue = true;

% try different number of diffusion iterations
% maxIter_bsp = [ 
% 0     1     2     3     4     5     6     7     8     9    10    11    12
% 13    14    15    16    17    18    19    20    21    22    23    25    27
% 29    31    33    35    37    39    41    43    50    60    70    80    90
% 100   200   300   400   500 ]
for I = [16 35 38 44] % [15 50 80 500]
    
    disp(['maxIter = ' num2str(maxIter_bsp(I))])
    optDiff_bsp.MaxIter = maxIter_bsp(I);

    % transfer diffused parameters to elastix structs
    for J = 1:N
        th2h_bsp(J) = infoReg_bsp.tp(1);
        th2h_bsp(J).TransformParameters = t_bsp{I}(J, :);
    end
    
    % compose transforms
    ttot = elastix_cat(th2h_bsp, th2h_rig, tCorrectHalvesMismatch, th2bf);
    
    % loop slices
    clear imh
    for J = 1:N
        
        % load histology
        [~, nameh] = fileparts(fileh(J).name);
        aux = scimat_load([HISTOLO_DIR filesep nameh '.tif']);
        
        % apply transforms to histology
        imh(J) = transformix(ttot(J), aux, opts);
        
    end
    
    % coordinates of image limits
    aux0 = scimat_index2world([1 1], imh(1));
    auxend = scimat_index2world([size(imh(1).data, 1) size(imh(1).data, 2)], imh(1));
    bbx = [aux0(1) auxend(1)]; % x-coordinates of the first and last corners of bounding box
    bby = [aux0(2) auxend(2)]; % y-coordinates of the first and last corners of bounding box
    bbz = [0 20e-6*(N-1)]; % z-coordinates of the first and last corners of bounding box
    
    % concatenate images slices into volume
    imhaux = cat(6, imh(:).data);
    
    % plot cross horizontal section
    hold off
    imagesc(bbx*1e3, bbz*1e3, permute(squeeze(imhaux(485, :, :, :, :, :)), [3 1 2]))
    if (I == 1)
        title('Histology pre-aligned to blockface')
    else
        title(['Histology after ' num2str(maxIter_bsp(I)) ' sweeps'])
    end
    set(gca, 'FontSize', 16)
    axis equal xy
    axis([2 12.5 0 20e-3*N])
    xlabel('x')
    ylabel('z')
    
    fig = gcf;
    set(fig, 'Color', 'white')
    
    export_fig(gca, [FIGDIR filesep 'stopping_criterion_' ...
        MOUSE '_histo_bsp_iter_' num2str(maxIter_bsp(I)) '_LA.png'])
    
    % plot cross vertical section
    hold off
    imagesc(bby*1e3, bbz*1e3, permute(squeeze(imhaux(:, 550, :, :, :, :)), [3 1 2]))
    set(gca, 'FontSize', 16)
    axis equal xy
    axis([2.5 10.75 0 20e-3*N])
    xlabel('y')
    ylabel('z')

    export_fig(gca, [FIGDIR filesep 'stopping_criterion_' ...
        MOUSE '_histo_bsp_iter_' num2str(maxIter_bsp(I)) '_SAX.png'])
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% B-spline intra-histology refinement, 27 registrations, 1 diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we assume that initial registrations have already been run in script
% histology_blockface_mouse.m
load(T_HISTOLO_BF, 'th2bf', 'tCorrectHalvesMismatch', 'th2h_rig', ...
    'th2h_bsp', ...
    'infoReg_rig', 'infoDiff_rig', 'optReg_rig', 'optDiff_rig', ...
    'infoReg_bsp_nodiff', 'infoDiff_bsp_nodiff', 'optReg_bsp_nodiff', 'optDiff_bsp_nodiff', ...
    'filehBadCombinedIdx', 'filebf', 'fileh')

N = length(fileh);

%% load histology after rigid refinement

% load histology slices and masks
clear imh mask
for I = 1:N
    
    disp(['I = ' num2str(I) ', Histology = ' fileh(I).name])
    imh(I) = scimat_load([HISTOLO_BF_P2_DIR filesep fileh(I).name]);
    mask(I) = scimat_load([HISTOLO_BF_P2_DIR filesep 'mask-' fileh(I).name]);
    
end

% B-spline transform diffusion registration of the histology volume
clear optReg_bsp optDiff_bsp
optReg_bsp_nodiff.verbose = true;
optReg_bsp_nodiff.SpatialSamplesRatio = 0.1;
optReg_bsp_nodiff.MaxIter = 27;
optReg_bsp_nodiff.mask = mask;
optReg_bsp_nodiff.t0 = [];
optReg_bsp_nodiff.RegParam = [SRCDIR '/regParam-bspline-gray-intrahisto.txt'];
optReg_bsp_nodiff.CacheFile = '/home/orie1416/Downloads/foo2.mat';
optReg_bsp_nodiff.MaxVal = 0.0; % don't stop before reaching maximum number of iterations
optDiff_bsp_nodiff.Alpha = 0.45;
optDiff_bsp_nodiff.MaxIter = 1;
tic
[th2h_bsp_nodiff, infoReg_bsp_nodiff, infoDiff_bsp_nodiff, imout] = transfdiffreg(...
    'BSplineTransform', imh, optReg_bsp_nodiff, optDiff_bsp_nodiff);
toc

% tic
% foo = optReg_bsp_nodiff;
% foo.mask = foo.mask(100:104);
% [th2h_bsp_nodiff, infoReg_bsp_nodiff, infoDiff_bsp_nodiff, imout] = transfdiffreg(...
%     'BSplineTransform', imh(100:104), foo, optDiff_bsp_nodiff);
% toc

% save result
save(T_HISTOLO_BF, 'th2h_bsp_nodiff', 'infoReg_bsp_nodiff', 'infoDiff_bsp_nodiff', ...
    'optReg_bsp_nodiff', 'optDiff_bsp_nodiff', '-append')


% load(T_HISTOLO_BF, 'th2h_bsp_nodiff', 'infoReg_bsp_nodiff', 'infoDiff_bsp_nodiff', ...
%     'optReg_bsp_nodiff', 'optDiff_bsp_nodiff')

foo = load('/home/orie1416/Downloads/foo2.mat')
th2h_bsp_nodiff = foo.ttot;
infoReg_bsp_nodiff = foo.info;
infoDiff_bsp_nodiff = foo.infoTransfdiff;

%% blockface-histology landmark error

% load landmarks
load(PTS_HIST2BF, 'ptsbf', 'ptsh', 'filebfidx')

% number of registration sweeps of the stack
NReg = elastix_length(th2h_bsp_nodiff);

% compute landmark error for each value of maximum number of diffusion
% iterations
errMed_bsp_nodiff_bf2h = zeros(1, NReg+1);
err95_bsp_nodiff_bf2h = zeros(1, NReg+1);
errMax_bsp_nodiff_bf2h = zeros(1, NReg+1);
for I = 0:NReg
    
    disp(['NReg = ' num2str(I)])

    if (I == 0) % alignment up to before B-spline refinement
        
        % compose transforms
        ttot = elastix_cat(th2h_rig, tCorrectHalvesMismatch, th2bf);
        
    else % add B-spline refinement
    
        % keep only the first I levels of B-splines
        th2h_bsp = elastix_colon(th2h_bsp_nodiff, (NReg-I+1):NReg);
        
        % compose transforms
        ttot = elastix_cat(th2h_bsp, th2h_rig, tCorrectHalvesMismatch, th2bf);
        
    end
    
    % compute landmark errors in each slice
    err = [];
    for J = filebfidx % J in [1,239], as there are 239 slices in this stack
        
        % apply transform to landmarks
        ptsbf_h = transformix_pts(ttot(J), ptsbf{J});
        
        % compute landmark errors
        err = [err ; sqrt(sum((ptsbf_h - ptsh{J}).^2, 2))];
        
    end
    
    % median and maximum error values
    errMed_bsp_nodiff_bf2h(I+1) = median(err);
    err95_bsp_nodiff_bf2h(I+1) = prctile(err, 95);
    errMax_bsp_nodiff_bf2h(I+1) = max(err);
    
end

save(ERR, 'errMed_bsp_nodiff_bf2h', 'err95_bsp_nodiff_bf2h', ...
    'errMax_bsp_nodiff_bf2h', '-append')

% plot errors
hold off
plot(0:NReg, errMed_bsp_nodiff_bf2h * 1e3, 'k', 'LineWidth', 2)
hold on
plot(0:NReg, err95_bsp_nodiff_bf2h * 1e3, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 18)
% set(gca, 'XTick', [1 10 100 500]);
xlabel('Number of stack sweeps')
ylabel('Landmark distance error (mm)')
title('Blockface-histology')

saveas(gca, [FIGDIR filesep 'errtot_' MOUSE '_h2bf_bsp_nodiff_' HISTO '.tif'])

%% intra-histology landmark error

% load landmarks
load(PTS_HIST2HIST, 'ptsh1', 'ptsh2', 'filehidx')

% load transformed slice so that we have the reference grid of the
% transformed histology
imh = scimat_load([HISTOLO_BF_P2_DIR filesep fileh(1).name]);

% interpolation grid to invert the refinement
R = size(imh.data, 1);
C = size(imh.data, 2);
[gx_t, gy_t] = scimat_ndgrid(imh, 1:8:R, 1:8:C);

% number of B-spline registration sweeps
NReg = elastix_length(th2h_bsp_nodiff);

% compute landmark error for each value of maximum number of diffusion
% iterations
errMed_bsp_nodiff_h2h = zeros(1, NReg+1);
err95_bsp_nodiff_h2h = zeros(1, NReg+1);
errMax_bsp_nodiff_h2h = zeros(1, NReg+1);
for I = 0:NReg
    
    disp(['NReg = ' num2str(I)])

    %% compute total transformations
    
    if (I == 0) % alignment up to before B-spline refinement
        
        % compose transforms
        ttot = elastix_cat(th2h_rig, tCorrectHalvesMismatch, th2bf);
        
    else % add B-spline refinement
    
        % keep only the first I levels of B-splines
        th2h_bsp = elastix_colon(th2h_bsp_nodiff, (NReg-I+1):NReg);
        
        % compose transforms
        ttot = elastix_cat(th2h_bsp, th2h_rig, tCorrectHalvesMismatch, th2bf);
        
    end
    
    % compute landmark errors between pairs of histology slices
    err = [];
    parfor K = 1:length(filehidx) 
        
        J = filehidx(K); % J in [1,239], as there are 239 slices in this stack
        
        fprintf('    J = %d... ', J)
        tic
        
        % check that we have landmarks
        if (isempty(ptsh1{J}))
            error(['No landmarks in ptsh1{' num2str(J) '}'])
        end
        % check that we have landmarks
        if (isempty(ptsh2{J}))
            error(['No landmarks in ptsh2{' num2str(J) '}'])
        end
        
        %% invert B-spline using interpolation
        
        % apply all steps of the reconstruction to the interpolation grid
        % Note: image transform original->reconstruction corresponds to
        % point transform reconstruction->original
        aux = transformix_pts(ttot(J), [gx_t(:) gy_t(:)]);
        gx1 = reshape(aux(:, 1), size(gx_t, 1), size(gx_t, 2));
        gy1 = reshape(aux(:, 2), size(gy_t, 1), size(gy_t, 2));
        
        aux = transformix_pts(ttot(J+1), [gx_t(:) gy_t(:)]);
        gx2 = reshape(aux(:, 1), size(gx_t, 1), size(gx_t, 2));
        gy2 = reshape(aux(:, 2), size(gy_t, 1), size(gy_t, 2));
        
        % map landmarks from the original slice to the transformed slice
        % using linear interpolation
        fx1 = scatteredInterpolant(gx1(:), gy1(:), gx_t(:), 'natural', 'linear');
        fy1 = scatteredInterpolant(gx1(:), gy1(:), gy_t(:), 'natural', 'linear');
        ptsh1_t = zeros(size(ptsh1{J}));
        ptsh1_t(:, 1) = fx1(ptsh1{J}(:, 1), ptsh1{J}(:, 2));
        ptsh1_t(:, 2) = fy1(ptsh1{J}(:, 1), ptsh1{J}(:, 2));

        fx2 = scatteredInterpolant(gx2(:), gy2(:), gx_t(:));
        fy2 = scatteredInterpolant(gx2(:), gy2(:), gy_t(:));
        ptsh2_t = zeros(size(ptsh2{J}));
        ptsh2_t(:, 1) = fx2(ptsh2{J}(:, 1), ptsh2{J}(:, 2));
        ptsh2_t(:, 2) = fy2(ptsh2{J}(:, 1), ptsh2{J}(:, 2));
       
        % compute landmark errors
        err = [err ; sqrt(sum((ptsh1_t - ptsh2_t).^2, 2))];

        fprintf('%.2f sec\n', toc)
        
    end
    
    % median and maximum error values
    errMed_bsp_nodiff_h2h(I+1) = median(err);
    err95_bsp_nodiff_h2h(I+1) = prctile(err, 95);
    errMax_bsp_nodiff_h2h(I+1) = max(err);
    
end

save(ERR, 'errMed_bsp_nodiff_h2h', 'err95_bsp_nodiff_h2h', ...
    'errMax_bsp_nodiff_h2h', '-append')

% plot errors
hold off
plot(0:NReg, errMed_bsp_nodiff_h2h * 1e3, 'k', 'LineWidth', 2)
hold on
plot(0:NReg, err95_bsp_nodiff_h2h * 1e3, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 18)
% set(gca, 'XTick', [1 10 100 500]);
xlabel('Number of stack sweeps')
ylabel('Landmark distance error (mm)')
title('Histology-histology')

saveas(gca, [FIGDIR filesep 'errtot_' MOUSE '_h2h_bsp_nodiff_' HISTO '.tif'])

%% generate histology virtual slices

load(T_HISTOLO_BF, 'th2h_bsp_nodiff', 'infoReg_bsp_nodiff', 'infoDiff_bsp_nodiff', ...
    'optReg_bsp_nodiff', 'optDiff_bsp_nodiff')

% number of registration sweeps of the stack
NReg = elastix_length(th2h_bsp_nodiff);

% extrapolate background pixels when applying transformations to histology
% slices
opts.AutoDefaultPixelValue = true;

% try different number of registration sweeps
for I = [5:2:15 20 27]
    
    disp(['NReg = ' num2str(I)])

    %% compute total transformations
    
    if (I == 0) % alignment up to before B-spline refinement
        
        % compose transforms
        ttot = elastix_cat(th2h_rig, tCorrectHalvesMismatch, th2bf);
        
    else % add B-spline refinement
    
        % keep only the first I levels of B-splines
        th2h_bsp = elastix_colon(th2h_bsp_nodiff, (NReg-I+1):NReg);
        
        % compose transforms
        ttot = elastix_cat(th2h_bsp, th2h_rig, tCorrectHalvesMismatch, th2bf);
        
    end

    % loop slices
    clear imh
    for J = 1:N
        
        % load histology
        [~, nameh] = fileparts(fileh(J).name);
        aux = scimat_load([HISTOLO_DIR filesep nameh '.tif']);
        
        % apply transforms to histology
        imh(J) = transformix(ttot(J), aux, opts);
        
    end
    
    % coordinates of image limits
    aux0 = scimat_index2world([1 1], imh(1));
    auxend = scimat_index2world([size(imh(1).data, 1) size(imh(1).data, 2)], imh(1));
    bbx = [aux0(1) auxend(1)]; % x-coordinates of the first and last corners of bounding box
    bby = [aux0(2) auxend(2)]; % y-coordinates of the first and last corners of bounding box
    bbz = [0 20e-6*(N-1)]; % z-coordinates of the first and last corners of bounding box
    
    % concatenate images slices into volume
    imhaux = cat(6, imh(:).data);
    
    % plot cross horizontal section
    figure
    hold off
    imagesc(bbx*1e3, bbz*1e3, permute(squeeze(imhaux(485, :, :, :, :, :)), [3 1 2]))
    if (I == 0)
        title('Histology after rigid refinement')
    else
        title(['Histology after ' num2str(I) ' sweeps'])
    end
    set(gca, 'FontSize', 16)
    axis equal xy
    axis([2 12.5 0 20e-3*N])
    xlabel('x')
    ylabel('z')
    
    fig = gcf;
    set(fig, 'Color', 'white')
    
    export_fig(gca, [FIGDIR filesep 'stopping_criterion_' ...
        MOUSE '_histo_bsp_nodiff_iter_' num2str(I+1) '_LA.png'])
    
    % plot cross vertical section
    figure
    hold off
    imagesc(bby*1e3, bbz*1e3, permute(squeeze(imhaux(:, 550, :, :, :, :)), [3 1 2]))
    set(gca, 'FontSize', 16)
    axis equal xy
    axis([2.5 10.75 0 20e-3*N])
    xlabel('y')
    ylabel('z')

    fig = gcf;
    set(fig, 'Color', 'white')
    
    export_fig(gca, [FIGDIR filesep 'stopping_criterion_' ...
        MOUSE '_histo_bsp_nodiff_iter_' num2str(I+1) '_SAX.png'])
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% combine diffusion and registration-only plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% plot errors
hold off
semilogx(maxIter_bsp, err95_bsp_bf2h * 1e3, 'k', 'LineWidth', 2)
hold on
semilogx(0:NReg, err95_bsp_nodiff_bf2h * 1e3, 'r', 'LineWidth', 2)
semilogx(maxIter_bsp, errMed_bsp_bf2h * 1e3, 'k', 'LineWidth', 2)
semilogx(0:NReg, errMed_bsp_nodiff_bf2h * 1e3, 'r', 'LineWidth', 2)
set(gca, 'FontSize', 18)
set(gca, 'XTick', [1 10 100 500]);
xlabel('Number of stack sweeps')
ylabel('Landmark distance error (mm)')
title('Blockface-histology')
legend('ATDR', 'Reg. only', 'Location', 'Best')

saveas(gca, [FIGDIR filesep 'fig8a-errtot_' MOUSE '_h2bf_bsp_combo_' HISTO '.tif'])

% cannot save as CSV for Elsevier's Interactive Plots different number of
% points for each plot. Instead, I've created a .csv file by hand using a
% spreadsheet

% plot errors
hold off
semilogx(maxIter_bsp, errMed_bsp_h2h * 1e3, 'k', 'LineWidth', 2)
hold on
semilogx(maxIter_bsp, err95_bsp_h2h * 1e3, 'k', 'LineWidth', 2)
semilogx(0:NReg, errMed_bsp_nodiff_h2h * 1e3, 'r', 'LineWidth', 2)
semilogx(0:NReg, err95_bsp_nodiff_h2h * 1e3, 'r', 'LineWidth', 2)
set(gca, 'FontSize', 18)
set(gca, 'XTick', [1 10 100 500]);
xlabel('Number of stack sweeps')
ylabel('Landmark distance error (mm)')
title('Histology-histology')

saveas(gca, [FIGDIR filesep 'errtot_' MOUSE '_h2h_bsp_combo_' HISTO '.tif'])

% cannot save as CSV for Elsevier's Interactive Plots different number of
% points for each plot. Instead, I've created a .csv file by hand using a
% spreadsheet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% High res B-spline intra-histology refinement, 1 registration, several 
%% diffusions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we assume the registrations have already been run in script
% histology_blockface_mouse.m
load(T_HISTOLO_BF, 'th2bf', 'tCorrectHalvesMismatch', 'th2h_rig', ...
    'th2h_bsp', ...
    'infoReg_rig', 'infoDiff_rig', 'optReg_rig', 'optDiff_rig', ...
    'infoReg_bsp', 'infoDiff_bsp', 'optReg_bsp', 'optDiff_bsp', ...
    'filehBadCombinedIdx', 'filebf', 'fileh')

load(ERR, 't_bsp', 'maxIter_bsp')

load(T_HISTO_BF, 'th2h_bsp_hi', 'infoReg_bsp_hi', 'infoDiff_bsp_hi', ...
    'optReg_bsp_hi', 'optDiff_bsp_hi')


N = length(fileh);

%% compute B-spline refinement for different values of diffusion iterations
%% (the transformations computed here are used in the other sections of 
%% this script)

% concatenate B-spline parameters as plain matrix
tpMat_hi = cat(1, infoReg_bsp_hi.tp(:).TransformParameters);

% try different number of diffusion iterations
clear t_bsp_hi infoTransfdiff_bsp_hi
maxIter_bsp_hi = [1 3 5 7 10 20 30 40 50 100 500 1000 10000]; % 100 iter -> 399 s, 500 iter -> 2000 s
for I = 1:length(maxIter_bsp_hi)

    tic
    disp(['maxIter = ' num2str(maxIter_bsp_hi(I))])
    
    % transform diffusion of B-spline parameters
    optDiff_bsp_hi.Transform = 'BSplineTransform';
    optDiff_bsp_hi.MaxIter = maxIter_bsp_hi(I);
    optDiff_bsp_hi.Epsilon = 0; % don't stop the diffusion because the update is too small
    optDiff_bsp_hi.tbsp = infoReg_bsp_hi.tp(1);
    [t_bsp_hi{I}, infoTransfdiff_bsp_hi{I}] = transfdiff(optDiff_bsp_hi, tpMat_hi);
    toc
    
end

save(ERRHI, 'maxIter_bsp_hi', 't_bsp_hi', 'infoTransfdiff_bsp_hi', '-v7.3')

