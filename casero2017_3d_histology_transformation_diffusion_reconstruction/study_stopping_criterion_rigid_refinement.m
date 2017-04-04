% study_stopping_criterion_rigid_refinement.m
%
% Script to check how landmark error changes with the number of stack
% sweeps for low-res rigid refinement.
%
% This script works with the output of histology_blockface_mouse.m, section
% "Intra-histology refinement: rigid diffusion".

% version: 0.1.7
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

%% variables

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

% transformations
T_HISTOLO_BF = [RESULT_DIR filesep 't_' HISTOLO '_to_Blockface.mat'];

% manual landmarks between blockface and histology
PTS_HIST2BF = [SRCDIR '/hand_tracing/pts_' HISTOLO '_to_Blockface.mat'];
PTS_HIST2HIST = [SRCDIR '/hand_tracing/pts_Intra_' HISTOLO '.mat'];

% landmark errors
ERR = [RESULT_DIR filesep 'err_' MOUSE '_' HISTOLO '.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% intra-histology refinement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we assume the regstrations have already been run in script
% histology_blockface_mouse.m
load(T_HISTOLO_BF, 'th2bf', 'tCorrectHalvesMismatch', 'th2h_rig', ...
    'infoReg_rig', 'infoDiff_rig', 'optReg_rig', 'optDiff_rig', ...
    'filehBadCombinedIdx', 'filebf', 'fileh')

N = length(fileh);

opts.AutoDefaultPixelValue = true;

%% compute affine refinement for different values of diffusion iterations
%% (the transformations computed here are used in the other sections of 
%% this script)

% convert elastix parameter vectors to affine matrices
tpMat = elastix_affine_struct2matrix(infoReg_rig.tp);

% contant parameters for the diffusion
optDiff_rig.Transform = 'EulerTransform';
optDiff_rig.Epsilon = 0.0; % don't stop before reaching maximum number of iterations
optDiff_rig.Alpha = 0.45;

% try different number of diffusion iterations
maxIter = [0 1:20 30:10:90 100:50:600 700:100:1600 2000:1000:5000];
for I = 2:length(maxIter)

    disp(['maxIter = ' num2str(maxIter(I))])
    
    % transform diffusion of affine parameters
    optDiff_rig.MaxIter = maxIter(I);
    [t_rig{I}, infoTransfdiff_rig{I}] = transfdiff(optDiff_rig, tpMat);
    
end

save(ERR, 'maxIter', 't_rig', 'infoTransfdiff_rig', 'optDiff_rig')


%% blockface-histology landmark error

% load landmarks
load(PTS_HIST2BF, 'ptsbf', 'ptsh', 'filebfidx')

% load transformations
load(ERR, 'maxIter', 't_rig')

% compute landmark error for each value of maximum number of diffusion
% iterations
errMed_rig_bf2h = zeros(1, length(maxIter));
err95_rig_bf2h = zeros(1, length(maxIter));
errMax_rig_bf2h = zeros(1, length(maxIter));
for I = 2:length(maxIter)
    
    disp(['maxIter = ' num2str(maxIter(I))])

    % transfer diffused parameters to elastix structs
    th2h_rig = elastix_affine_matrix2struct(t_rig{I}, infoReg_rig.tp(1));
    
    % compose transforms
    ttot = elastix_cat(th2h_rig, tCorrectHalvesMismatch, th2bf);
    
    % compute landmark errors in each slice
    err = [];
    for J = filebfidx % J in [1,239], as there are 239 slices in this stack
        
        % apply transform to landmarks
        ptsbf_h = transformix_pts(ttot(J), ptsbf{J});
        
        % compute landmark errors
        err = [err ; sqrt(sum((ptsbf_h - ptsh{J}).^2, 2))];
        
    end
    
    % median and maximum error values
    errMed_rig_bf2h(I) = median(err);
    err95_rig_bf2h(I) = prctile(err, 95);
    errMax_rig_bf2h(I) = max(err);
    
end

% reference error before alignment
t_rig{1} = repmat(eye(3, 3), 1, 1, N);
err = [];
for J = filebfidx
    err = [err; sqrt(sum((ptsbf{J} - ptsh{J}).^2, 2))];
end
errMed_rig_bf2h(1) = median(err);
err95_rig_bf2h(1) = prctile(err, 95);
errMax_rig_bf2h(1) = max(err);

save(ERR, 'errMed_rig_bf2h', 'err95_rig_bf2h', 'errMax_rig_bf2h', '-append')

% plot errors
hold off
semilogx(maxIter, errMed_rig_bf2h * 1e3, 'k', 'LineWidth', 2)
hold on
semilogx(maxIter, err95_rig_bf2h * 1e3, 'k', 'LineWidth', 2)
% semilogx(maxIter, errMax_rig_bf2h, 'k')
set(gca, 'FontSize', 18)
set(gca, 'XTick', [1 10 100 1000 5000]);
xlabel('Number of stack sweeps')
ylabel('Landmark distance error (mm)')
title('Blockface-histology')

saveas(gca, [FIGDIR filesep 'fig7a-errtot_' MOUSE '_h2bf_rig_' HISTO '.tif'])

% save as CSV for Elsevier's Interactive Plots
% too many rows for Elsevier's visualizer
csvfile = [FIGDIR filesep 'fig7a-errtot_' MOUSE '_h2bf_rig_' HISTO '.csv'];
fid = fopen(csvfile, 'w');
fprintf(fid, '%s,%s,%s\n', ...
    'stack sweeps', 'Median landmark distance error (mm)', '95% landmark distance error (mm)');
fclose(fid);   
dlmwrite(csvfile, ...
    [maxIter', errMed_rig_bf2h' * 1e3, err95_rig_bf2h' * 1e3], '-append', 'delimiter', ',','precision', '%.3E');
system(['sed -i ''s/E+/E/g'' ' csvfile]);


%% intra-histology landmark error

% load landmarks
load(PTS_HIST2HIST, 'ptsh1', 'ptsh2', 'filehidx')

% load transformations
load(ERR, 'maxIter', 't_rig')

% compute landmark error for each value of maximum number of diffusion
% iterations
errMed_rig_h2h = zeros(1, length(maxIter));
err95_rig_h2h = zeros(1, length(maxIter));
errMax_rig_h2h = zeros(1, length(maxIter));
for I = 2:length(maxIter)
    
    disp(['maxIter = ' num2str(maxIter(I))])

    % transfer diffused parameters to elastix structs
    th2h_rig = elastix_affine_matrix2struct(t_rig{I}, infoReg_rig.tp(1));
    
    % compose transforms
    ttot = elastix_cat(th2h_rig, tCorrectHalvesMismatch, th2bf);
    
    % compute landmark errors between pairs of histology slices
    err = [];
    for J = filehidx % J in [1,239], as there are 239 slices in this stack
        
        % 1st slice of the pair
        
        % merge the three transformations into single transformation
        ttotComb = elastix_compose_afftransf(...
            ttot(J).InitialTransformParametersFileName.InitialTransformParametersFileName, ...
            ttot(J).InitialTransformParametersFileName);
        ttot(J).InitialTransformParametersFileName = 'NoInitialTransform';
        ttotComb = elastix_compose_afftransf(ttotComb, ttot(J));
        
        % compute inverse transformation, because this one maps
        % histology->blockface, i.e. blockface points->histology points,
        % but we have histology points
        A = elastix_affine_struct2matrix(ttotComb);
        ttotCombInv1 = elastix_affine_matrix2struct(inv(A), ttotComb);
        
        % 2nd slice of the pair
        
        % merge the three transformations into single transformation
        ttotComb = elastix_compose_afftransf(...
            ttot(J+1).InitialTransformParametersFileName.InitialTransformParametersFileName, ...
            ttot(J+1).InitialTransformParametersFileName);
        ttot(J+1).InitialTransformParametersFileName = 'NoInitialTransform';
        ttotComb = elastix_compose_afftransf(ttotComb, ttot(J+1));
        
        % compute inverse transformation, because this one maps
        % histology->blockface, i.e. blockface points->histology points,
        % but we have histology points
        A = elastix_affine_struct2matrix(ttotComb);
        ttotCombInv2 = elastix_affine_matrix2struct(inv(A), ttotComb);
        
        % check that we have landmarks
        if (isempty(ptsh1{J}))
            error(['No landmarks in ptsh1{' num2str(J) '}'])
        end
        % check that we have landmarks
        if (isempty(ptsh2{J}))
            error(['No landmarks in ptsh2{' num2str(J) '}'])
        end
        
        % apply inverse transform to histology landmarks
        ptsh1_mapped = transformix_pts(ttotCombInv1, ptsh1{J});
        ptsh2_mapped = transformix_pts(ttotCombInv2, ptsh2{J});
        
        % compute landmark errors
        err = [err ; sqrt(sum((ptsh1_mapped - ptsh2_mapped).^2, 2))];
        
    end
    
    % median and maximum error values
    errMed_rig_h2h(I) = median(err);
    err95_rig_h2h(I) = prctile(err, 95);
    errMax_rig_h2h(I) = max(err);
    
end

% reference error before alignment
err = [];
for J = filehidx
    err = [err; sqrt(sum((ptsh1{J} - ptsh2{J}).^2, 2))];
end
errMed_rig_h2h(1) = median(err);
err95_rig_h2h(1) = prctile(err, 95);
errMax_rig_h2h(1) = max(err);

save(ERR, 'errMed_rig_h2h', 'err95_rig_h2h', 'errMax_rig_h2h', '-append')

% plot errors
hold off
semilogx(maxIter, errMed_rig_h2h * 1e3, 'k', 'LineWidth', 2)
hold on
semilogx(maxIter, err95_rig_h2h * 1e3, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 18)
set(gca, 'XTick', [1 10 100 1000 5000]);
xlabel('Number of stack sweeps')
ylabel('Landmark distance error (mm)')
title('Histology-histology')

saveas(gca, [FIGDIR filesep 'fig7a-errtot_' MOUSE '_h2h_rig_' HISTO '.tif'])

% save as CSV for Elsevier's Interactive Plots
% too many rows for Elsevier's visualizer
csvfile = [FIGDIR filesep 'fig7a-errtot_' MOUSE '_h2h_rig_' HISTO '.csv'];
fid = fopen(csvfile, 'w');
fprintf(fid, '%s,%s,%s\n', ...
    'stack sweeps', 'Median landmark distance error (mm)', '95% landmark distance error (mm)');
fclose(fid);   
dlmwrite(csvfile, ...
    [maxIter', errMed_rig_h2h' * 1e3, err95_rig_h2h' * 1e3], '-append', 'delimiter', ',','precision', '%.3E');
system(['sed -i ''s/E+/E/g'' ' csvfile]);

%% generate histology virtual slices

% load infoReg_rig for transformix structs
load(T_HISTOLO_BF, 'infoReg_rig')

% load transformations
load(ERR, 'maxIter', 't_rig', 'optDiff_rig')

% try different number of diffusion iterations
% maxIter = [ 0     1     2     3     4     5     6     7     8     9    10
% 11    12    13    14    15    16    17    18    19    20    30    40
% 50    60    70    80    90   100   150   200   250   300   350   400
% 450   500   550   600   700   800   900  1000  1100  1200  1300  1400 
% 1500  1600  2000  3000  4000  5000 ]
for I = [1 4 30 53] % [0, 3, 150, 5000]
    
    disp(['maxIter = ' num2str(maxIter(I))])
    optDiff_rig.MaxIter = maxIter(I);

    % transfer diffused parameters to elastix structs
    th2h_rig = elastix_affine_matrix2struct(t_rig{I}, infoReg_rig.tp(1));
    
    % compose transforms
    ttot = elastix_cat(th2h_rig, tCorrectHalvesMismatch, th2bf);
    
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
        title(['Histology after ' num2str(maxIter(I)) ' sweeps'])
    end
    set(gca, 'FontSize', 16)
    axis equal xy
    axis([2 12.5 0 20e-3*N])
    xlabel('x')
    ylabel('z')
    
    fig = gcf;
    set(fig, 'Color', 'white')
    
    export_fig(gca, [FIGDIR filesep 'stopping_criterion_' ...
        MOUSE '_histo_iter_' num2str(maxIter(I)) '_LA.png'])
    
    % plot cross vertical section
    hold off
    imagesc(bby*1e3, bbz*1e3, permute(squeeze(imhaux(:, 550, :, :, :, :)), [3 1 2]))
    set(gca, 'FontSize', 16)
    axis equal xy
    axis([2.5 10.75 0 20e-3*N])
    xlabel('y')
    ylabel('z')

    export_fig(gca, [FIGDIR filesep 'stopping_criterion_' ...
        MOUSE '_histo_iter_' num2str(maxIter(I)) '_SAX.png'])
    
end

%% generate blockface virtual slices

N = length(filebf);

% loop slices
imbf = zeros(971, 1099, N, 'uint8');
for J = 1:N
    
    % load histology
    aux = scimat_load([BFP_DIR filesep filebf(J).name]);
    
    imbf(:, :, J) = aux.data;
    
end

% Note: we are going to use the axes from 

% plot cross horizontal section
hold off
imagesc(bbx*1e3, bbz*1e3, permute(squeeze(imbf(485, :, :)), [2 1]))
title('Blockface')
set(gca, 'FontSize', 16)
axis equal tight xy
axis([2 12.5 0 20e-3*N])
xlabel('x (mm)')
ylabel('z (mm)')

export_fig(gca, [FIGDIR filesep 'virtual_slice_' ...
    MOUSE '_bf_LA_R485.png'])

% plot cross vertical section
hold off
imagesc(bby*1e3, bbz*1e3, permute(squeeze(imbf(:, 550, :)), [2 1]))
set(gca, 'FontSize', 16)
axis equal tight xy
axis([2.5 10.75 0 20e-3*N])
xlabel('y (mm)')
ylabel('z (mm)')

export_fig(gca, [FIGDIR filesep 'virtual_slice_' ...
    MOUSE '_bf_SAX_C550.png'])

