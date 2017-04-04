% histology_noblockface_mouse.m
%
% Script to reconstruct Mouse Q53 heart's low-res histology, without an
% external reference, to compare to "histology_blockface_mouse.m". This
% corresponds to experiment "2.2.1.	Low resolution reconstruction without
% external blockface reference" in the paper.

% version: 0.2.1
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

% histology directories for data
% HISTO = 'HistologyTrichrome';
HISTO = 'HistologySiriusRed';
HISTOLO = [HISTO 'LoRes'];
HISTO_DIR = [DATADIR filesep HISTO];
HISTOLO_DIR = [DATADIR filesep HISTOLO];

% pre-processing of the histology so that we can register the slices
% between them
HISTOLO_P_DIR = [IMPROC_DIR filesep HISTOLO '_Preprocessed'];

% transformations
T_HISTOLO_BF = [RESULT_DIR filesep 't_' HISTOLO '_to_Blockface.mat'];
T_HISTOLO_NOBF = [RESULT_DIR filesep 't_' HISTOLO '_without_Blockface.mat'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial rigid reconstruction. We start from the middle of the heart 
%% outwards
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% pre-processing of histology
%% It's faster to pre-process than transform 

% list of blockface and histology files
load(T_HISTOLO_BF, 'filebf', 'fileh')

N = length(fileh);

% use auto estimation of background colour in transformix
optTransformix.AutoDefaultPixelValue = true;

% compute a mean histogram from all the histology images. This will be the
% reference histogram for colour correction
colref = zeros(3, 256);
for I = 1:N

    disp(['I = ' num2str(I) ', Histology = ' fileh(I).name])
    
    % load original histology slice
    imh0 = scimat_load([HISTOLO_DIR filesep fileh(I).name]);
    
    % compute bin count histogram of image, separated by channel
    for CH = 1:3
        aux = imh0.data(:, :, :, :, CH);
        colref(CH, :) = colref(CH, :) + hist(aux(:), 0:255);
    end

end

% preprocess all images for registration
clear opts
opts.Grayscale = true;

% save to directory
for I = 1:N
    
    disp(['I = ' num2str(I) ', Histology = ' fileh(I).name])
    
    % load original histology slice
    imh0 = scimat_load([HISTOLO_DIR filesep fileh(I).name]);
    
    % preprocess the image
    [~, imh, mask] = histology_preprocessing(colref, imh0, opts);
    
    [~, filename] = fileparts(fileh(I).name);
    scimat_save([HISTOLO_P_DIR filesep filename '.mha'], imh);
    scimat_save([HISTOLO_P_DIR filesep 'mask-' filename '.mha'], mask);
    
end

%% sequential alignment of stack (here only compute registration of each 
%% slice to its one neighbour in the sequence)

% use auto estimation of background colour in transformix
optTransformix.AutoDefaultPixelValue = true;

% our center slices is going to be I = 121
ICENT = 121;

% align slice 1->2, 2->3, ..., 120->121
parfor I = 1:(ICENT-1)
    
    disp(['I = ' num2str(I)])
    
    % load two adjacent histology slices
    [~, filename0] = fileparts(fileh(I+1).name);
    imh0 = scimat_load([HISTOLO_P_DIR filesep filename0 '.mha']);
    [~, filename] = fileparts(fileh(I).name);
    imh = scimat_load([HISTOLO_P_DIR filesep filename '.mha']);
    
    subplot(2, 1, 1)
    hold off
    scimat_imagesc(imh0)
    subplot(2, 1, 2)
    hold off
    scimat_imagesc(imh)
    
    % rigid registration of current slice to next one
    th2h(I) = regmatchedfilt(imh0, imh, (-45:45)/180*pi);

%     % apply rigid transform to original histology
%     imh2 = transformix(th2h(I), imh, optTransformix);
%     
%     subplot(2, 1, 2)
%     imagesc(imfuse(imh0.data, imh2.data))

end

% load reference slice
imref = scimat_load([HISTOLO_DIR filesep fileh(ICENT).name]);

% no transformation for the central slice, and adjust the transformation to
% the pixel size and field of view of the reference slice
th2h(ICENT) = th2h(ICENT-1);
th2h(ICENT).TransformParameters = [0 0 0];
th2h(ICENT).Origin = scimat_index2world([1 1], imref);
th2h(ICENT).Size = [size(imref.data, 2) size(imref.data, 1)];
th2h(ICENT).Spacing([2 1]) = [imref.axis.spacing];

% align slice 121<-122, 122<-123, ..., 238<-239
parfor I = (ICENT+1):N
    
    disp(['I = ' num2str(I)])
    
    % load two adjacent histology slices
    [~, filename0] = fileparts(fileh(I-1).name);
    imh0 = scimat_load([HISTOLO_P_DIR filesep filename0 '.mha']);
    [~, filename] = fileparts(fileh(I).name);
    imh = scimat_load([HISTOLO_P_DIR filesep filename '.mha']);
    
    subplot(2, 1, 1)
    hold off
    scimat_imagesc(imh0)
    subplot(2, 1, 2)
    hold off
    scimat_imagesc(imh)
    
    % rigid registration of current slice to previous one
    th2h(I) = regmatchedfilt(imh0, imh, (-45:45)/180*pi);

%     % apply rigid transform to original histology
%     imh2 = transformix(th2h(I), imh, optTransformix);
%     
%     subplot(2, 1, 2)
%     imagesc(imfuse(imh0.data, imh2.data))

end

% save result
save(T_HISTOLO_NOBF, 'th2h')

%% accumulate registrations down the stack on both sides of reference slice

% init accumulate transformation on reference slices as no transformation
th2hAcc(ICENT) = th2h(ICENT);

% accumulate transforms 120->(121), 119->(120->121), ...
for I = (ICENT-1):-1:1
    
    disp(['I = ' num2str(I) ', Histology = ' fileh(I).name])
    
    % load original histology slice
    imh0 = scimat_load([HISTOLO_DIR filesep fileh(I).name]);

    % cumulative transform
    th2hAcc(I) = elastix_compose_afftransf(th2hAcc(I+1), th2h(I));

    % correct output field of view to the first frame
    th2hAcc(I).Origin = th2h(ICENT).Origin;
    th2hAcc(I).Spacing = th2h(ICENT).Spacing;
    th2hAcc(I).Size = th2h(ICENT).Size;
    
    imh(I) = transformix(th2hAcc(I), imh0, optTransformix);

end

% accumulate transforms (121)<-122, (121<-122)<-123,...
for I = (ICENT+1):N
    
    disp(['I = ' num2str(I) ', Histology = ' fileh(I).name])
    
    % load original histology slice
    imh0 = scimat_load([HISTOLO_DIR filesep fileh(I).name]);

    % cumulative transform
    th2hAcc(I) = elastix_compose_afftransf(th2hAcc(I-1), th2h(I));

    % correct output field of view to the first frame
    th2hAcc(I).Origin = th2h(ICENT).Origin;
    th2hAcc(I).Spacing = th2h(ICENT).Spacing;
    th2hAcc(I).Size = th2h(ICENT).Size;
    
    imh(I) = transformix(th2hAcc(I), imh0, optTransformix);

end

% save result
save(T_HISTOLO_NOBF, 'th2hAcc', '-append')


%% virtual slices

% coordinates of image limits
aux0 = scimat_index2world([1 1], imh(1));
auxend = scimat_index2world([size(imh(1).data, 1) size(imh(1).data, 2)], imh(1));
bbx = [aux0(1), auxend(1)]; % x-coordinates of the first and last corners of bounding box
bby = [aux0(2), auxend(2)]; % y-coordinates of the first and last corners of bounding box
bbz = 20e-6*[0 N-1]; % z-coordinates of each pixel

aux = cat(3, imh(:).data);

% plot cross horizontal section
% subplot(2, 1, 1)
hold off
imagesc(bbx * 1e3, bbz * 1e3, permute(aux(395, :, :, :, :, :), [3 2 5 1 4]))
set(gca, 'FontSize', 16)
title('Histology sequential pre-alignment')
xlabel('x (mm)')
ylabel('z (mm)')
axis equal tight xy

set(gcf, 'Color', 'white')
export_fig(gca, [FIGDIR filesep 'virtual_slice_' MOUSE '_histo_lores_nobf_LA_R395.png'])

% plot cross vertical section
% subplot(2, 1, 2)
hold off
imagesc(bby * 1e3, bbz * 1e3, permute(aux(:, 460, :, :, :, :), [3 1 5 2 4]))
set(gca, 'FontSize', 16)
xlabel('y (mm)')
ylabel('z (mm)')
axis equal tight xy

set(gcf, 'Color', 'white')
export_fig(gca, [FIGDIR filesep 'virtual_slice_' MOUSE '_histo_lores_nobf_SAX_C460.png'])

% for I = 140:20:780
%     imagesc(bby * 1e3, bbz * 1e3, permute(aux(:, I, :, :, :, :), [3 1 5 2 4]))
%     title(['I= ' num2str(I)])
%     axis xy equal
%     pause
% end
% 
% for I = 140:20:710
%     imagesc(bbx * 1e3, bbz * 1e3, permute(aux(I, :, :, :, :, :), [3 2 5 1 4]))
%     title(['I= ' num2str(I)])
%     axis xy equal
%     pause
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WE ARE HERE. Below this, code that needs rewriting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Intra-histology registration correction: rigid diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% use auto estimation of background colour in transformix
optTransformix.AutoDefaultPixelValue = true;

%% intra-histology correction

for I = 1:N
    
    % load all preprocessed histology and mask
    imh(I) = scimat_load([HISTLOPDIR filesep fileh(I).name]);
    
    % apply initial transformation
    imh(I) = transformix(th2h(I), imh(I), opts);
    
end

% rigid transform diffusion registration of the histology volume
clear optReg optDiff
optReg.Angle = (-10:.25:10) / 180 * pi;
optReg.verbose = true;
optDiff.MaxIter = 50;
tic
[taux, infoReg2, infoRigDiff2, imout] = transfdiffreg(...
    'regmatchedfilt', imh, optReg, optDiff);
toc

% concatenate transforms
th2h = elastix_cat(taux, th2h);

save([IMPROCDIR filesep T_NOBF_INTRAHISTLO], ...
    'th2h', 'infoReg2', 'infoRigDiff2', '-append')


%% visually check histology to histology alignment

clear imh
opts.AutoDefaultPixelValue = true;

% load transforms
load([IMPROCDIR filesep T_NOBF_INTRAHISTLO], 'th2h')
ttot = th2h;

% init memory to store transformed histology
imh = zeros([ttot(1).Size([2 1]) 1 1 3 N], 'uint8');

for I = 1:N
    
    % load histology
    [~, nameh] = fileparts(fileh(I).name);
    aux = scimat_load([HISTLODIR filesep nameh '.png']);
    
    % apply transforms to histology
    aux = transformix(ttot(I), aux, opts);
    imh(:, :, :, :, :, I) = aux.data;
    
end

% plot cross horizontal section
subplot(2, 1, 1)
hold off 
imagesc(permute(squeeze(imh(600, :, :, :, :, :)), [3 1 2]))
set(gca, 'FontSize', 16)

% plot cross vertical section
subplot(2, 1, 2)
hold off
imagesc(permute(squeeze(imh(:, 967, :, :, :, :)), [3 1 2]))
set(gca, 'FontSize', 16)

saveas(gca, [FIGDIR filesep INTRAHISTLORIGID '-cross-sections.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Intra-histology registration correction: B-spline diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear imh
opts.AutoDefaultPixelValue = true;

% compose transforms
load([IMPROCDIR filesep T_NOBF_INTRAHISTLO], 'th2h')
ttot = th2h;

% init memory to store transformed histology
[~, nameh] = fileparts(fileh(1).name);
imh = scimat_load([HISTLOPDIR filesep nameh '.tif']);
imh(1:N) = imh;

% apply previous transforms to pre-processed histology
for I = 1:N
    
    % load histology
    [~, nameh] = fileparts(fileh(I).name);
    imh(I) = scimat_load([HISTLOPDIR filesep nameh '.tif']);
    mask(I) = scimat_load([HISTLOPDIR filesep 'mask-' nameh '.tif']);
    
    % apply transforms to histology
    imh(I) = transformix(ttot(I), imh(I), opts);
    mask(I) = transformix(ttot(I), mask(I), opts);
    
end


%% intra-histology correction

% B-spline transform diffusion registration of the histology volume
clear optReg optDiff
optReg.verbose = true;
optReg.SpatialSamplesRatio = 0.1;
optReg.MaxIter = 4;
%optReg.MaxVal = 1;
optReg.mask = mask;
optReg.t0 = [];
optReg.RegParam = [SRCDIR '/regParam-bspline-gray-intrahisto.txt'];
optDiff.MaxIter = 5;
tic
[th2h_bsp, infoBspReg, infoBspDiff, imout] = transfdiffreg(...
    'BSplineTransform', imh, optReg, optDiff);
toc

% save result
save([IMPROCDIR filesep T_NOBF_INTRAHISTLO], ...
    'th2h_bsp', 'infoBspReg', 'infoBspDiff', '-append')

%% visually check histology to histology alignment

clear imh
opts.AutoDefaultPixelValue = true;

% compose transforms
load([IMPROCDIR filesep T_NOBF_INTRAHISTLO], 'th2h', 'th2h_bsp')
ttot = elastix_cat(th2h_bsp, th2h);

% init memory to store transformed histology
imh = zeros([ttot(1).Size([2 1]) 1 1 3 N], 'uint8');

for I = 1:N
    
    % load histology
    [~, nameh] = fileparts(fileh(I).name);
    aux = scimat_load([HISTLODIR filesep nameh '.png']);
    
    % apply transforms to histology
    aux = transformix(ttot(I), aux, opts);
    imh(:, :, :, :, :, I) = aux.data;
    
end

% plot cross horizontal section
subplot(2, 1, 1)
hold off 
imagesc(permute(squeeze(imh(600, :, :, :, :, :)), [3 1 2]))
set(gca, 'FontSize', 16)

% plot cross vertical section
subplot(2, 1, 2)
hold off
imagesc(permute(squeeze(imh(:, 967, :, :, :, :)), [3 1 2]))
set(gca, 'FontSize', 16)

saveas(gca, [FIGDIR filesep INTRAHISTLOBSP '-cross-sections.png'])
