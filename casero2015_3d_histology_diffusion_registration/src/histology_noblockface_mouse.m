% histology_noblockface_mouse.m
%
% Script to reconstruct one half of Mouse Q53 heart's low-res histology,
% without an external reference, to compare to
% "histology_blockface_mouse.m".
%
%   - Sequential pre-alignment of low-res histology.
%   - Rigid diffusion refinement of low-res intra-histology.
%   - B-spline diffusion refinement of low-res intra-histology 
%     (4 registration sweeps, 5 neighbor updates each).

% version: 0.1.3
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

% enable access to elastix shared object
libpath = getenv('LD_LIBRARY_PATH');
setenv('LD_LIBRARY_PATH', [ '/usr/local/bin/elastix/lib:' libpath]);

PROJDIR = '/home/orie1416/Software/private-gerardus-papers/casero2015_3d_histology_diffusion_registration';
DOCDIR = [PROJDIR '/doc'];
SRCDIR = [PROJDIR '/src'];
FIGDIR = [DOCDIR '/figures'];

MOUSE = 'Q53';

DATADIR = ['/data2/Mouse/' MOUSE];
IMPROCDIR = [DATADIR '/ImageProcessing'];
% HISTO = 'HistologyTrichrome';
HISTO = 'HistologySiriusRed';
HISTOLO = [HISTO 'LoRes']; % low-res histology
HISTLODIR = [DATADIR filesep HISTOLO];

T_NOBF_INTRAHISTLO = ['t_NoBlockfaceIntra_' HISTOLO '.mat'];

% pre-processed histology. Original histology slices, converted to
% grayscale, inverted, histograms matched...
HISTLOP = [HISTOLO 'Preprocessed'];
HISTLOPDIR = [IMPROCDIR filesep HISTLOP];

INTRAHISTLORIGID = ['NoBlockfaceIntra_' HISTOLO];
INTRAHISTLOBSP = ['NoBlockfaceIntra_' HISTOLO '_Bsp'];

% HISTREC = [HISTO 'Reconstructed'];
% HISTDIR = ['/data2/Mouse/' MOUSE filesep HISTO];
% HISTRECDIR = [IMPROCDIR filesep HISTREC];
% HISTLOREGDIR = [IMPROCDIR filesep HISTOLO '_to_BlockfaceCorrectedRegistration'];
% 
% HISTPDIR = [HISTRECDIR 'Preprocessed'];
% 
% HISTRECDIR_LVFREE = [HISTRECDIR '_LV_FreeWall'];
% HISTREFDIR_LVFREE = [HISTRECDIR 'Refined_LV_FreeWall'];
% 
% T_HISTLO2BF = ['t_' HISTOLO '_to_Blockface.mat'];
% T_TOTHIST_LVFREE = ['t_Total_' HISTO '_LV_FreeWall.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pre-processing of histology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% list of histology files
fileh = dir([HISTLODIR filesep MOUSE '*.png']);
N = length(fileh);

% convert to grayscale
opts.Grayscale = true;

% pre-process using I+1 as the reference to I
for I = N-1:-1:1
    
    % load two adjacent histology slices
    imh0 = scimat_load([HISTLODIR filesep fileh(I+1).name]);
    imh = scimat_load([HISTLODIR filesep fileh(I).name]);
    
    % preprocess histology to prepare it for registration
    [imh0, imh, mask] = histology_preprocessing(imh0, imh, opts);

    % save preprocessed images
    scimat_save([HISTLOPDIR filesep fileh(I).name], imh);
    scimat_save([HISTLOPDIR filesep 'mask-' fileh(I).name], mask);
    
end

% create a pre-processing for the last slice, but using N-1 as the
% reference
% load two adjacent histology slices
imh0 = scimat_load([HISTLODIR filesep fileh(end-1).name]);
imh = scimat_load([HISTLODIR filesep fileh(end).name]);

% preprocess histology to prepare it for registration
[imh0, imh, mask] = histology_preprocessing(imh0, imh, opts);

% save preprocessed images
scimat_save([HISTLOPDIR filesep fileh(end).name], imh);
scimat_save([HISTLOPDIR filesep 'mask-' fileh(end).name], mask);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial rigid reconstruction. We start from the last histology slice, N.
%% Then register N-1 to N. Then N-2 to N-1, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% sequential rigid registration starting from the largest slice

% use auto estimation of background colour in transformix
optTransformix.AutoDefaultPixelValue = true;

imh = scimat_load([HISTLOPDIR filesep fileh(1).name]);
imh(1:N) = imh;

for I = N-1:-1:1
    
    % load two adjacent histology slices
    imh0 = scimat_load([HISTLOPDIR filesep fileh(I+1).name]);
    imh = scimat_load([HISTLOPDIR filesep fileh(I).name]);
    
    subplot(2, 1, 1)
    hold off
    scimat_imagesc(imh0)
    subplot(2, 1, 2)
    hold off
    scimat_imagesc(imh)
    
    % apply transform from previous iteration to the reference slice
    if (I < N-1)
        
        imh0 = transformix(th2h(I+1), imh0, optTransformix);
        
    end
    
    % rigid registration of current slice to next one
    th2h(I) = regmatchedfilt(imh0, imh, (-45:45)/180*pi);

    % apply rigid transform to original histology
    imh2 = transformix(th2h(I), imh, optTransformix);
    
    subplot(2, 1, 2)
    imagesc(imfuse(imh0.data, imh2.data))

end

% no transformation for the biggest slice
th2h(N) = th2h(N-1);
th2h(N).TransformParameters = [0 0 0];

% save result
save([IMPROCDIR filesep T_NOBF_INTRAHISTLO], 'th2h')

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
