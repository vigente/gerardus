% manual_landmarks.m
%
% Script to hand trace landmarks in pairs of
%   - blockface / histology slices
%   - histology / histology slices
% to validate the reconstruction algorithms.

% version: 0.3.3
% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2016-2017 University of Oxford
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

%% constants, paths, filenames...

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

% % binary masks for blockface stacks
% BF_MASKS = [IMPROC_DIR filesep 'Blockface_Masks.mat'];

% blockface external reference (jumps corrected, scratches corrected,
% perspective corrected, illumination corrected)
BFP_DIR = [IMPROC_DIR '/Blockface_Preprocessed'];

% histology directories for data
% HISTO = 'HistologyTrichrome';
HISTO = 'HistologySiriusRed';
HISTOLO = [HISTO 'LoRes'];
HISTO_DIR = [DATADIR filesep HISTO];
HISTOLO_DIR = [DATADIR filesep HISTOLO];

% files to save the landmarks
PTS_HIST2BF = [SRCDIR '/hand_tracing/pts_' HISTOLO '_to_Blockface.mat'];
PTS_HIST2HIST = [SRCDIR '/hand_tracing/' 'pts_Intra_' HISTOLO '.mat'];

% transformations
T_HISTOLO_BF = [RESULT_DIR filesep 't_' HISTOLO '_to_Blockface.mat'];

%% images where we are going to place landmarks on

load(T_HISTOLO_BF, 'fileAllIdx', 'filehBadCombinedIdx', 'filebf', 'fileh')

N = length(filebf);

%% generate histo-blockface landmarks by hand

% loop every 10-th slice
ptsh = cell(1, N);
ptsbf = cell(1, N);
filebf2hLandmarksIdx = nan(1, N);
filebfidx = [1:10:51 60:10:70 79:10:149 158:10:168 177 186 188 190 198 ...
    202 210 220 230 239];
for I = filebfidx
    
    % real slice index (it's useful to save them)
    filebf2hLandmarksIdx(I) = str2num(fileh(I).name(5:7));
    
    % load histology and blockface images. Histology has been downsized to
    % blockface size, so they have similar pixel size
    imh = scimat_load([HISTOLO_DIR filesep fileh(I).name]);
    imbf = scimat_im2scimat(...
        imread([BFP_DIR filesep filebf(I).name]), ...
        [imh.axis.spacing], ... % resolution
        [0 0] ... % spacing
        );

    % select matching landmarks
    [ptsh{I}, ptsbf{I}] = cpselect(squeeze(imh.data), imbf.data, 'Wait', true);
    
    % convert coordinates from row, colum to real world coordinates
    ptsh{I} = scimat_index2world(ptsh{end}(:, [2 1]), imh);
    ptsbf{I} = scimat_index2world(ptsbf{end}(:, [2 1]), imbf);
    
    % plot image with landmarks
    subplot(2, 1, 1)
    hold off
    scimat_imagesc(imbf);
    hold on
    plot(ptsbf{I}(:, 1), ptsbf{I}(:, 2), 'wo', 'LineWidth', 3)
    plot(ptsbf{I}(:, 1), ptsbf{I}(:, 2), 'kx', 'LineWidth', 2)
    
    subplot(2, 1, 2)
    hold off
    scimat_imagesc(imh);
    hold on
    plot(ptsh{I}(:, 1), ptsh{I}(:, 2), 'ko', 'LineWidth', 3)
    plot(ptsh{I}(:, 1), ptsh{I}(:, 2), 'wx', 'LineWidth', 2)
    
end

% visually check landmarks
for I = filebfidx
   
    % load histology and blockface images. Histology has been downsized to
    % blockface size, so they have similar pixel size
    imh = scimat_load([HISTOLO_DIR filesep fileh(I).name]);
    aux = scimat_load([BFP_DIR filesep filebf(I).name]);
    imbf = scimat_im2scimat(...
        aux.data, ...
        [imh.axis.spacing], ... % resolution
        [0 0] ... % spacing
        );

    % plot image with landmarks
    subplot(2, 1, 1)
    hold off
    scimat_imagesc(imbf);
    hold on
    plot(ptsbf{I}(:, 1), ptsbf{I}(:, 2), 'wo', 'LineWidth', 3)
    plot(ptsbf{I}(:, 1), ptsbf{I}(:, 2), 'kx', 'LineWidth', 2)
    title(['I = ' num2str(I)])
    
    subplot(2, 1, 2)
    hold off
    scimat_imagesc(imh);
    hold on
    plot(ptsh{I}(:, 1), ptsh{I}(:, 2), 'ko', 'LineWidth', 3)
    plot(ptsh{I}(:, 1), ptsh{I}(:, 2), 'wx', 'LineWidth', 2) 
    
    pause
    
end

% save landmarks
save(PTS_HIST2BF, 'ptsbf', 'ptsh', 'filebf2hLandmarksIdx', 'filebfidx')


%% generate histo-histo landmarks by hand

% loop every 5-th slice
ptsh1 = cell(1, N);
ptsh2 = cell(1, N);
fileh2hLandmarks1Idx = nan(1, N);
fileh2hLandmarks2Idx = nan(1, N);
filehidx = [1:5:56 60:5:70 74:5:149 158:10:168 177 185:10:N];
for I = filehidx
    
    % real slice index (it's useful to save them)
    fileh2hLandmarks1Idx(I) = str2num(fileh(I).name(5:7));
    fileh2hLandmarks2Idx(I) = str2num(fileh(I+1).name(5:7));

    % load two consecutive histology slices
    imh1 = scimat_load([HISTOLO_DIR filesep fileh(I).name]);
    imh2 = scimat_load([HISTOLO_DIR filesep fileh(I+1).name]);

    % select matching landmarks
    [ptsh1{I}, ptsh2{I}] = cpselect(squeeze(imh1.data), ...
        squeeze(imh2.data), 'Wait', true);
    
    % convert coordinates from row, colum to real world coordinates
    ptsh1{I} = scimat_index2world(ptsh1{end}(:, [2 1]), imh1);
    ptsh2{I} = scimat_index2world(ptsh2{end}(:, [2 1]), imh2);
    
    % plot image with landmarks
    subplot(2, 1, 1)
    hold off
    scimat_imagesc(imh1);
    hold on
    plot(ptsh1{I}(:, 1), ptsh1{I}(:, 2), 'wo', 'LineWidth', 3)
    plot(ptsh1{I}(:, 1), ptsh1{I}(:, 2), 'kx', 'LineWidth', 2)
    
    subplot(2, 1, 2)
    hold off
    scimat_imagesc(imh2);
    hold on
    plot(ptsh2{I}(:, 1), ptsh2{I}(:, 2), 'ko', 'LineWidth', 3)
    plot(ptsh2{I}(:, 1), ptsh2{I}(:, 2), 'wx', 'LineWidth', 2)
    
end

% visually check landmarks
for I = filehidx
   
    % load two consecutive histology slices
    imh1 = scimat_load([HISTOLO_DIR filesep fileh(I).name]);
    imh2 = scimat_load([HISTOLO_DIR filesep fileh(I+1).name]);
    
    % plot image with landmarks
    subplot(2, 1, 1)
    hold off
    scimat_imagesc(imh1);
    hold on
    plot(ptsh1{I}(:, 1), ptsh1{I}(:, 2), 'wo', 'LineWidth', 3)
    plot(ptsh1{I}(:, 1), ptsh1{I}(:, 2), 'kx', 'LineWidth', 2)
    title(['I = ' num2str(I)])
    
    subplot(2, 1, 2)
    hold off
    scimat_imagesc(imh2);
    hold on
    plot(ptsh2{I}(:, 1), ptsh2{I}(:, 2), 'ko', 'LineWidth', 3)
    plot(ptsh2{I}(:, 1), ptsh2{I}(:, 2), 'wx', 'LineWidth', 2) 
    
    pause
    
end

% save landmarks
save(PTS_HIST2HIST, 'ptsh1', 'ptsh2', 'fileh2hLandmarks1Idx', ...
    'fileh2hLandmarks2Idx', 'filehidx')
