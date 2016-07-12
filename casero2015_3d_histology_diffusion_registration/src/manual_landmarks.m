% manual_landmarks.m
%
% Script to hand trace landmarks in pairs of
%   - blockface / histology slices
%   - histology / histology slices
% to validate the reconstruction algorithms.

% version: 0.2.2
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

%% constants, paths, filenames...

DOCDIR = '/home/orie1416/Software/private-gerardus-papers/casero2015_3d_histology_diffusion_registration/doc';
SRCDIR = [DOCDIR '/../src'];
FIGDIR = [DOCDIR '/figures'];

MOUSE = 'Q53';

HISTO = 'HistologySiriusRed';
HISTOLO = [HISTO 'LoRes'];
HISTDIR = ['/data2/Mouse/' MOUSE filesep HISTO];
HISTLODIR = ['/data2/Mouse/' MOUSE filesep HISTOLO];
IMPROCDIR = ['/data2/Mouse/' MOUSE '/ImageProcessing'];
BFCDATADIR = [IMPROCDIR '/BlockfaceCorrected']; % alignment only
BFPDATADIR = [IMPROCDIR '/BlockfacePreprocessed']; % cropped, intensity equalisation, inversion

PTS_HIST2BF = ['pts_' HISTOLO '_to_Blockface.mat'];
PTS_HIST2HIST = ['pts_Intra_' HISTOLO '.mat'];

%% images where we are going to place landmarks on

% list of blockface and histology files
[filebf, fileh] = match_file_lists(...
    dir([BFCDATADIR filesep '*.png']), ...
    dir([HISTLODIR filesep '*.tif']), ...
    9:11, ...
    5:7 ...
    );
N = length(filebf);

%% generate histo-blockface landmarks by hand

% crop factors for blockface images
crop_row = 284:1254;
crop_col = 766:1864;

% size scaling = histology / blockface
% this is estimated by measuring the same antomical distance in a blockface
% image and the corresponding histology image
sc = 866.5806 / 594.4855;

% loop every 5-th slice
for I = 1:5:N
    
    % load histology and blockface images. The blockface images are going
    % to be scaled to histology's dimensions, so we are going to assume
    % that they have the same pixel size (even if they don't)
    imh = scimat_load([HISTLODIR filesep fileh(I).name]);
    imbf = scimat_im2scimat(...
        imread([BFCDATADIR filesep filebf(I).name]), ...
        [imh.axis.spacing], ... % resolution
        [0 0] ... % spacing
        );

    % crop blockface image
    imbf.data = imbf.data(crop_row, crop_col);

    % upscale blockface so that it's in the same scale as histology
    imbf.data = imresize(imbf.data, sc);
    
    % correct blockface image size in the metadata
    imbf.axis(1).size = size(imbf.data, 1);
    imbf.axis(2).size = size(imbf.data, 2);
    
    % select matching landmarks
    [ptsh{I}, ptsbf{I}] = cpselect(squeeze(imh.data), imbf.data, 'Wait', true);
    
    % convert coordinates from row, colum to real world coordinates
    ptsh{I} = scimat_index2world(ptsh{I}(:, [2 1]), imh);
    ptsbf{I} = scimat_index2world(ptsbf{I}(:, [2 1]), imbf);
    
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

% save landmarks
save([IMPROCDIR filesep PTS_HIST2BF], 'ptsbf', 'ptsh')

%% generate histo-histo landmarks by hand

clear ptsh ptsbf

% loop every 5-th slice
for I = 1:5:N-1
    
    % load histology and blockface images. The blockface images are going
    % to be scaled to histology's dimensions, so we are going to assume
    % that they have the same pixel size (even if they don't)
    imh1 = scimat_load([HISTLODIR filesep fileh(I).name]);
    imh2 = scimat_load([HISTLODIR filesep fileh(I+1).name]);

    % select matching landmarks
    [ptsh1{I}, ptsh2{I}] = cpselect(squeeze(imh1.data), ...
        squeeze(imh2.data), 'Wait', true);
    
    % convert coordinates from row, colum to real world coordinates
    ptsh1{I} = scimat_index2world(ptsh1{I}(:, [2 1]), imh1);
    ptsh2{I} = scimat_index2world(ptsh2{I}(:, [2 1]), imh2);
    
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

% save landmarks
save([IMPROCDIR filesep PTS_HIST2HIST], 'ptsh1', 'ptsh2')
