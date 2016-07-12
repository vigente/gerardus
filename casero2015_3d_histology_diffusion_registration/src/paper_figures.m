% paper_figures.m
%
% Script to generate the plots/figures used in the paper.
%
% Version: 0.2.18
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

% laptop at home
BASEDIR = '/media/rcasero/mvldata/data';
PROJDIR = '/home/rcasero/Software/casero2015_3d_histology_diffusion_registration';

% workstation at the office
BASEDIR = '/data2';
PROJDIR = '/home/orie1416/Software/casero2015_3d_histology_diffusion_registration';

% documentation directories
DOCDIR = [PROJDIR filesep 'doc'];
FIGDIR = [DOCDIR filesep 'figures'];
SRCDIR = [PROJDIR filesep 'src'];

MOUSE = 'Q53';

IMPROCDIR = [BASEDIR '/Mouse/' MOUSE '/ImageProcessing'];
BFDATADIR = [BASEDIR '/Mouse/' MOUSE '/Blockface'];
BFSDATADIR = [IMPROCDIR '/BlockfaceStabilised'];
BFCDATADIR = [IMPROCDIR '/BlockfaceCorrected'];
BFPDATADIR = [IMPROCDIR '/BlockfacePreprocessed'];
% HISTO = 'HistologyTrichrome';
HISTO = 'HistologySiriusRed';
HISTOLO = [HISTO 'LoRes'];
HISTLODIR = [BASEDIR '/Mouse/' MOUSE filesep HISTOLO];
HISTDIR = [BASEDIR '/Mouse/' MOUSE filesep HISTO];
HISTLOREGDIR = [BASEDIR '/Mouse/' MOUSE filesep HISTOLO '_to_BlockfaceCorrectedRegistration'];
%HISTLO2BFREGSIMILARITYDIR = [IMPROCDIR filesep HISTOLO 'GrayToBlockfaceRegistrationSimilarity'];
HISTLO2BFREGRIGID = [HISTOLO '_to_BlockfaceRegistrationRigid'];
HISTLO2BFREGRIGIDDIR = [IMPROCDIR filesep HISTLO2BFREGRIGID];
INTRAHISTLORIGID = ['Intra_' HISTOLO];
INTRAHISTLOBSP = ['Intra_' HISTOLO '_Bsp'];
INTRANOBFHISTLORIGID = ['NoBlockfaceIntra_' HISTOLO];

HISTLOPDIR = [HISTLO2BFREGRIGIDDIR 'Preprocessed'];

HISTREC = [HISTO 'Reconstructed'];
HISTREC_LVFREE = [HISTREC '_LV_FreeWall'];
HISTREF_LVFREE = [HISTREC 'Refined_LV_FreeWall'];
HISTRECDIR = [IMPROCDIR filesep HISTREC];
HISTRECDIR_LVFREE = [HISTRECDIR '_LV_FreeWall'];
HISTREFDIR_LVFREE = [HISTRECDIR 'Refined_LV_FreeWall'];

% initial alignment of histology to blockface, common to all subsequent
% refinements
T_HISTLO2BF = ['t_' HISTOLO '_to_Blockface.mat'];
T_HISTLO2BF_FILE = [IMPROCDIR filesep T_HISTLO2BF];
T_INTRAHISTLO = ['t_Intra_' HISTOLO '.mat'];
T_NOBF_INTRAHISTLO = ['t_NoBlockfaceIntra_' HISTOLO '.mat'];
T_TOTHIST_LVFREE = ['t_Total_' HISTO '_LV_FreeWall.mat'];

% manually traced landmarks
PTS_HIST2BF = ['pts_' HISTOLO '_to_Blockface.mat'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% repeat convolution of solution kernels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of repeat convolutions
N = 42;

% values of constant alpha
for a = [0.50 0.49 0.45];
    
    % create the local neighborhood kernel
    h = [a 1-2*a a];
    
    % repeat convolution
    h2 = h;
    for I = 1:N; h2 = conv(h2, h); end
    
    % in frequency domain
    H2 = abs(fft(h2));
    
    % plot repeat convolution
    hold off
    N2 = length(h2);
    x2 = (1:N2)-(N2+1)/2;
    plot(x2, h2, 'k', 'LineWidth', 3)
    xlabel('Relative slice index')
    ylabel('\Psi')
    set(gca, 'FontSize', 18)
    axis([min(x2) max(x2) 0 max(h2)*1.1])
    text(length(x2)/2*.4, max(h2)*.9, ['\alpha = ' num2str(a)], 'FontSize', 18)
    
    saveas(gca, [FIGDIR filesep 'kernel-spatial-alpha-' num2str(a) '.tif'])
    
    % plot frequency domain
    M2 = length(H2);
    theta2 = linspace(0, 2, M2+1);
    theta2(end) = [];
    plot(theta2, H2, 'k', 'LineWidth', 3)
    xlabel('Angular frequency \Omega/ \pi')
    ylabel('|DFT(\Psi)|')
    set(gca, 'FontSize', 18)
    text(1.3, max(H2)*.9, ['\alpha = ' num2str(a)], 'FontSize', 18)
    
    saveas(gca, [FIGDIR filesep 'kernel-freq-alpha-' num2str(a) '.tif'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% visualization of transformation diffusion with vertical translation
%% (solution is sine wave)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coefficient
alpha = 0.49;

% number of slices
N = 100;

% place slices along true shape as sine wave
x = (0:(N-1))/N;
t = sin(2*pi*x);

% keep a copy of the true shape
t0 = t;

% random initialization of slice positions away from their perfect
% alignment
rng(10)
noise = rand(1, N);
noise = noise - mean(noise);
t = t + 2.4 * noise;

% plot slice positions
hold off
plot(0:N-1, t, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 16)
hold on
plot(0:N-1, t0, ':k', 'LineWidth', 2)
axis([0 100 -3 2])
xlabel('slice index')
ylabel('slice position')

% save figure
saveas(gca, [FIGDIR filesep 'registration-diffusion-toy-translation-ground-truth.png'])

% compute registration between slices
fprev = -diff(t);
fpos = -diff(t(end:-1:1));
fpos = fpos(end:-1:1);

% pad with values so that computing the slice update step is more simple.
% After padding, for the i-th slice we can do alpha*(fprev(i)+fpos(i))
fprev = [NaN fprev];
fpos = [fpos NaN];

% init accumulated transform to apply to each slice
ftot = zeros(1, N);

% NOTE: the optimal reconstruction is ITER = 42

% allocate space for reconstruction error
err = zeros(1, 7000);

%% reduce amount of registration noise, but don't get to optimal reconstruction
for I = 1:5

    % impose Neumann conditions at the extremes, i.e. 2*delta(0,1) and
    % 2*delta(N-1,N-2) respectively
    fprev(1) = fpos(1);
    fpos(end) = fprev(end);
    
    % compute the transform to apply to each slice in this iteration
    f = alpha * (fprev + fpos);
    
    % compose with the previous transforms to obtain a total transform
    ftot = ftot + f;
    
    % update neighbour transformation step (note that we ignore the element
    % for the ghost slices)
    fprev(2:end) = fprev(2:end) - f(2:end) + f(1:end-1);
    fpos(1:end-1) = fpos(1:end-1) - f(1:end-1) + f(2:end);
    
    % compute the error between the reconstruction and the true shape
    err(I) = norm(t + ftot - t0);
    
   
end

% plot slice positions
hold off
plot(0:N-1, t + ftot, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 16)
hold on
plot(0:N-1, t0, ':k', 'LineWidth', 2)
axis([0 100 -3 2])
xlabel('slice index')
ylabel('slice position')

% save figure
saveas(gca, [FIGDIR filesep 'registration-diffusion-toy-translation-iter-5.png'])

%% run more iterations until achieving optimal reconstruction
for I = 6:42

    % impose Neumann conditions at the extremes, i.e. 2*delta(0,1) and
    % 2*delta(N-1,N-2) respectively
    fprev(1) = fpos(1);
    fpos(end) = fprev(end);
    
    % compute the transform to apply to each slice in this iteration
    f = alpha * (fprev + fpos);
    
    % compose with the previous transforms to obtain a total transform
    ftot = ftot + f;
    
    % update neighbour transformation step (note that we ignore the element
    % for the ghost slices)
    fprev(2:end) = fprev(2:end) - f(2:end) + f(1:end-1);
    fpos(1:end-1) = fpos(1:end-1) - f(1:end-1) + f(2:end);
    
    % compute the error between the reconstruction and the true shape
    err(I) = norm(t + ftot - t0);
    
   
end

% plot slice positions
hold off
plot(0:N-1, t + ftot, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 16)
hold on
plot(0:N-1, t0, ':k', 'LineWidth', 2)
axis([0 100 -3 2])
xlabel('slice index')
ylabel('slice position')

% save figure
saveas(gca, [FIGDIR filesep 'registration-diffusion-toy-translation-iter-42.png'])

%% run more iterations until achieving optimal reconstruction
for I = 43:1000

    % impose Neumann conditions at the extremes, i.e. 2*delta(0,1) and
    % 2*delta(N-1,N-2) respectively
    fprev(1) = fpos(1);
    fpos(end) = fprev(end);
    
    % compute the transform to apply to each slice in this iteration
    f = alpha * (fprev + fpos);
    
    % compose with the previous transforms to obtain a total transform
    ftot = ftot + f;
    
    % update neighbour transformation step (note that we ignore the element
    % for the ghost slices)
    fprev(2:end) = fprev(2:end) - f(2:end) + f(1:end-1);
    fpos(1:end-1) = fpos(1:end-1) - f(1:end-1) + f(2:end);
    
    % compute the error between the reconstruction and the true shape
    err(I) = norm(t + ftot - t0);
    
   
end

% plot slice positions
hold off
plot(0:N-1, t + ftot, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 16)
hold on
plot(0:N-1, t0, ':k', 'LineWidth', 2)
axis([0 100 -3 2])
xlabel('slice index')
ylabel('slice position')

% save figure
saveas(gca, [FIGDIR filesep 'registration-diffusion-toy-translation-iter-1000.png'])

%% run more iterations until achieving optimal reconstruction
for I = 1001:7000

    % impose Neumann conditions at the extremes, i.e. 2*delta(0,1) and
    % 2*delta(N-1,N-2) respectively
    fprev(1) = fpos(1);
    fpos(end) = fprev(end);
    
    % compute the transform to apply to each slice in this iteration
    f = alpha * (fprev + fpos);
    
    % compose with the previous transforms to obtain a total transform
    ftot = ftot + f;
    
    % update neighbour transformation step (note that we ignore the element
    % for the ghost slices)
    fprev(2:end) = fprev(2:end) - f(2:end) + f(1:end-1);
    fpos(1:end-1) = fpos(1:end-1) - f(1:end-1) + f(2:end);
    
    % compute the error between the reconstruction and the true shape
    err(I) = norm(t + ftot - t0);
    
   
end

% plot slice positions
hold off
plot(0:N-1, t + ftot, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 16)
hold on
plot(0:N-1, t0, ':k', 'LineWidth', 2)
axis([0 100 -3 2])
xlabel('slice index')
ylabel('slice position')

% save figure
saveas(gca, [FIGDIR filesep 'registration-diffusion-toy-translation-iter-7000.png'])

%% plot reconstruction error
hold off
semilogx(err, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 16)
xlabel('diffusion iteration')
ylabel('error')

% save figure
saveas(gca, [FIGDIR filesep 'registration-diffusion-toy-translation-err.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% common code to all virtual slices sections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% list of blockface and histology files
[filebf, fileh] = match_file_lists(...
    dir([BFPDATADIR filesep '*.png']), ...
    dir([HISTLO2BFREGRIGIDDIR filesep '*.mha']), ...
    9:11, ...
    5:7 ...
    );
N = length(fileh);

% load histology slice, and change units to mm so that the plots look nicer
imh = scimat_load([HISTLO2BFREGRIGIDDIR filesep fileh(140).name]);
imh.axis(1).spacing = imh.axis(1).spacing * 1e3;
imh.axis(2).spacing = imh.axis(2).spacing * 1e3;

% coordinates of the planes we are going to cut through
R1 = 515;
R2 = 612;
C1 = 974;
y1 = scimat_index2world([R1 1], imh);
y1 = y1(2);
y2 = scimat_index2world([R2 1], imh);
y2 = y2(2);
x1 = scimat_index2world([1 C1], imh);
x1 = x1(1);

% coordinates of image limits
aux0 = scimat_index2world([1 1], imh);
auxend = scimat_index2world([size(imh.data, 1) size(imh.data, 2)], imh);
bbx = [aux0(1) auxend(1)]; % x-coordinates of the first and last corners of bounding box
bby = [aux0(2) auxend(2)]; % y-coordinates of the first and last corners of bounding box
bbz = [0 10e-6*1e3*(N-1)]; % z-coordinates of the first and last corners of bounding box

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% low-res virtual slices at different stages of refinement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% histology real slice example and reference

% load histology slice, and change units to mm so that the plots look nicer
imh = scimat_load([HISTLO2BFREGRIGIDDIR filesep fileh(140).name]);
imh.axis(1).spacing = imh.axis(1).spacing * 1e3;
imh.axis(2).spacing = imh.axis(2).spacing * 1e3;

% plot horizontal slice to use as reference
hold off
scimat_imagesc(imh)
set(gca, 'FontSize', 16)
hold on
rectangle('Position', [0.007945e3 0.0043e3 ...
    3579/20*imh.axis(2).spacing 4500/20*imh.axis(1).spacing], ...
    'LineStyle', '--', 'LineWidth', 2)
plot([0 13], [y1 y1], 'k', 'LineWidth', 4)
plot([0 13], [y2 y2], 'k', 'LineWidth', 4)
plot([x1 x1], [3 13], 'k', 'LineWidth', 4)
plot([0 13], [y1 y1], 'b', 'LineWidth', 2)
plot([0 13], [y2 y2], 'r', 'LineWidth', 2)
plot([x1 x1], [3 13], 'g', 'LineWidth', 2)
axis ij equal tight
axis([2 12 3 12])

% save image
export_fig([FIGDIR filesep MOUSE '_' HISTLO2BFREGRIGID '_slice_140.tif'])

%% real blockface example image

% load blockface
[~, namebf] = fileparts(filebf(140).name);
imbf = scimat_load([BFPDATADIR filesep namebf '.png']);
imbf.axis(1).spacing = imh.axis(1).spacing;
imbf.axis(2).spacing = imh.axis(2).spacing;
imbf.axis(1).min = imh.axis(1).min;
imbf.axis(2).min = imh.axis(2).min;

% plot horizontal slice to use as reference
hold off
scimat_imagesc(imbf)
set(gca, 'FontSize', 16)
set(gcf, 'color', 'w');
hold on
rectangle('Position', [0.007945e3 0.0043e3 ...
    3579/20*imh.axis(2).spacing 4500/20*imh.axis(1).spacing], ...
    'LineStyle', '--', 'LineWidth', 2, 'EdgeColor', 'k')
plot([0 13], [y1 y1], 'k', 'LineWidth', 4)
plot([0 13], [y2 y2], 'k', 'LineWidth', 4)
plot([x1 x1], [3 13], 'k', 'LineWidth', 4)
plot([0 13], [y1 y1], 'b', 'LineWidth', 2)
plot([0 13], [y2 y2], 'r', 'LineWidth', 2)
plot([x1 x1], [3 13], 'g', 'LineWidth', 2)
axis ij equal tight
axis([2 12 3 12])

% save image
export_fig([FIGDIR filesep MOUSE '_blockface_slice_140.tif'])

%% histology detail

% load high resolution slice and apply alignment to blockface so that we
% can plot the detail
load([IMPROCDIR filesep T_HISTLO2BF], 'th2bf')
opts.AutoDefaultPixelValue = true;
[~, nameh] = fileparts(fileh(140).name);
aux = transformix(th2bf(140), [HISTDIR filesep nameh '.tif'], opts);
imh = scimat_load(aux);
imh.axis(1).spacing = imh.axis(1).spacing * 1e3;
imh.axis(2).spacing = imh.axis(2).spacing * 1e3;

% plot histology detail
hold off
scimat_imagesc(imh)
set(gca, 'FontSize', 16)
hold on
rectangle('Position', [0.007945e3 0.0043e3 ...
    3579/20*imh.axis(2).spacing 4500/20*imh.axis(1).spacing], ...
    'LineStyle', '--', 'LineWidth', 2)
plot([0 13], [y1 y1], 'k', 'LineWidth', 4)
plot([0 13], [y2 y2], 'k', 'LineWidth', 4)
plot([x1 x1], [3 13], 'k', 'LineWidth', 4)
plot([0 13], [y1 y1], 'b', 'LineWidth', 2)
plot([0 13], [y2 y2], 'r', 'LineWidth', 2)
plot([x1 x1], [3 13], 'g', 'LineWidth', 2)
axis ij equal tight
axis([0.007945e3+[0 3579/20*imh.axis(2).spacing]+[-.05 .05]...
    0.0043e3+[0 4500/20*imh.axis(1).spacing]+[-.05 .05]])

% save image
export_fig([FIGDIR filesep MOUSE '_' HISTLO2BFREGRIGID '_slice_140_detail.tif'])

%% blockface virtual slice

clear imbfall

for I = 1:N
    
    % load blockface
    [~, namebf] = fileparts(filebf(I).name);
    imbfall(:, :, I, 1, :) = imread([BFPDATADIR filesep namebf '.png']);
    
end
imbfall = imbfall(:, :, end:-1:1, :, :);

% plot virtual slice 1
hold off
im = squeeze(imbfall(R1, :, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bbx, bbz, im)
axis xy equal tight
axis([3.5 10.5 0 1.1])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_R' num2str(R1) '_' MOUSE ...
    '_' 'Blockface' '.tif'])

% plot virtual slice 2
hold off
im = squeeze(imbfall(R2, :, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bbx, bbz, im)
axis xy equal tight
axis([3 11.5 0 1.5])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_R' num2str(R2) '_' MOUSE ...
    '_' 'Blockface' '.tif'])

% plot virtual slice 3
hold off
im = squeeze(imbfall(:, C1, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bby, bbz, im)
axis xy equal tight
axis([4 11.5 0 1.5])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_C' num2str(C1) '_' MOUSE ...
    '_' 'Blockface' '.tif'])


%% blockface virtual slice (detail)

% load the blockface using the previous section

% plot virtual slice 1
hold off
im = squeeze(imbfall(R1, :, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bbx, bbz, im)
axis xy equal tight
axis([7.9450 9.5790 0 0.9])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_R' num2str(R1) '_' MOUSE ...
    '_' 'Blockface_detail' '.tif'])

% plot virtual slice 2
hold off
im = squeeze(imbfall(R2, :, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bbx, bbz, im)
axis xy equal tight
axis([7.9450 9.5790 0 1.2])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_R' num2str(R2) '_' MOUSE ...
    '_' 'Blockface_detail' '.tif'])

% plot virtual slice 3
hold off
im = squeeze(imbfall(:, C1, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bby, bbz, im)
axis xy equal tight
axis([4.3000 6.3654 0 1.3])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_C' num2str(C1) '_' MOUSE ...
    '_' 'Blockface_detail' '.tif'])


%% histology alignment without blockface

% load and compose transforms
load([IMPROCDIR filesep T_NOBF_INTRAHISTLO], 'th2h', 'th2h_bsp')
ttot = elastix_cat(th2h_bsp, th2h);

% init memory to store transformed histology
imhall = zeros([ttot(1).Size([2 1]) N 1 3], 'uint8');

for I = 1:N
    
    % load histology
    [~, nameh] = fileparts(fileh(I).name);
    aux = scimat_load([HISTLODIR filesep nameh '.png']);
    
    % apply transforms to histology
    aux = transformix(ttot(I), aux, opts);
    imhall(:, :, I, 1, :) = aux.data;
    
end
imhall = imhall(:, :, end:-1:1, :, :);

% plot virtual slice 1
hold off
im = squeeze(imhall(285, :, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bbx, bbz, im)
axis xy equal tight
axis([3.25 11.5 0 1.5])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_R' num2str(285) '_' MOUSE ...
    '_' INTRANOBFHISTLORIGID '.tif'])

% plot virtual slice 2
hold off
im = squeeze(imhall(367, :, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bbx, bbz, im)
axis xy equal tight
axis([2.5 12.5 0 1.5])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_R' num2str(367) '_' MOUSE ...
    '_' INTRANOBFHISTLORIGID '.tif'])

% plot virtual slice 3
hold off
im = squeeze(imhall(:, 888, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bby, bbz, im)
axis xy equal tight
axis([2.5 11 0 1.5])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_C' num2str(888) '_' MOUSE ...
    '_' INTRANOBFHISTLORIGID '.tif'])

%% histology to blockface rigid registration

opts.AutoDefaultPixelValue = true;

% load and compose transforms
load([IMPROCDIR filesep T_HISTLO2BF], 'th2bf')
ttot = th2bf;

% init memory to store transformed histology
imhall = zeros([ttot(1).Size([2 1]) N 1 3], 'uint8');

for I = 1:N
    
    % load histology
    [~, nameh] = fileparts(fileh(I).name);
    aux = scimat_load([HISTLODIR filesep nameh '.png']);
    
    % apply transforms to histology
    aux = transformix(ttot(I), aux, opts);
    imhall(:, :, I, 1, :) = aux.data;
    
end
imhall = imhall(:, :, end:-1:1, :, :);

% plot virtual slice 1
hold off
im = squeeze(imhall(R1, :, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bbx, bbz, im)
axis xy equal tight
axis([3.5 10.5 0 1.1])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_R' num2str(R1) '_' MOUSE ...
    '_' HISTLO2BFREGRIGID '.tif'])

% plot virtual slice 2
hold off
im = squeeze(imhall(R2, :, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bbx, bbz, im)
axis xy equal tight
axis([3 11.5 0 1.5])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_R' num2str(R2) '_' MOUSE ...
    '_' HISTLO2BFREGRIGID '.tif'])

% plot virtual slice 3
hold off
im = squeeze(imhall(:, C1, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bby, bbz, im)
axis xy equal tight
axis([4 11.5 0 1.5])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_C' num2str(C1) '_' MOUSE ...
    '_' HISTLO2BFREGRIGID '.tif'])

%% intra-histology rigid diffusion

clear imhall
opts.AutoDefaultPixelValue = true;

% load and compose transforms up to an including this refinement
load([IMPROCDIR filesep T_HISTLO2BF], 'th2bf')
load([IMPROCDIR filesep T_INTRAHISTLO], 'th2h')
ttot = elastix_cat(th2h, th2bf);

% init memory to store transformed histology
imhall = zeros([ttot(1).Size([2 1]) N 1 3], 'uint8');

for I = 1:N
    
    % load histology
    [~, nameh] = fileparts(fileh(I).name);
    aux = scimat_load([HISTLODIR filesep nameh '.png']);
    
    % apply transforms to histology
    aux = transformix(ttot(I), aux, opts);
    imhall(:, :, I, 1, :) = aux.data;
    
end
imhall = imhall(:, :, end:-1:1, :, :);

% plot virtual slice 1
hold off
im = squeeze(imhall(R1, :, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bbx, bbz, im)
axis xy equal tight
axis([3.5 10.5 0 1.1])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_R' num2str(R1) '_' MOUSE ...
    '_' INTRAHISTLORIGID '.tif'])

% plot virtual slice 2
hold off
im = squeeze(imhall(R2, :, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bbx, bbz, im)
axis xy equal tight
axis([3 11.5 0 1.5])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_R' num2str(R2) '_' MOUSE ...
    '_' INTRAHISTLORIGID '.tif'])

% plot virtual slice 3
hold off
im = squeeze(imhall(:, C1, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bby, bbz, im)
axis xy equal tight
axis([4 11.5 0 1.5])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_C' num2str(C1) '_' MOUSE ...
    '_' INTRAHISTLORIGID '.tif'])

%% intra-histology B-spline diffusion

clear imh
opts.AutoDefaultPixelValue = true;

% compose transforms
load([IMPROCDIR filesep T_HISTLO2BF], 'th2bf')
load([IMPROCDIR filesep T_INTRAHISTLO], 'th2h', 'th2h_bsp')
ttot = elastix_cat(th2h_bsp, th2h, th2bf);

% init memory to store transformed histology
imhall = zeros([ttot(1).Size([2 1]) N 1 3], 'uint8');

for I = 1:N
    
    % load histology
    [~, nameh] = fileparts(fileh(I).name);
    aux = scimat_load([HISTLODIR filesep nameh '.tif']);
    
    % apply transforms to histology
    aux = transformix(ttot(I), aux, opts);
    imhall(:, :, I, 1, :) = aux.data;
    
end
imhall = imhall(:, :, end:-1:1, :, :);

% plot virtual slice 1
hold off
im = squeeze(imhall(R1, :, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bbx, bbz, im)
axis xy equal tight
axis([3.5 10.5 0 1.1])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_R' num2str(R1) '_' MOUSE ...
    '_' INTRAHISTLOBSP '.tif'])

% plot virtual slice 2
hold off
im = squeeze(imhall(R2, :, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bbx, bbz, im)
axis xy equal tight
axis([3 11.5 0 1.5])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_R' num2str(R2) '_' MOUSE ...
    '_' INTRAHISTLOBSP '.tif'])

% plot virtual slice 3
hold off
im = squeeze(imhall(:, C1, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bby, bbz, im)
axis xy equal tight
axis([4 11.5 0 1.5])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_C' num2str(C1) '_' MOUSE ...
    '_' INTRAHISTLOBSP '.tif'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% hi-res virtual slices at different stages of refinement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% compute the rows and column for the virtual slices

% load example image with low-res metainformation
imh = scimat_load([HISTLO2BFREGRIGIDDIR filesep fileh(140).name]);

aux = scimat_index2world([R1, 1], imh);
Y1 = aux(2);
aux = scimat_index2world([R2, 1], imh);
Y2 = aux(2);
aux = scimat_index2world([1, C1], imh);
X1 = aux(1);

% load example image with hi-res metainformation
imh = scimat_load([HISTRECDIR filesep fileh(140).name]);

aux = scimat_world2index([1, Y1], imh);
RH1 = round(aux(1));
aux = scimat_world2index([1, Y2], imh);
RH2 = round(aux(1));
aux = scimat_world2index([X1, 1], imh);
CH1 = round(aux(2));

% coordinates of image limits
aux0 = scimat_index2world([1 1], imh);
auxend = scimat_index2world([size(imh.data, 1) size(imh.data, 2)], imh);
bbx = [aux0(1) auxend(1)]*1e3; % x-coordinates of the first and last corners of bounding box
bby = [aux0(2) auxend(2)]*1e3; % y-coordinates of the first and last corners of bounding box
bbz = [0 10e-6*1e3*(N-1)]; % z-coordinates of the first and last corners of bounding box


%% histology cropped volume after lo-res reconstruction

% load total transforms: low-res reconstruction + cropping of the window of
% interest
load([HISTREFDIR_LVFREE filesep T_TOTHIST_LVFREE], 'ttotlo')

% init memory to store transformed histology
imhall = zeros([ttotlo(1).Size([2 1]) N 3], 'uint8');

for I = 1:N
    
    % load histology
    [~, nameh] = fileparts(fileh(I).name);
    aux = scimat_load([HISTRECDIR filesep fileh(I).name]);
    imhall(:, :, I, :, :) = aux.data;
    
end
imhall = imhall(:, :, end:-1:1, :, :);

% plot virtual slice 1
hold off
im = squeeze(imhall(RH1, :, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bbx, bbz, im)
axis xy equal tight
axis([bbx 0 0.9])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_R' num2str(RH1) '_' MOUSE ...
    '_' HISTREC_LVFREE '.tif'])

% plot virtual slice 2
hold off
im = squeeze(imhall(RH2, :, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bbx, bbz, im)
axis xy equal tight
axis([bbx 0 1.2])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_R' num2str(RH2) '_' MOUSE ...
    '_' HISTREC_LVFREE '.tif'])

% plot virtual slice 3
hold off
im = squeeze(imhall(:, CH1, :, :, :));
im = permute(im, [2 1 3]);
imagesc(bby, bbz, im)
axis xy equal tight
axis([bby 0 1.3])
set(gcf, 'color', 'white');
set(gca, 'FontSize', 16)

% save figure
export_fig([FIGDIR filesep 'virtual_slice_C' num2str(CH1) '_' MOUSE ...
    '_' HISTREC_LVFREE '.tif'])


%% apply all refinements to hi-res histology. 

% compose transforms: low-res reconstruction, cropping of window of
% interest, and first hi-res B-spline refinement
load([HISTREFDIR_LVFREE filesep T_TOTHIST_LVFREE], 'ttotlo', 'th2h_bsp')

% loop number of refinement sweeps
for NREF = 1:elastix_length(th2h_bsp)
    
    disp(['NREF = ' num2str(NREF)])
    
    % how many refinement sweeps do we want?
    t_aux = elastix_colon(th2h_bsp, ...
        5-NREF:elastix_length(th2h_bsp));
    
    % compose low-res refinement with the hi-res refinements
    ttot = elastix_cat(t_aux, ttotlo);
    
    % init memory to store transformed histology
    imhall = zeros([ttot(1).Size([2 1]) N 1 3], 'uint8');
    
    for I = 1:N
        
        disp(['I = ' num2str(I)])
        
        tic
        % apply transforms to histology
        [~, nameh] = fileparts(fileh(I).name);
        auxfile = transformix(ttot(I), [HISTDIR filesep nameh '.png'], opts);
        aux = scimat_load(auxfile);
        imhall(:, :, I, 1, :) = aux.data;
        delete(auxfile)
        disp(['time: ' num2str(toc) ' sec'])
        
    end
    imhall = imhall(:, :, end:-1:1, :, :);
    
    save(['foo' num2str(NREF) '.mat'], 'imhall', '-v7.3')
    
    % plot virtual slice 1
    hold off
    im = squeeze(imhall(RH1, :, :, :, :));
    im = permute(im, [2 1 3]);
    imagesc(bbx, bbz, im)
    axis xy equal tight
    axis([bbx 0 0.9])
    set(gcf, 'color', 'white');
    set(gca, 'FontSize', 16)
    
    % save figure
    export_fig([FIGDIR filesep 'virtual_slice_R' num2str(RH1) '_' MOUSE ...
        '_' HISTREF_LVFREE '_nref_' num2str(NREF) '.tif'])
    
    % plot virtual slice 2
    hold off
    im = squeeze(imhall(RH2, :, :, :, :));
    im = permute(im, [2 1 3]);
    imagesc(bbx, bbz, im)
    axis xy equal tight
    axis([bbx 0 1.3])
    set(gcf, 'color', 'white');
    set(gca, 'FontSize', 16)
    
    % save figure
    export_fig([FIGDIR filesep 'virtual_slice_R' num2str(RH2) '_' MOUSE ...
        '_' HISTREF_LVFREE '_nref_' num2str(NREF) '.tif'])
    
    % plot virtual slice 3
    hold off
    im = squeeze(imhall(:, CH1, :, :, :));
    im = permute(im, [2 1 3]);
    imagesc(bby, bbz, im)
    axis xy equal tight
    axis([bby 0 1.3])
    set(gcf, 'color', 'white');
    set(gca, 'FontSize', 16)
    
    % save figure
    export_fig([FIGDIR filesep 'virtual_slice_C' num2str(CH1) '_' MOUSE ...
        '_' HISTREF_LVFREE '_nref_' num2str(NREF) '.tif'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reconstruction validation (no blockface as external reference)
%% slice by slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load error results:

% low-res 4 B-spline registrations, 5 diffusions each
err_h2bf_nobf = load([SRCDIR filesep 'err_Q53_blockface_histology_t_NoBlockfaceIntra_HistologySiriusRedLoRes.mat']);
err_h2h_nobf = load([SRCDIR filesep 'err_Q53_intra_histology_t_NoBlockfaceIntra_HistologySiriusRedLoRes.mat']);

% tidy up variables
err_h2bf_nobf = err_h2bf_nobf.err;
err_h2h_nobf = err_h2h_nobf.err;

% remove wrong landmark due to human error when hand tracing them
err_h2bf_nobf{126}(13, :) = [];
err_h2h_nobf{11}(end, :) = [];

%% plot errors, but each slice becomes a point in the plot. This way we can
%% see that the error increases as we move left away from the slice we
%% started from

%% histo-blockface

% indices of slices with landmarks
idx = find(~cellfun(@isempty, err_h2bf_nobf));

close
hold off
J = 5; % coarsest level of refinement
plot(idx, ...
    cellfun(@(x) median(x(:, J), 1), err_h2bf_nobf(idx)) * 1e3, ...
    'k-', 'LineWidth', 2)
hold on
J = 1; % finest level of refinement
plot(idx, ...
    cellfun(@(x) median(x(:, J), 1), err_h2bf_nobf(idx)) * 1e3, ...
    'k--', 'LineWidth', 2)
J = 5;
plot(idx, ...
    cellfun(@(x) prctile(x(:, J), 100, 1), err_h2bf_nobf(idx)) * 1e3, ...
    'r', 'LineWidth', 2)
J = 1;
plot(idx, ...
    cellfun(@(x) prctile(x(:, J), 100, 1), err_h2bf_nobf(idx)) * 1e3, ...
    'r--', 'LineWidth', 2)
set(gca, 'FontSize', 18)
xlabel('Slice index')
ylabel('Landmark distance error (mm)')
text(75, 2.1, 'direction of registration', 'FontSize', 18)
ah = annotation('textarrow', [.85 .5], [.8 .8]);
axis([0 150 0 2.4])

saveas(gca, [FIGDIR filesep 'err_Q53_blockface_histology_t_NoBlockfaceIntra_HistologySiriusRedLoRes_by_slice.tif'])

%% histo-histo

% indices of slices with landmarks
idx = find(~cellfun(@isempty, err_h2h_nobf));

close
hold off
J = 5;
plot(idx, ...
    cellfun(@(x) median(x(:, J), 1), err_h2h_nobf(idx)) * 1e3, ...
    'k', 'LineWidth', 2)
hold on
J = 1;
plot(idx, ...
    cellfun(@(x) median(x(:, J), 1), err_h2h_nobf(idx)) * 1e3, ...
    'k--', 'LineWidth', 2)
J = 5;
plot(idx, ...
    cellfun(@(x) prctile(x(:, J), 100, 1), err_h2h_nobf(idx)) * 1e3, ...
    'r', 'LineWidth', 2)
J = 1;
plot(idx, ...
    cellfun(@(x) prctile(x(:, J), 100, 1), err_h2h_nobf(idx)) * 1e3, ...
    'r--', 'LineWidth', 2)
set(gca, 'FontSize', 18)
xlabel('Slice index')
ylabel('Landmark distance error (mm)')
text(75, 0.4, 'direction of registration', 'FontSize', 18)
ah = annotation('textarrow', [.85 .5], [.65 .65]);
axis([0 150 0 0.575])

saveas(gca, [FIGDIR filesep 'err_Q53_intra_histology_t_NoBlockfaceIntra_HistologySiriusRedLoRes_by_slice.tif'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reconstruction validation (using blockface as external reference)
%% slice by slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% low-res 4 B-spline registrations, 5 diffusions
err_h2bf = load([SRCDIR filesep 'err_Q53_blockface_histology_t_Intra_HistologySiriusRedLoRes.mat']);
err_h2h = load([SRCDIR filesep 'err_Q53_intra_histology_t_Intra_HistologySiriusRedLoRes.mat']);

% tidy up variables
err_h2bf = err_h2bf.err;
err_h2h = err_h2h.err;

% remove wrong landmark due to human error when hand tracing them
err_h2bf{126}(13, :) = [];
err_h2h{11}(end, :) = [];

%% plot errors, but each slice becomes a point in the plot. This way we can
%% see that the error increases as we move left away from the slice we
%% started from

%% histo-blockface

% indices of slices with landmarks
idx = find(~cellfun(@isempty, err_h2bf));

close
hold off
J = 5; % coarsest level of refinement
plot(idx, ...
    cellfun(@(x) median(x(:, J), 1), err_h2bf(idx)) * 1e3, ...
    'k-', 'LineWidth', 2)
hold on
J = 1; % finest level of refinement
plot(idx, ...
    cellfun(@(x) median(x(:, J), 1), err_h2bf(idx)) * 1e3, ...
    'k--', 'LineWidth', 2)
J = 5;
plot(idx, ...
    cellfun(@(x) prctile(x(:, J), 100, 1), err_h2bf(idx)) * 1e3, ...
    'r', 'LineWidth', 2)
J = 1;
plot(idx, ...
    cellfun(@(x) prctile(x(:, J), 100, 1), err_h2bf(idx)) * 1e3, ...
    'r--', 'LineWidth', 2)
set(gca, 'FontSize', 18)
xlabel('Slice index')
ylabel('Landmark distance error (mm)')
axis([0 150 0 2.4])

saveas(gca, [FIGDIR filesep 'err_Q53_blockface_histology_t_Intra_HistologySiriusRedLoRes_by_slice.tif'])

%% histo-histo

% indices of slices with landmarks
idx = find(~cellfun(@isempty, err_h2h));

close
hold off
J = 5;
plot(idx, ...
    cellfun(@(x) median(x(:, J), 1), err_h2h(idx)) * 1e3, ...
    'k', 'LineWidth', 2)
hold on
J = 1;
plot(idx, ...
    cellfun(@(x) median(x(:, J), 1), err_h2h(idx)) * 1e3, ...
    'k--', 'LineWidth', 2)
J = 5;
plot(idx, ...
    cellfun(@(x) prctile(x(:, J), 100, 1), err_h2h(idx)) * 1e3, ...
    'r', 'LineWidth', 2)
J = 1;
plot(idx, ...
    cellfun(@(x) prctile(x(:, J), 100, 1), err_h2h(idx)) * 1e3, ...
    'r--', 'LineWidth', 2)
set(gca, 'FontSize', 18)
xlabel('Slice index')
ylabel('Landmark distance error (mm)')
axis([0 150 0 0.575])

saveas(gca, [FIGDIR filesep 'err_Q53_intra_histology_t_Intra_HistologySiriusRedLoRes_by_slice.tif'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reconstruction validation (using blockface as external reference)
%% aggregate slices, show each step of refinement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load error results:

% low-res 4 B-spline registrations, 5 diffusions
err_h2bf = load([SRCDIR filesep 'err_Q53_blockface_histology_t_Intra_HistologySiriusRedLoRes.mat']);
err_h2h = load([SRCDIR filesep 'err_Q53_intra_histology_t_Intra_HistologySiriusRedLoRes.mat']);

% low-res 20 B-spline registrations, 1 diffusion (i.e. no diffusion)
err_h2bf_nodiff = load([SRCDIR filesep 'err_Q53_blockface_histology_t_IntraNoDiff_HistologySiriusRedLoRes.mat']);
err_h2h_nodiff = load([SRCDIR filesep 'err_Q53_intra_histology_t_IntraNoDiff_HistologySiriusRedLoRes.mat']);

% low-res 1 B-spline registration, several levels of diffusion
err_h2bf_onlydiff = load([SRCDIR filesep 'err_Q53_blockface_histology_t_IntraOnlyDiff_HistologySiriusRedLoRes.mat']);
err_h2h_onlydiff = load([SRCDIR filesep 'err_Q53_intra_histology_t_IntraOnlyDiff_HistologySiriusRedLoRes.mat']);

% % hi-res 4 B-spline registrations, 100 diffusions each
% err_h2bf_hires = load([SRCDIR filesep 'err_Q53_blockface_histology_t_Total_HistologySiriusRed_LV_FreeWall.mat']);
% err_h2h_hires = load([SRCDIR filesep 'err_Q53_intra_histology_t_Total_HistologySiriusRed_LV_FreeWall.mat']);

% tidy up variables
MaxDiffIter_onlydiff = err_h2h_onlydiff.MaxDiffIter;
err_h2bf = err_h2bf.err;
err_h2bf_nodiff = err_h2bf_nodiff.err;
err_h2bf_onlydiff = err_h2bf_onlydiff.err;
err_h2h = err_h2h.err;
err_h2h_nodiff = err_h2h_nodiff.err;
err_h2h_onlydiff = err_h2h_onlydiff.err;
% err_h2bf_hires = err_h2bf_hires.err;
% err_h2h_hires = err_h2h_hires.err;

% concatenate error by landmarks so that we can see how it changes from
% level to level of refinement. 
errtot_h2bf = cat(1, err_h2bf{:});
errtot_h2bf_nodiff = cat(1, err_h2bf_nodiff{:});
errtot_h2bf_onlydiff = cat(1, err_h2bf_onlydiff{:});
errtot_h2h = cat(1, err_h2h{:});
errtot_h2h_nodiff = cat(1, err_h2h_nodiff{:});
errtot_h2h_onlydiff = cat(1, err_h2h_onlydiff{:});
% errtot_h2bf_hires = cat(1, err_h2bf_hires{:});
% errtot_h2h_hires = cat(1, err_h2h_hires{:});

% remove outlier, wrong placement of landmarks by user
errtot_h2bf(327, :) = [];
errtot_h2bf_nodiff(327, :) = [];
errtot_h2bf_onlydiff(327, :) = [];
errtot_h2h(39, :) = [];
errtot_h2h_nodiff(39, :) = [];
errtot_h2h_onlydiff(39, :) = [];
% errtot_h2bf_hires(327, :) = [];
% errtot_h2h_hires(39, :) = [];

% flip error left to right, so that we start with the coarsest level of
% reconstruction, and move towards finer reconstruction
errtot_h2bf(:, end:-1:1) = errtot_h2bf;
errtot_h2bf_nodiff(:, end:-1:1) = errtot_h2bf_nodiff;
errtot_h2h(:, end:-1:1) = errtot_h2h;
errtot_h2h_nodiff(:, end:-1:1) = errtot_h2h_nodiff;
% errtot_h2bf_hires(:, end:-1:1) = errtot_h2bf_hires;
% errtot_h2h_hires(:, end:-1:1) = errtot_h2h_hires;

% % the high resolution errors have lots of NaN because most of the landmarks
% % fall outside the ROI. We have to remove those NaNs
% errtot_h2bf_hires(any(isnan(errtot_h2bf_hires), 2), :) = [];
% errtot_h2h_hires(any(isnan(errtot_h2h_hires), 2), :) = [];

%% plot error

% plot lo-res blockface-histology landmark error (curves)
hold off
plot([1 2 7+(0:5:15)], median(errtot_h2bf, 1) * 1e3, 'k', 'LineWidth', 2)
hold on
plot([1 2 7+(0:5:15)], median(errtot_h2bf, 1) * 1e3, 'ok', 'LineWidth', 2)
set(gca, 'XTick', 1:22)
set(gca, 'XTickLabel', ...
    {'Bf', 'Rig', '1', '', '', '', '5', '', '', '', '', '10', '', '', '', '', '15', '', '', '', '', '20'})
set(gca, 'FontSize', 16)
rotateXLabels(gca, 45)
plot(2:22, median(errtot_h2bf_nodiff(:, 2:end), 1) * 1e3, 'k--', 'LineWidth', 2)
plot(2+MaxDiffIter_onlydiff(1:2), median(errtot_h2bf_onlydiff(:, 1:2) * 1e3, 1), 'k:', 'LineWidth', 2)
plot(2+MaxDiffIter_onlydiff(1:2), median(errtot_h2bf_onlydiff(:, 1:2) * 1e3, 1), 'kx', 'LineWidth', 2)
plot([1 2 7+(0:5:15)], prctile(errtot_h2bf, 95, 1) * 1e3, 'r', 'LineWidth', 2)
plot([1 2 7+(0:5:15)], prctile(errtot_h2bf, 95, 1) * 1e3, 'ro', 'LineWidth', 2)
plot(2:22, prctile(errtot_h2bf_nodiff(:, 2:end), 95, 1) * 1e3, 'r--', 'LineWidth', 2)
plot(2+MaxDiffIter_onlydiff(1:2), prctile(errtot_h2bf_onlydiff(:, 1:2) * 1e3, 95, 1), 'r:', 'LineWidth', 2)
plot(2+MaxDiffIter_onlydiff(1:2), prctile(errtot_h2bf_onlydiff(:, 1:2) * 1e3, 95, 1), 'rx', 'LineWidth', 2)
xlabel('Number of registration sweeps or equivalent neighbor transformation updates')
ylabel('Landmark distance error (mm)')

% save figure
saveas(gca, [FIGDIR filesep 'errtot_h2bf_' MOUSE '-' HISTO '.tif'])

% plot lo-res histology-histology landmark error
hold off
plot([1 2 7+(0:5:15)], median(errtot_h2h, 1) * 1e3, 'k', 'LineWidth', 2)
hold on
plot([1 2 7+(0:5:15)], median(errtot_h2h, 1) * 1e3, 'ok', 'LineWidth', 2)
set(gca, 'XTick', 1:22)
set(gca, 'XTickLabel', ...
    {'Bf', 'Rig', '1', '', '', '', '5', '', '', '', '', '10', '', '', '', '', '15', '', '', '', '', '20'})
set(gca, 'FontSize', 16)
rotateXLabels(gca, 45)
plot(2:22, median(errtot_h2h_nodiff(:, 2:end), 1) * 1e3, 'k--', 'LineWidth', 2)
plot(2+MaxDiffIter_onlydiff(1:2), median(errtot_h2h_onlydiff(:, 1:2) * 1e3, 1), 'k:', 'LineWidth', 2)
plot(2+MaxDiffIter_onlydiff(1:2), median(errtot_h2h_onlydiff(:, 1:2) * 1e3, 1), 'kx', 'LineWidth', 2)
plot([1 2 7+(0:5:15)], prctile(errtot_h2h, 95, 1) * 1e3, 'r', 'LineWidth', 2)
plot([1 2 7+(0:5:15)], prctile(errtot_h2h, 95, 1) * 1e3, 'ro', 'LineWidth', 2)
plot(2:22, prctile(errtot_h2h_nodiff(:, 2:end), 95, 1) * 1e3, 'r--', 'LineWidth', 2)
plot(2+MaxDiffIter_onlydiff(1:2), prctile(errtot_h2h_onlydiff(:, 1:2) * 1e3, 95, 1), 'r:', 'LineWidth', 2)
plot(2+MaxDiffIter_onlydiff(1:2), prctile(errtot_h2h_onlydiff(:, 1:2) * 1e3, 95, 1), 'rx', 'LineWidth', 2)
xlabel('Number of registration sweeps or equivalent neighbor transformation updates')
ylabel('Landmark distance error (mm)')

% save figure
saveas(gca, [FIGDIR filesep 'errtot_h2h_' MOUSE '-' HISTO '.tif'])

