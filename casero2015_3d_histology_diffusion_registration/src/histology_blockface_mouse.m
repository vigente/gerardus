% histology_blockface_mouse.m
%
% Script to reconstruct one half of Mouse Q53 heart's low-res histology,
% using blockface as external reference, and detail of Left Ventricle (LV)
% free wall.
%
%   - Correct alignment of blockface stack.
%   - Pre-process blockface images.
%   - Rigid alignment of low-res histology to blockface.
%   - Rigid diffusion refinement of low-res intra-histology.
%   - B-spline diffusion refinement of low-res intra-histology.
%     (4 registration sweeps, 5 neighbor updates each).
%   - B-spline diffusion refinement of hi-res intra-histology, 
%     detail of LV free wall (4 registration sweeps, 
%     100 neighbor updates each).

% version: 0.4.7
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

DOCDIR = '/home/orie1416/Software/casero2015_3d_histology_diffusion_registration/doc';
SRCDIR = [DOCDIR '/../src'];
FIGDIR = [DOCDIR '/figures'];

MOUSE = 'Q53';

IMPROCDIR = ['/data2/Mouse/' MOUSE '/ImageProcessing'];
BFDATADIR = ['/data2/Mouse/' MOUSE '/Blockface'];
BFSDATADIR = [IMPROCDIR '/BlockfaceStabilised'];
BFCDATADIR = [IMPROCDIR '/BlockfaceCorrected'];
BFPDATADIR = [IMPROCDIR '/BlockfacePreprocessed'];
% HISTO = 'HistologyTrichrome';
HISTO = 'HistologySiriusRed';
HISTREC = [HISTO 'Reconstructed'];
HISTOLO = [HISTO 'LoRes'];
HISTLODIR = ['/data2/Mouse/' MOUSE filesep HISTOLO];
HISTDIR = ['/data2/Mouse/' MOUSE filesep HISTO];
HISTRECDIR = [IMPROCDIR filesep HISTREC];
HISTLOREGDIR = [IMPROCDIR filesep HISTOLO '_to_BlockfaceCorrectedRegistration'];
HISTLO2BFREGRIGID = [HISTOLO '_to_BlockfaceRegistrationRigid'];
HISTLO2BFREGRIGIDDIR = [IMPROCDIR filesep HISTLO2BFREGRIGID];
INTRAHISTLORIGID = ['Intra_' HISTOLO];
INTRAHISTLOBSP = ['Intra_' HISTOLO '_Bsp'];

HISTLOPDIR = [HISTLO2BFREGRIGIDDIR 'Preprocessed'];
HISTPDIR = [HISTRECDIR 'Preprocessed'];

HISTRECDIR_LVFREE = [HISTRECDIR '_LV_FreeWall'];
HISTREFDIR_LVFREE = [HISTRECDIR 'Refined_LV_FreeWall'];

T_HISTLO2BF = ['t_' HISTOLO '_to_Blockface.mat'];
T_INTRAHISTLO = ['t_Intra_' HISTOLO '.mat'];
T_TOTHIST_LVFREE = ['t_Total_' HISTO '_LV_FreeWall.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Blockface intraframe registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find intraframe transformations
tic
[t55, tParam55, iterInfo55, regParam55] ...
    = blockface_intraframe_reg(BFDATADIR, dir([BFDATADIR '/*_55_*.png']));
toc
tic
[t90, tParam90, iterInfo90, regParam90] ...
    = blockface_intraframe_reg(BFDATADIR, dir([BFDATADIR '/*_90_*.png']));
toc

% save transforms
save([IMPROCDIR filesep 'Blockface55_IntraframeReg.mat'], ...
    't55', 'tParam55', 'iterInfo55', 'regParam55');
save([IMPROCDIR filesep 'Blockface90_IntraframeReg.mat'], ...
    't90', 'tParam90', 'iterInfo90', 'regParam90');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Eyeballing confirmation of camera shift frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load transforms
load([IMPROCDIR filesep 'Blockface55_IntraframeReg.mat'], ...
    't55', 'tParam55', 'iterInfo55', 'regParam55');
load([IMPROCDIR filesep 'Blockface90_IntraframeReg.mat'], ...
    't90', 'tParam90', 'iterInfo90', 'regParam90');

switch MOUSE

    case 'Q53'
        
        % shifts that propagate to subsequent frames
        idx55prop = [68 94 176 238 285];
        idx90prop = [202 221 222 231 232 241 242];
        
        % shifts that happen in an isolated frame, and don't need to
        % propagate
        idx55noprop = [221 231 241];
        idx90noprop = [];
    
    case 'Q62'
        
        % shifts that propagate to subsequent frames
        idx55prop = [3 14 15 18 20 22 33 45 58 59 75 76 82 83 85 86 ...
            89 90 96 97 106 107 108 109 137 141 144 145 151 152 157 ...
            158 189 190 195 197 198 205 206 219 223 224];
        idx90prop = [14 15 18 20 22 23 25 26 33 43 44 45 46 47 48 58 59 ...
            64 66 68 69 70 71 75 76 82 83 84 85 86 89 90 96 97 102 103 ...
            104 106 107 108 109 131 132 137 141 144 145 148 149 150 151 ...
            152 153 154 155 158 168 169 170 171 172 189 190 192:195 197 ...
            198 200:202 204:208 210:212 223 224 234 235];
        
        % shifts that happen in an isolated frame, and don't need to
        % propagate
        idx55noprop = [];
        idx90noprop = [];
    
end

% save indices of frames to correct
save([IMPROCDIR filesep 'Blockface55_IntraframeReg.mat'], ...
    'idx55prop', 'idx55noprop', '-append')
save([IMPROCDIR filesep 'Blockface90_IntraframeReg.mat'], ...
    'idx90prop', 'idx90noprop', '-append')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stabilise frame shifts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load translations
load([IMPROCDIR filesep 'Blockface55_IntraframeReg.mat'], ...
    't55', 'tParam55', 'iterInfo55', 'regParam55', 'idx55prop', 'idx55noprop');
load([IMPROCDIR filesep 'Blockface90_IntraframeReg.mat'], ...
    't90', 'tParam90', 'iterInfo90', 'regParam90', 'idx90prop', 'idx90noprop');

% correct frame shifts
blockface_correct_frame_shifts(BFDATADIR, ...
    dir([BFDATADIR '/*_55_*.png']), tParam55, idx55prop, idx55noprop, ...
    BFSDATADIR);
blockface_correct_frame_shifts(BFDATADIR, ...
    dir([BFDATADIR '/*_90_*.png']), tParam90, idx90prop, idx90noprop, ...
    BFSDATADIR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find corresponding points to correct perspective of 55º frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% list files 
files55 = dir([BFSDATADIR filesep '*_55_*.png']);
files90 = dir([BFSDATADIR filesep '*_90_*.png']);
if (length(files55) ~= length(files90))
    error('Different number of frames for 55º and 90º images')
end

% pick up a few of the frames to find the map from 55º to 90º
pts55 = cell(1, length(files55));
pts90 = cell(1, length(files90));
for I = 1:60:length(files55)
    
    % load images
    im55 = imread([BFSDATADIR filesep files55(I).name]);
    im90 = imread([BFSDATADIR filesep files90(I).name]);

    % extend histogram to dynamic range for better visualization
    min55 = double(min(im55(:)));
    max55 = double(max(im55(:)));
    im55 = uint8((double(im55) - min55) / (max55 - min55) * 255);
    min90 = double(min(im90(:)));
    max90 = double(max(im90(:)));
    im90 = uint8((double(im90) - min90) / (max90 - min90) * 255);
    
    % find matched landmarks in both images
    [pts55{I}, pts90{I}] = cpselect(im55, im90, 'Wait', true);
    
end

% compute projective transformation to correct perspective
tformPerspective = fitgeotrans(cat(1, pts55{:}), cat(1, pts90{:}), 'Projective');

% save matched landmarks and projective transform
save([IMPROCDIR filesep 'BlockfaceStabilisedLandmarks_55_90.mat'], ...
    'pts55', 'pts90', 'tformPerspective')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correct perspective of 45º frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load projective transform
load([IMPROCDIR filesep 'BlockfaceStabilisedLandmarks_55_90.mat'], ...
    'tformPerspective')

% list files 
files55 = dir([BFSDATADIR filesep '*_55_*.png']);

% perspective correction
blockface_correct_perspective(BFSDATADIR, files55, tformPerspective, ...
    BFCDATADIR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find angle of wax block scratches by eyeballing scratches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% list of blockface images
files = dir([BFCDATADIR filesep '*_55_*.png']);

% click a pair of points on the same scratch, and repeat for several
% scratches

x = cell(1, length(files));
y = cell(1, length(files));
for I = 1:60:length(files)
    
    % load blockface image
    im = imread([BFCDATADIR filesep files(I).name]);
    
    % convert to colourmap
    imagesc(im)
    
    % click pairs of points on a few scratchs
    [x{I}, y{I}] = getpts;

end

% concatenate all the points that have been clicked on scratches
pts = [cat(1, x{:}) cat(1, y{:})];

% split into points on the left and points on the right
ptsl = pts(1:2:end, :);
ptsr = pts(2:2:end, :);

% angle of the scratches (rad)
alpha = atan2(ptsr(:, 2)-ptsl(:, 2), ptsr(:, 1)-ptsl(:, 1));

% median angle of the scratches
alphamed = median(alpha);

% load projective transform
load([IMPROCDIR filesep 'BlockfaceStabilisedLandmarks_55_90.mat'], ...
    'tformPerspective');

% center of the heart in X, Y coordinates
c = [1344 831];

% create tform to correct the scratch angulation, so that scratches become
% horizontal
R = [cos(-alphamed) sin(-alphamed); ...
    -sin(-alphamed) cos(-alphamed)];
tformScratch = affine2d([...
    R [0; 0];
    c*(eye(2)-R) 1]);
    
% save scratch landmarks and rotation transform
save([IMPROCDIR filesep 'BlockfaceScratchAngle_55_90.mat'], ...
    'x', 'y', 'alphamed', 'tformScratch')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correct perspective 55º->90º and scratch angle in one go
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load projective transform
load([IMPROCDIR filesep 'BlockfaceStabilisedLandmarks_55_90.mat'], ...
    'tformPerspective')

% load scratch angle transform
load([IMPROCDIR filesep 'BlockfaceScratchAngle_55_90.mat'], ...
    'tformScratch')

% create combined correction
tformCombined = projective2d(tformPerspective.T * tformScratch.T);
tformCombined = projective2d(tformScratch.T * tformPerspective.T);
tformCombined = tformPerspective;
tformCombined.T(1:2, 1:2) = tformCombined.T(1:2, 1:2) * tformPerspective.T(1:2, 1:2)

% list files 
files55 = dir([BFSDATADIR filesep '*_55_*.png']);

% perspective and scratch angle correction
blockface_correct_perspective(BFSDATADIR, files55, tformCombined, ...
    BFCDATADIR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create masks to correct illumination of 55º->90º frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create blockface masks

% create masks for the blockface volume
[ellipmask, polymask] = blockface_create_masks(BFCDATADIR, ...
    dir([BFCDATADIR filesep '*_55_*.png']));

% save the masks
save([IMPROCDIR filesep 'BlockfaceMasks.mat'], 'ellipmask', 'polymask')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rigid registration of histology to blockface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% list of blockface and histology files
[filebf, fileh] = match_file_lists(...
    dir([BFCDATADIR filesep '*.png']), ...
    dir([HISTLODIR filesep '*.png']), ...
    9:11, ...
    5:7 ...
    );
N = length(filebf);

%% cropping and scaling of blockface

% crop factors for blockface images
crop_row = 284:1254;
crop_col = 766:1864;

% size scaling = histology / blockface
% this is estimated by measuring the same antomical distance in a blockface
% image and the corresponding histology image
sc = 866.5806 / 594.4855;

% load mask to equalise blockface illumination
load([IMPROCDIR filesep 'BlockfaceMasks.mat'], 'ellipmask', 'polymask')

% load blockface images
info = imfinfo([BFCDATADIR filesep filebf(1).name]);
imbf = zeros(info.Height, info.Width, N, ['uint' num2str(info.BitDepth)]);
for I = 1:N
    imbf(:, :, I) = imread([BFCDATADIR filesep filebf(I).name]);
end

% crop blockface data
ellipmask = ellipmask(crop_row, crop_col);
polymask = polymask(crop_row, crop_col);
imbf = imbf(crop_row, crop_col, :);

% upscale masks so that they are in the same scale as histology
ellipmask = imresize(ellipmask, sc);
polymask = imresize(polymask, sc);

% upscale blockface so that it's in the same scale as histology
imbf = imresize(imbf, sc);

% save masks
save([BFPDATADIR filesep 'BlockfaceMasks.mat'], 'ellipmask', 'polymask')

%% pre-processing of blockface images

% correction algorithm parameters
ratio = 1/16;
thr = 45;
radheart = 10;
radpoly = 50;

for I = 1:N
    disp(['I = ' num2str(I)])
    
    % load histology so that we have pixel size values
    scimath = scimat_load([HISTLODIR filesep fileh(I).name]);

    % equalise blockface image
    scimath.data = blockface_equalise_illumination(...
        imbf(:, :, I), polymask, ellipmask, ratio, thr, radheart, radpoly);
    scimath.axis(1).size = size(scimath.data, 1);
    scimath.axis(2).size = size(scimath.data, 2);
    scimath.axis(1).min = -scimath.axis(1).spacing / 2;
    scimath.axis(2).min = -scimath.axis(2).spacing / 2;

    % save blockface
    [~, auxfilename] = fileparts(filebf(I).name);
    scimat_save([BFPDATADIR filesep auxfilename '.mha'], scimath);
end

%% matched filter and rigid registration of histology to blockface

% list of blockface and histology files
[filebf, fileh] = match_file_lists(...
    dir([BFPDATADIR filesep '*.mha']), ...
    dir([HISTLODIR filesep '*.png']), ...
    9:11, ...
    5:7 ...
    );
N = length(filebf);

% iterate blockface images
delete(gcp)
mypool = parpool('local', 3);
parfor I = 1:N

    disp(['I = ' num2str(I) ', Histology = ' fileh(I).name])
    th2bf(I) = aux_histology2blockface([HISTLODIR filesep fileh(I).name], ...
        [BFPDATADIR filesep filebf(I).name], HISTLO2BFREGRIGIDDIR, ...
        polymask, ellipmask);

end

% save transforms
save([IMPROCDIR filesep T_HISTLO2BF], 'th2bf')

%% visually check histology to blockface alignment

% list of blockface and histology files
[filebf, fileh] = match_file_lists(...
    dir([BFPDATADIR filesep '*.mha']), ...
    dir([HISTLO2BFREGRIGIDDIR filesep '*.mha']), ...
    9:11, ...
    5:7 ...
    );
N = length(filebf);

aux = scimat_load([HISTLO2BFREGRIGIDDIR filesep fileh(1).name]);
imh = aux.data;
imh = zeros([size(imh) N], 'like', imh);

for I = 1:N
    
    % load blockface and registered histology
    imbf = scimat_load([BFPDATADIR filesep filebf(I).name]);
    aux = scimat_load([HISTLO2BFREGRIGIDDIR filesep fileh(I).name]);
    imh(:, :, :, :, :, I) = aux.data;
    
    % plot alignment
    subplot(2, 2, 1)
    imagesc(imbf.data)
    subplot(2, 2, 2)
    imagesc(squeeze(imh(:, :, :, :, :, I)))
    subplot(2, 1, 2)
    imagesc(imfuse(imbf.data, squeeze(imh(:, :, :, :, :, I))))
%     pause
    
end

% plot cross horizontal section
subplot(2, 1, 1)
hold off 
imagesc(permute(squeeze(imh(974, :, :, :, :, :)), [3 1 2]))
set(gca, 'FontSize', 16)

% plot cross vertical section
subplot(2, 1, 2)
hold off
imagesc(permute(squeeze(imh(:, 967, :, :, :, :)), [3 1 2]))
set(gca, 'FontSize', 16)

saveas(gca, [FIGDIR filesep HISTLO2BFREGRIGID '-cross-sections.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Intra-histology registration correction: rigid diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% list of histology slices that are good enough for intra-histology
% registration
idxintra = [12:32 35:41 43:55 59 61 63 66:72 74:77 79:82 85 87:90 92:95 ...
    97:104 106:115 121 123 127 128 130:143 145 146];

%% pre-processing of histology
%% It's faster to pre-process than transform 

% list of histology files
fileh = dir([HISTLO2BFREGRIGIDDIR filesep MOUSE '*.mha']);
N = length(fileh);

% allocate memory
imh = scimat_load([HISTLO2BFREGRIGIDDIR filesep fileh(1).name]);
imh(1:N) = imh;

for I = 1:N
    
    % load histology registered to blockface
    imh(I) = scimat_load([HISTLO2BFREGRIGIDDIR filesep fileh(I).name]);
    
end

% preprocess all images for registration
opts.Grayscale = true;
[~, imh, mask] = histology_preprocessing(imh(end), imh, opts);

for I = 1:N
    
    scimat_save([HISTLOPDIR filesep fileh(I).name], imh(I));
    scimat_save([HISTLOPDIR filesep 'mask-' fileh(I).name], mask(I));
    
end

%% intra-histology correction

% load all preprocessed histology and mask
for I = 1:N
    
    imh(I) = scimat_load([HISTLOPDIR filesep fileh(I).name]);
    mask(I) = scimat_load([HISTLOPDIR filesep 'mask-' fileh(I).name]);
    
end

% rigid transform diffusion registration of the histology volume
clear optReg optDiff
optReg.Angle = (-10:.25:10) / 180 * pi;
optReg.verbose = true;
optDiff.MaxIter = 50;
tic
[th2h, infoReg, infoRigDiff, imout] = transfdiffreg(...
    'regmatchedfilt', imh, optReg, optDiff);
toc

save([IMPROCDIR filesep T_INTRAHISTLO], ...
    'th2h', 'infoReg', 'infoRigDiff')


%% visually check histology to histology alignment

clear imh
opts.AutoDefaultPixelValue = true;

% compose transforms
load([IMPROCDIR filesep T_HISTLO2BF], 'th2bf')
load([IMPROCDIR filesep T_INTRAHISTLO], 'th2h')
ttot = elastix_cat(th2h, th2bf);

% init memory to store transformed histology
imh = zeros([ttot(1).Size([2 1]) 1 1 3 N], 'like', imh);

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
imagesc(permute(squeeze(imh(974, :, :, :, :, :)), [3 1 2]))
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

% load all of the histology images

clear imh
opts.AutoDefaultPixelValue = true;

% compose transforms
load([IMPROCDIR filesep T_HISTLO2BF], 'th2bf')
load([IMPROCDIR filesep T_INTRAHISTLO], 'th2h')
ttot = elastix_cat(th2h, th2bf);

% init memory to store transformed histology
[~, nameh] = fileparts(fileh(1).name);
imh = scimat_load([HISTLODIR filesep nameh '.tif']);
imh(1:N) = imh;

% apply previous transforms to original histology
for I = 1:N
    
    % load histology
    [~, nameh] = fileparts(fileh(I).name);
    imh(I) = scimat_load([HISTLODIR filesep nameh '.tif']);
    
    % apply transforms to histology
    imh(I) = transformix(ttot(I), imh(I), opts);
    
end

% preprocess all images for registration
opts.Grayscale = true;
mask = imh(1);
mask.data = uint8(mask.data(:, :, 1, 1, 1));
mask(1:N) = mask;
for I = 1:N
    
    [~, imh(I), mask(I)] ...
        = histology_preprocessing(imh(end), imh(I), opts);
    
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
save([IMPROCDIR filesep T_INTRAHISTLO], ...
    'th2h_bsp', 'infoBspReg', 'infoBspDiff', '-append')

%% visually check histology to histology alignment

clear imh
opts.AutoDefaultPixelValue = true;

% compose transforms
load([IMPROCDIR filesep T_HISTLO2BF], 'th2bf')
load([IMPROCDIR filesep T_INTRAHISTLO], 'th2h', 'th2h_bsp')
ttot = elastix_cat(th2h_bsp, th2h, th2bf);

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
imagesc(permute(squeeze(imh(974, :, :, :, :, :)), [3 1 2]))
set(gca, 'FontSize', 16)

% plot cross vertical section
subplot(2, 1, 2)
hold off
imagesc(permute(squeeze(imh(:, 967, :, :, :, :)), [3 1 2]))
set(gca, 'FontSize', 16)

saveas(gca, [FIGDIR filesep INTRAHISTLOBSP '-cross-sections.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Review 3D reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use auto estimation of background colour in transformix
optTransformix.AutoDefaultPixelValue = true;

% load blockface to histology
load([IMPROCDIR filesep T_HISTLO2BF], 'th2bf')

% load images and transforms
load([IMPROCDIR filesep T_INTRAHISTLO], 'th2h')
for I = 1:N
    im(:, :, :, I) = imread([HISTLOPDIR filesep fileh(I).name]);
end

% number of transform levels
L = elastix_length(th2h);

% compute measures of coefficient magnitude
clear tMax tMed
for L = 1:elastix_length(th2h)
    taux = elastix_colon(th2h, L);
    tMat = cat(1, taux(:).TransformParameters);
    tMed(:, L) = median(abs(tMat), 2);
    tMax(:, L) = max(abs(tMat), [], 2);
end

% loop levels of B-spline diffusion
for J = 1:L
    
    % keep the affine alignment (L) and the B-splines from J to L-1
    t = elastix_colon(th2h, J:L);
    
    % apply transform to the images
    for I = 1:N
        
        aux = imread([HISTLODIR filesep fileh(I).name]);
        aux = transformix(th2bf(I), aux, optTransformix);
        imout(:, :, :, I) = transformix(t(I), aux, optTransformix);
        
    end
    
    subplot(2, 2, 1)
    hold off
    imagesc(permute(squeeze(imout(619, :, :, :)), [3 1 2]))
    set(gca, 'FontSize', 16)
    title('Long axis free wall')
    
    % plot cross vertical section
    subplot(2, 2, 2)
    hold off
    imagesc(permute(squeeze(imout(840, :, :, :)), [3 1 2]))
    set(gca, 'FontSize', 16)
    title('Long axis apex')

    subplot(2, 2, 3)
    hold off
    imagesc(permute(squeeze(imout(1100, :, :, :)), [3 1 2]))
    set(gca, 'FontSize', 16)
    title('Long axis septum')

    subplot(2, 2, 4)
    hold off
    imagesc(permute(squeeze(imout(:, 800, :, :)), [3 1 2]))
    set(gca, 'FontSize', 16)
    title('Short axis mid-level')
    
    drawnow
    
    saveas(gca, [FIGDIR filesep HISTOLO '-reconstruct-diff-' ...
        num2str(optDiff.MaxIter) '-lvl-' num2str(L-J+1) '.png'])

end

% view imfusion of adjacent slices
for I = 1:N-1
    
    imagesc(imfuse(imout(:, :, :, I), imout(:, :, :, I+1)))
    pause
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply transforms to a small region of the high resolution histology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% propagate the reconstruction computed for the lo-res to the hi-res

% list of blockface and histology files
[filebf, fileh] = match_file_lists(...
    dir([BFPDATADIR filesep '*.mha']), ...
    dir([HISTLO2BFREGRIGIDDIR filesep '*.mha']), ...
    9:11, ...
    5:7 ...
    );
N = length(filebf);

clear imh
opts.AutoDefaultPixelValue = true;

% compose transforms
load([IMPROCDIR filesep T_HISTLO2BF], 'th2bf')
load([IMPROCDIR filesep T_INTRAHISTLO], 'th2h', 'th2h_bsp')
ttotlo = elastix_cat(th2h_bsp, th2h, th2bf);

for I = 1:N
    
    % convert transform to high resolution (undo the 5%=1/20 downsample)
%     ttot(I).Size = ttot(I).Size * 20;
    ttotlo(I).Size = [3579 4500];
    ttotlo(I).Spacing = ttotlo(I).Spacing / 20;
    
    % bottom corner of area we want to plot
    ttotlo(I).Origin = [0.007945 0.0043];

    % don't use compression, otherwise the java virtual machine runs out of
    % memory when we try to read the image
    ttotlo(I).CompressResultImage = 'false';
    
end

% save transforms that we apply to small portion of the histology
% save([IMPROCDIR filesep T_TOTHIST_LVFREE], 'ttotlo')

% load total transforms
load([HISTREFDIR_LVFREE filesep T_TOTHIST_LVFREE], 'ttotlo')

tic
for I = 1:N
    
    disp(['I = ' num2str(I)])
    
    % load histology
    [~, nameh] = fileparts(fileh(I).name);
    
    % apply transforms to histology using files (too large to store in
    % memory)
    opts.outfile = [HISTRECDIR filesep fileh(I).name];
    opts.AutoDefaultPixelValue = true;
    transformix(ttotlo(I), [HISTDIR filesep nameh '.tif'], opts);
%     imh(:, :, :, :, :, I) = aux.data;
    toc
end

% init memory to store transformed histology
imh = zeros([ttotlo(1).Size([2 1]) N 3], 'uint8');

for I = 1:N
    
    % load histology
    [~, nameh] = fileparts(fileh(I).name);
    aux = scimat_load([HISTRECDIR filesep fileh(I).name]);
    imh(:, :, I, :, :) = aux.data;
    
end

imagesc(squeeze(imh(:, :, 148, :)))

% plot cross horizontal section
hold off 
imagesc(permute(squeeze(imh(2500, :, :, :)), [2 1 3]))
set(gca, 'FontSize', 16)

% plot cross vertical section
imagesc(permute(squeeze(imh(:, 2000, :, :)), [2 1 3]))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Refine reconstruction in high resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% pre-process the hi-res histology so that we can refine the reconstruction

% matrix with a list of all the histology files
imfiles = [repmat([HISTRECDIR filesep], N, 1) cat(1, fileh.name)];

% convert each row into a cell
imfiles = mat2cell(imfiles, ones(1, N));

% preprocess all histology slices
clear opts
opts.Grayscale = true;
opts.outdir = HISTPDIR;
[~, imfiles, mask] = histology_preprocessing(imfiles{end}, imfiles, opts);

%% intra-histology correction

% in case we start from here, get list of preprocessed histology
imfiles = dir([HISTPDIR filesep 'Q*.mha']);
maskfiles = dir([HISTPDIR filesep 'mask*.mha']);
N = length(imfiles);
imfiles = [repmat([HISTPDIR filesep], N, 1) cat(1, imfiles.name)];
imfiles = mat2cell(imfiles, ones(1, N));
maskfiles = [repmat([HISTPDIR filesep], N, 1) cat(1, maskfiles.name)];
maskfiles = mat2cell(maskfiles, ones(1, N));

% remove the first 25 slices, because in the cropped histology they are
% empty or almost empty
maskfiles = maskfiles(25:end);
imfiles = imfiles(25:end);

% B-spline transform diffusion registration of the histology volume
clear optRegBsp optDiffBsp
optRegBsp.verbose = true;
optRegBsp.SpatialSamplesRatio = 0.1;
optRegBsp.MaxIter = 4;
%optReg.MaxVal = 1;
optRegBsp.mask = maskfiles;
optRegBsp.t0 = [];
optRegBsp.RegParam = [DOCDIR '/regParam-bspline-gray-intrahisto.txt'];
optDiffBsp.MaxIter = 100;
tic
[th2h_bsp, infoBspReg, infoBspDiff] = transfdiffreg(...
    'BSplineTransform', imfiles, optRegBsp, optDiffBsp);
toc

% add zero transforms for the slices we ignored
aux(1:N) = th2h_bsp(1);
for I = 1:24
    aux(I).TransformParameters(:) = 0;
end
aux(25:end) = th2h_bsp;
th2h_bsp = aux;

% save result
save([HISTREFDIR_LVFREE filesep T_TOTHIST_LVFREE], ...
    'th2h_bsp', 'infoBspReg', 'infoBspDiff', ...
    'optRegBsp', 'optDiffBsp', '-append')

%% apply transform to images, and visually check reconstruction refinement

clear imh opts
opts.AutoDefaultPixelValue = true;

% compose transforms
load([HISTREFDIR_LVFREE filesep T_TOTHIST_LVFREE], 'ttotlo', 'th2h_bsp')
% remove two registration iterations
th2h_bsp = elastix_colon(th2h_bsp, 3:4);
ttot = elastix_cat(th2h_bsp, ttotlo);

for I = 1:N
    
    tic
    disp(['I = ' num2str(I) '/' num2str(N)])
    
    % name of output file
    [~, nameh] = fileparts(fileh(I).name);
    opts.outfile = [HISTREFDIR_LVFREE filesep nameh '.mha'];
    
    % apply transforms to histology
    transformix(ttot(I), [HISTDIR filesep nameh '.png'], opts);
    
    disp(['(' num2str(toc) ' sec)'])
    
end

% init memory to store transformed histology
imh = zeros([ttot(1).Size([2 1]) N 1 3], 'uint8');

for I = 1:N
    
    % load transformed slice
    scimat = scimat_load([HISTREFDIR_LVFREE filesep fileh(I).name]);
    imh(:, :, I, :, :) = scimat.data;
    
end

% plot cross horizontal section
subplot(2, 1, 1)
hold off 
ra = imref2d([size(imh, 3) size(imh, 2)], scimat.axis(2).spacing, 10e-6);
imshow(permute(squeeze(imh(2000, :, :, :)), [2 1 3]), ra)
set(gca, 'FontSize', 16)

% plot cross vertical section
subplot(2, 1, 2)
hold off
ra = imref2d([size(imh, 3) size(imh, 1)], scimat.axis(1).spacing, 10e-6);
figure
imshow(permute(squeeze(imh(:, 3000, :, :)), [2 1 3]), ra)
set(gca, 'FontSize', 16)

% saveas(gca, [FIGDIR filesep INTRAHISTLOBSP '-cross-sections.png'])

%% apply only the lo res reconstruction to high resolution, without 
%% refinement, and crop to keep only the small block of the LV free wall


load([HISTREFDIR_LVFREE filesep T_TOTHIST_LVFREE], 'ttotlo')
ttot = ttotlo;

for I = 1:N
    
    tic
    disp(['I = ' num2str(I) '/' num2str(N)])
    
    % name of output file
    [~, nameh] = fileparts(fileh(I).name);
    opts.outfile = [HISTRECDIR_LVFREE filesep nameh '.mha'];
    
    % apply transforms to histology
    transformix(ttot(I), [HISTDIR filesep nameh '.png'], opts);
    
    disp(['(' num2str(toc) ' sec)'])
    
end

% init memory to store transformed histology
imh0 = zeros([ttot(1).Size([2 1]) N 1 3], 'uint8');

for I = 1:N
    
    % load transformed slice
    scimat = scimat_load([HISTRECDIR_LVFREE filesep fileh(I).name]);
    imh0(:, :, I, :, :) = scimat.data;
    
end

% plot cross horizontal section
subplot(2, 1, 1)
hold off 
ra0 = imref2d([size(imh0, 3) size(imh0, 2)], scimat.axis(2).spacing, 10e-6);
imshow(permute(squeeze(imh0(2000, :, :, :)), [2 1 3]), ra0)
set(gca, 'FontSize', 16)

% plot cross vertical section
subplot(2, 1, 2)
hold off
ra0 = imref2d([size(imh0, 3) size(imh0, 1)], scimat.axis(1).spacing, 10e-6);
imshow(permute(squeeze(imh0(:, 2000, :, :)), [2 1 3]), ra0)
set(gca, 'FontSize', 16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create video of the histology registered to the blockface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert -quality 100 -crop 891x733+885+475 *.png HistologyTrichromeLoRes.mpg

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
