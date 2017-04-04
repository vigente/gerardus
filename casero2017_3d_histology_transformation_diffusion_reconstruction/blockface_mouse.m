% blockface_mouse.m
%
% version 0.3.14
%
% Prepare the blockface volume of mouse whole heart.

% % enable access to elastix shared object
% libpath = getenv('LD_LIBRARY_PATH');
% setenv('LD_LIBRARY_PATH', [ '/usr/local/bin/elastix/lib:' libpath]);

% select mouse to process
MOUSE = 'Q53';

% laptop at home
TOPDIR = '/home/rcasero/Software/private-gerardus-papers/casero2015_3d_histology_reconstruction_diffusion_velocity_field version';
DATADIR = ['/media/rcasero/mvldata/data/Mouse/' MOUSE];

% workstation at the office
TOPDIR = '/home/orie1416/Software/private-gerardus-papers/casero2015_3d_histology_reconstruction_diffusion_velocity_field version';
DATADIR = ['/data2/Mouse/' MOUSE];

% common directories for documentation
DOCDIR = [TOPDIR '/doc'];
SRCDIR = [TOPDIR '/src'];
FIGDIR = [TOPDIR '/doc/figures'];

% original blockface pictures
BF_DIR = [DATADIR '/Blockface'];

% directory to save processed images to
IMPROC_DIR = [DATADIR '/Image_Processing'];

% output blockface directories (in processing order)
BFS_DIR = [IMPROC_DIR '/Blockface_Stabilised'];
BFH_DIR = [IMPROC_DIR '/Blockface_Horizontal_Scratches'];
BFHC_DIR = [IMPROC_DIR '/Blockface_Horizontal_Scratches_Corrected'];
BFC_DIR = [IMPROC_DIR '/Blockface_Perspective_Corrected']; % perspective corrected and horizontal scratches
BFP_DIR = [IMPROC_DIR '/Blockface_Preprocessed'];

% binary masks for blockface stacks
BF_MASKS = [IMPROC_DIR filesep 'Blockface_Masks.mat'];

% transforms for a sequential pre-alignment of the blockface. Not all the
% transforms will be used, only those where there's an obvious jump. This
% is to avoid introducing registration noise between pairs of slices that
% are already well aligned
%
% Note: Slices 1-300 of the 55º block show horizontal heart while slice
% 301-640 show the apex pointing down. (The wax block was rotated 90º
% halfway while slicing).
%
% indices of frames to correct:
%   'idx90prop': (slices 1-300) jumps that propagate to the rest of the stack
%   'idx90prop2': (slices 301-640) ditto
%   'idx90noprop': (slices 1-300) jumps that only affect that slice
%   'idx90noprop2': (slices 301-640) ditto
%   'idx55prop'    |
%   'idx55prop2'   |  ditto for 55º
%   'idx55noprop'  |
%   'idx55noprop2' |
% indices of frames corrected manually:
%   'idx90prop2_manual': (slices 1-300)
%   'idx55prop_manual': (slices 1-640)
%   'idx55noprop_manual': (slices 1-640)
% transformation to correct 90º and 55º stack (help
% blockface_intraframe_reg for meaning of variables)
%   't90', 'tParam90', 'iterInfo90', 'regParam90'
%   't55', 'tParam55', 'iterInfo55', 'regParam55'
% projective transform to map the 55º stack onto the 90º stack coordinate
% frame
%   'tformPerspective'
% correct scratch angle
%   'scratchX', 'scratchY': scratchX(I), scratchY(I) are two points clicked
%                           on the same scratch, in a few frames of the
%                           (slices 1-300) stack. A line through this two
%                           points approximately follows the scratch
%   'alphamed': median angle of all the scratches
%   'tformScratch': transformation to make the scratches in the (slices
%                   1-300) stack horizontal
%   'scratchX2', 'scratchY2', 'alphamed2', 'tformScratch2': ditto for (301-640)
BF_INTRA_T = [IMPROC_DIR '/Blockface_Intra_Prealignment_Transform.mat'];

% histology directories for data. We need these just to read the pixel size
% from histology images at the end
% HISTO = 'HistologyTrichrome';
HISTO = 'HistologySiriusRed';
HISTOLO = [HISTO 'LoRes'];
HISTO_DIR = [DATADIR filesep HISTO];
HISTOLO_DIR = [DATADIR filesep HISTOLO];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Eyeballing detection of camera shift frames in blockface 55º and 90º
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open some image viewer application that allows you to scroll the
% blockface block, and make a note of the slices where there's a visible
% jump. For example, if the jump is between slice 53 and 54, then write
% down 54.
%
% This list of slices is in the "switch" block below:

switch MOUSE

    case 'Q53'
        
        % 90º stack
        % shifts that propagate to subsequent frames
        idx90prop = [202 221 222 231 232 241 242];
        idx90prop2 = [303:305, 323:327, 332, 334, 343:345, 347:349, ...
            392, 464:466, 469, 506:509, 586];
        % shifts that happen in an isolated frame, and don't need to
        % propagate
        idx90noprop = [266];
        idx90noprop2 = [301, 311, 316, 318, 320, 337, 340, 378, 383, ...
            456, 511, 596];
        
        % 55º stack
        % shifts that propagate to subsequent frames
        idx55prop = [66 68 94 176 238 285];
        idx55prop2 = [303 304 305 311 313 316 317 318 319 323 324 325 ...
            326 343 344 345 347 348 349 392 395 464 465 466 469 506 507 ...
            508 509 512 514 586 597];
        
        % shifts that happen in an isolated frame, and don't need to
        % propagate
        idx55noprop = [221 231 241 266];
        idx55noprop2 = [337 378 383];
    
    case 'Q62'
        
%         % only done up to frame 300
%         
%         % shifts that propagate to subsequent frames
%         idx55prop = [3 14 15 18 20 22 33 45 58 59 75 76 82 83 85 86 ...
%             89 90 96 97 106 107 108 109 137 141 144 145 151 152 157 ...
%             158 189 190 195 197 198 205 206 219 223 224];
%         idx90prop = [14 15 18 20 22 23 25 26 33 43 44 45 46 47 48 58 59 ...
%             64 66 68 69 70 71 75 76 82 83 84 85 86 89 90 96 97 102 103 ...
%             104 106 107 108 109 131 132 137 141 144 145 148 149 150 151 ...
%             152 153 154 155 158 168 169 170 171 172 189 190 192:195 197 ...
%             198 200:202 204:208 210:212 223 224 234 235];
%         
%         % shifts that happen in an isolated frame, and don't need to
%         % propagate
%         idx55noprop = [];
%         idx90noprop = [];
    
end

% save indices of frames to correct
save(BF_INTRA_T, 'idx90prop', 'idx90prop2', ...
    'idx90noprop', 'idx90noprop2')
save(BF_INTRA_T, 'idx55prop', 'idx55noprop', ...
    'idx55prop2', 'idx55noprop2', '-append')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stabilise blockface (remove jumps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% automatic alignment of 90º blockface

% we start aligning the 90º blockface, because the 55º blockface rotates
% 90º between slices 300 and 301, and that's going to need some special
% steps

% find intraframe transformations for the places where there are jumps
% (right now we don't worry whether the jumps are isolated or they transmit
% to the rest of the stack)
% 23434.722410 seconds (16 Jul 2016) for full stack
% 2883.161011 seconds (20 Jul 2016) selected slices
tic
[t90, tParam90, iterInfo90, regParam90] ...
    = blockface_intraframe_reg(BF_DIR, dir([BF_DIR '/*_90_*.bmp']), ...
    sort([idx90prop idx90prop2 idx90noprop idx90noprop2]));
toc 

% save transforms
save(BF_INTRA_T, ...
    't90', 'tParam90', 'iterInfo90', 'regParam90', '-append');

%% manual correction of selected slices in 90º blockface

% big changes in the microtome setting has created big jumps in image
% between slices:
%  (596 -> 597 is just a fail in segmentation)
%
% we correct those by selecting corresponding landmarks by hand between the
% bad slice and the previous one, and computing a similarity transformion
% to replace the automatic one
idx90prop_manual = [];
idx90prop2_manual = [301 303 469 509 596 597];

for I = unique([idx90prop_manual idx90prop2_manual])
    
    % open interface so that we can click corresponding points on slice(I)
    % and slice(I-1)
    [movingPoints, fixedPoints] = cpselect(...
        imread([BF_DIR '/Q53_90_0' num2str(I) '.bmp']), ... % moving
        imread([BF_DIR '/Q53_90_0' num2str(I-1) '.bmp']), ... % fixed
        'Wait', true);
    
    % compute similarity transformation from the hand segmented points. We
    % want similarity to account for rotations, zoom-ins and zoom-outs, and
    % translations
    tform = fitgeotrans(movingPoints, fixedPoints, 'NonreflectiveSimilarity');
    
    % elastix struct corresponding to tform transformation
    tElx = elastix_affine_matrix2struct(inv(tform.T), tParam90(I));
    
    if (DEBUG)
        
        % apply the computed transformation to the moving image
        moving = scimat_load([BF_DIR '/Q53_90_0' num2str(I) '.bmp']);
        moving_reg = transformix(tElx, moving);
        
        % load the fixed image and plot the overlap
        fixed = scimat_load([BF_DIR '/Q53_90_0' num2str(I-1) '.bmp']);
        close
        subplot(2, 1, 1)
        imshowpair(squeeze(fixed.data), squeeze(moving.data))
        subplot(2, 1, 2)
        imshowpair(squeeze(fixed.data), squeeze(moving_reg.data))
        
    end
    
    % if we accept the transformation, replace the automatically computed
    % by registration
    tParam90(I) = tElx;

end

% save updated transforms
save(BF_INTRA_T, 'tParam90', '-append');

% remove slices that we have registered by hand from the list of
% automatically registered ones. Also, because we have added slices to be
% registered, a couple have gone from being 
idx90prop2 = setdiff(idx90prop2, [303 469 509]);
idx90noprop2 = setdiff(idx90noprop2, [301 596]);

save(BF_INTRA_T, ...
    'idx90prop2', 'idx90prop_manual', 'idx90prop2_manual', ...
    'idx90noprop', 'idx90noprop2', '-append')

% apply transforms to create intermediate volume
blockface_correct_frame_shifts(BF_DIR, ...
    dir([BF_DIR '/*_90_*.bmp']), tParam90, ...
    unique([idx90noprop idx90noprop2]), ...
    BFS_DIR);


%% automatic alignment of 55º blockface

tic
[t55, tParam55, iterInfo55, regParam55] ...
    = blockface_intraframe_reg(BF_DIR, dir([BF_DIR '/*_55_*.bmp']), ...
    sort([idx55prop idx55prop2 idx55noprop idx55noprop2]));
toc % 1994.400186 seconds.

% save transforms
save(BF_INTRA_T, 't55', 'tParam55', 'iterInfo55', '-append');

%% manual correction of selected slices in 55º blockface (similarity transformation)

% apply transforms to create intermediate volume
blockface_correct_frame_shifts(BF_DIR, ...
    dir([BF_DIR '/*_55_*.bmp']), tParam55, ...
    unique([idx55noprop idx55noprop2]), ...
    BFS_DIR);

% we correct those by hand
idx55prop_manual = [312 320 321 325 327 332 340 341 343 344 392 470 471 514 586];
idx55noprop_manual = [383 456 511 596];

for I = unique([idx55prop_manual idx55noprop_manual])
    
    % load images
    moving = scimat_load([BF_DIR '/Q53_55_0' num2str(I) '.bmp']);
    fixed = scimat_load([BF_DIR '/Q53_55_0' num2str(I-1) '.bmp']);
    
    % open interface so that we can click corresponding points on slice(I)
    % and slice(I-1)
    [movingPoints, fixedPoints] = cpselect(...
        squeeze(moving.data), ... % moving
        squeeze(fixed.data), ... % fixed
        'Wait', true);
    
    % compute similarity transformation from the hand segmented points. We
    % want similarity to account for rotations, zoom-ins and zoom-outs, and
    % translations
    tform = fitgeotrans(movingPoints, fixedPoints, 'nonreflectivesimilarity');
    
    % elastix struct corresponding to tform transformation
    tElx = elastix_affine_matrix2struct(inv(tform.T), tParam55(I));
    
    if (DEBUG)
        
        % apply the computed transformation to the moving image
        moving_reg = transformix(tElx, moving);
        
        % load the fixed image and plot the overlap
        close
        subplot(2, 1, 1)
        imshowpair(squeeze(fixed.data), squeeze(moving.data))
        subplot(2, 1, 2)
        imshowpair(squeeze(fixed.data), squeeze(moving_reg.data))
        
    end
    
    % if we accept the transformation, replace the automatically computed
    % by registration
    tParam55(I) = tElx;
    
    pause

end

% because of hand corrections, some transforms that were non-propagating
% need to be relabelled as propagating

% save updated transforms
save(BF_INTRA_T, 'tParam55', 'idx55prop_manual', ...
    'idx55noprop_manual', '-append');

% apply transforms to create intermediate volume
blockface_correct_frame_shifts(BF_DIR, ...
    dir([BF_DIR '/*_55_*.bmp']), tParam55, ...
    unique([idx55noprop idx55noprop2 idx55noprop_manual]), ...
    BFS_DIR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find angle of horizontal wax block scratches by eyeballing scratches 
%% (slices 1-300), with the heart horizontally placed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% list of blockface images
files = dir([BFS_DIR filesep '*_55_*.bmp']);

% click a pair of points on the same scratch, and repeat for several
% scratches and several slices

x = cell(1, length(files));
y = cell(1, length(files));
for I = 1:60:300
    
    % load blockface image
    im = imread([BFS_DIR filesep files(I).name]);
    
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

% center of the heart in X, Y coordinates
% c = [1344 831];
c = [1339 960];

% create tform to correct the scratch angulation, so that scratches become
% horizontal
R = [cos(-alphamed) sin(-alphamed); ...
    -sin(-alphamed) cos(-alphamed)];
tformScratch = affine2d([...
    R [0; 0];
    c*(eye(2)-R) 1]);
    
% save scratch landmarks and rotation transform
scratchX = x;
scratchY = y;
save(BF_INTRA_T, 'scratchX', 'scratchY', 'alphamed', 'tformScratch', '-append')

N = length(dir([BFS_DIR '/*_55_*.bmp']));

% apply the transformation
for I = 1:300

    % load images
    moving = scimat_load([BFS_DIR '/Q53_55_0' num2str(I, '%03.0f') '.bmp']);
    
    % create reference frame
    rmoving = imref2d(size(moving.data));
    
    % apply the transformation
    moving_reg = moving;
    [moving_reg.data, rfixed] = imwarp(...
        squeeze(moving.data), rmoving, ...
        tformScratch, ...
        'OutputView', imref2d(size(fixed.data)));

    % save image
    imwrite(moving_reg.data, [BFH_DIR '/Q53_55_0' num2str(I, '%03.0f') '.bmp']);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find angle of vertical wax block scratches by eyeballing scratches 
%% (slices 301-640), with the heart vertically placed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% click a pair of points on the same scratch, and repeat for several
% scratches, and several slices
x = cell(1, N);
y = cell(1, N);
for I = 301:60:N
    
    % load blockface image
    im = imread([BFS_DIR filesep files(I).name]);
    
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
scratchX2 = x;
scratchY2 = y;
alphamed2 = alphamed;
tformScratch2 = tformScratch;
save(BF_INTRA_T, 'scratchX2', 'scratchY2', 'alphamed2', 'tformScratch2', '-append')

% apply the transformation
for I = 301:640

    % load images
    moving = scimat_load([BFS_DIR '/Q53_55_0' num2str(I, '%03.0f') '.bmp']);
    
    % create reference frame
    rmoving = imref2d(size(moving.data));
    
    % apply the transformation
    moving_reg = moving;
    [moving_reg.data, rfixed] = imwarp(...
        squeeze(moving.data), rmoving, ...
        tformScratch, ...
        'OutputView', imref2d(size(fixed.data)));

    % save image
    imwrite(moving_reg.data, [BFH_DIR '/Q53_55_0' num2str(I, '%03.0f') '.bmp']);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create masks to correct scratches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create blockface masks for the 55º stack

file = dir([BFH_DIR filesep '*_55_*.bmp']);

% create masks for the blockface volume (slices 1-300)
[ellipmask, polymask] = blockface_create_masks(BFH_DIR, file(1:300));

% create masks for the blockface volume (slices 301-640)
[ellipmask2, polymask2] = blockface_create_masks(BFH_DIR, file(301:640));

% save the masks
save(BF_MASKS, 'ellipmask', 'polymask', 'ellipmask2', 'polymask2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correct scratches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply the transformation to (sections 1-300)
for I = 1:300

    % load image
    im = imread([BFH_DIR filesep file(I).name]);
    
    % typical intensity value of the wax
    wax = double(im);
    wax(~polymask | ellipmask) = nan;
    imed = median(wax(~isnan(wax)));
    
    % correct scratches
    for J = 1:size(im, 1)
        
        % typical intensity value of the current row
        jmed = median(wax(J, ~isnan(wax(J, :))));
        
        % scale intensity value so that all rows have the same typical
        % intensity
        if (~isnan(jmed))
            im(J, :) = double(im(J, :)) - jmed + imed;
        end
        
    end

    % save image
    imwrite(im, [BFHC_DIR filesep file(I).name]);
    
end

% apply the transformation to (sections 301-640)
for I = 301:640

    % load image
    im = imread([BFH_DIR filesep file(I).name]);
    
    % typical intensity value of the wax
    wax = double(im);
    wax(~polymask2 | ellipmask2) = nan;
    imed = median(wax(~isnan(wax)));
    
    % correct scratches
    for J = 1:size(im, 1)
        
        % typical intensity value of the current row
        jmed = median(wax(J, ~isnan(wax(J, :))));
        
        % scale intensity value so that all rows have the same typical
        % intensity
        if (~isnan(jmed))
            im(J, :) = double(im(J, :)) - jmed + imed;
        end
        
    end

    % save image
    imwrite(im, [BFHC_DIR filesep file(I).name]);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute perspective correction of 55º block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% list of files
file55 = dir([BFHC_DIR filesep '*_55_*.bmp']);
file90 = dir([BFS_DIR filesep '*_90_*.bmp']);


%% (slices 1-300)

% slice to compute the projective transformation from
I = 150;

% load corresponding 55º and 90º slices
fixed = scimat_load([BFS_DIR filesep file90(I).name]);
moving = scimat_load([BFHC_DIR filesep file55(I).name]);

% create reference frame
rmoving = imref2d(size(moving.data));

% increase contrast
imax = double(max(moving.data(:)));
imin = double(min(moving.data(:)));
moving.data = uint8((double(moving.data) - imin) / (imax - imin) * 255);

imax = double(max(fixed.data(:)));
imin = double(min(fixed.data(:)));
fixed.data = uint8((double(fixed.data) - imin) / (imax - imin) * 255);

% open interface so that we can click corresponding points on slice(I)
% and slice(I-1)
[movingPoints, fixedPoints] = cpselect(...
    squeeze(moving.data), ... % moving
    squeeze(fixed.data), ... % fixed
    'Wait', true);

% compute a projective transformation to correct the rotation (with
% perspective)
tformPerspective = fitgeotrans(movingPoints, fixedPoints, 'projective');

% apply the transformation
moving_reg = moving;
[moving_reg.data, rfixed] = imwarp(moving.data, rmoving, tformPerspective, ...
    'OutputView', imref2d(size(fixed.data)));

% visual check that the perspective correction is fine
if (DEBUG)
    
    % load the fixed image and plot the overlap
    close
    subplot(2, 1, 1)
    imshowpair(squeeze(fixed.data), squeeze(moving.data))
    subplot(2, 1, 2)
    imshowpair(squeeze(fixed.data), squeeze(moving_reg.data))
    
end

% save projective transforms
movingPointsPerspective = movingPoints;
fixedPointsPerspective = fixedPoints;
save(BF_INTRA_T, 'tformPerspective', 'movingPointsPerspective', ...
    'fixedPointsPerspective', '-append');

% apply the transformation
for I = 1:300

    % load images
    moving = scimat_load([BFHC_DIR filesep file55(I).name]);
    
    % create reference frame
    rmoving = imref2d(size(moving.data));
    
    % apply the transformation
    moving_reg = moving;
    [moving_reg.data, rfixed] = imwarp(...
        squeeze(moving.data), rmoving, ...
        tformPerspective, ...
        'OutputView', imref2d(size(fixed.data)));

    % save image
    imwrite(moving_reg.data, [BFC_DIR filesep file55(I).name]);
    
end


%% (slices 301-640)

% we are going to match slice 301 at 55º to slice 300 with perspective
% corrected. This way, we also ensure continuity

% load corresponding 55º and 90º slices
fixed = scimat_load([BFC_DIR filesep file55(300).name]);
moving = scimat_load([BFHC_DIR filesep file55(301).name]);

% create reference frame
rmoving = imref2d(size(moving.data));

% increase contrast
imax = double(max(moving.data(:)));
imin = double(min(moving.data(:)));
moving.data = uint8((double(moving.data) - imin) / (imax - imin) * 255);

imax = double(max(fixed.data(:)));
imin = double(min(fixed.data(:)));
fixed.data = uint8((double(fixed.data) - imin) / (imax - imin) * 255);

% open interface so that we can click corresponding points on slice(I)
% and slice(I-1)
[movingPoints, fixedPoints] = cpselect(...
    squeeze(moving.data), ... % moving
    squeeze(fixed.data), ... % fixed
    'Wait', true);

% compute a projective transformation to correct the rotation (with
% perspective)
tformPerspective = fitgeotrans(movingPoints, fixedPoints, 'projective');

% apply the transformation
moving_reg = moving;
[moving_reg.data, rfixed] = imwarp(moving.data, rmoving, tformPerspective, ...
    'OutputView', imref2d(size(fixed.data)));

% visual check that the perspective correction is fine
if (DEBUG)
    
    % load the fixed image and plot the overlap
    close
    subplot(2, 1, 1)
    imshowpair(squeeze(fixed.data), squeeze(moving.data))
    subplot(2, 1, 2)
    imshowpair(squeeze(fixed.data), squeeze(moving_reg.data))
    
end

% save projective transforms
tformPerspective2 = tformPerspective;
movingPointsPerspective2 = movingPoints;
fixedPointsPerspective2 = fixedPoints;
save(BF_INTRA_T, 'tformPerspective2', 'movingPointsPerspective2', ...
    'fixedPointsPerspective2', '-append');

% apply the transformation
for I = 301:640

    % load images
    moving = scimat_load([BFHC_DIR filesep file55(I).name]);
    
    % create reference frame
    rmoving = imref2d(size(moving.data));
    
    % apply the transformation
    moving_reg = moving;
    [moving_reg.data, rfixed] = imwarp(...
        squeeze(moving.data), rmoving, ...
        tformPerspective, ...
        'OutputView', imref2d(size(fixed.data)));

    % save image
    imwrite(moving_reg.data, [BFC_DIR filesep file55(I).name]);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Downsample the histology to the size of the blockface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This needs to be done before cropping and illumination correction of the
% blockface, because cropping and illumination read the pixel information
% from the downsampled histology

% we want to downsample the histology to the size of the blockface so that
% both have the same pixel size. We open the same slice in both modalities,
% and measure the length between two easily identifiable landmarks
%
% The first half of the histology stack (1-300) was scanned at a different
% resolution than the second (301-640), so they need to be corrected
% separatedly

% list of blockface and histology files
[filebf, fileh] = match_file_lists(...
    dir([BFC_DIR filesep '*.bmp']), ...
    dir([HISTO_DIR filesep '*.tif']), ...
    9:11, ...
    5:7 ...
    );
N = length(filebf);

% list of files
filehcell = [];
filehlocell = [];
for I = 1:length(fileh)
    
    filehcell{I} = [HISTO_DIR filesep fileh(I).name];
    filehlocell{I} = [HISTOLO_DIR filesep fileh(I).name];
    
end

%% (slices 1-300)
imbf = imread([BFC_DIR filesep filebf(147).name]);
imh = imread([HISTO_DIR filesep fileh(147).name]);

subplot(2, 1, 1)
hold off
imagesc(imbf)
subplot(2, 1, 2)
hold off
imagesc(imh)

% longitudinal and width measures from the images
% blockface = [589.8551, 353.0127]
% histology = [1.616867641459869e+04, 9839]

% average K from the value computed from the histology/blockface pair
% above, and another one
K = mean([589.8551/1.616867641459869e+04, 353.0127/9839]);

%% (slices 301-640)
imbf = imread([BFC_DIR filesep filebf(152).name]);
imh = imread([HISTO_DIR filesep fileh(152).name]);

subplot(2, 1, 1)
hold off
imagesc(imbf)
subplot(2, 1, 2)
hold off
imagesc(imh)

% longitudinal and width measures from the images, slice 152
% blockface = [587.9124, 523.0468]
% histology = [9.1994e+03, 7.7570e+03]
% longitudinal and width measures from the images, slice 250
% blockface = [533.8352, 411.0195]
% histology = [8.5072e+03, 5.8856e+03]

K2 = mean([587.9124, 523.0468; 533.8352, 411.0195] ...
    ./ [9.1994e+03, 7.7570e+03; 8.5072e+03, 5.8856e+03]);

%% resize the histology files
im_resize(filehcell(1:148), filehlocell(1:148), K);
im_resize(filehcell(149:end), filehlocell(149:end), K2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cropping and illumination correction of blockface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% list of blockface and histology files
[filebf, fileh] = match_file_lists(...
    dir([BFC_DIR filesep '*.bmp']), ...
    dir([HISTOLO_DIR filesep '*.tif']), ...
    9:11, ...
    5:7 ...
    );
N = length(filebf);

%% Create blockface masks for the perspective corrected 55º stack

% create masks for the blockface volume
[ellipmask90, polymask90] = blockface_create_masks(BFC_DIR, filebf);

% save the masks
save(BF_MASKS, 'ellipmask90', 'polymask90', '-append')

%% cropping and scaling of blockface

% crop factors for blockface images
crop_row = 828-485:828+485;
crop_col = 1302-549:1302+549;

% load mask to equalise blockface illumination
load(BF_MASKS, 'ellipmask90', 'polymask90')

% load blockface images
info = imfinfo([BFC_DIR filesep filebf(1).name]);
imbf = zeros(info.Height, info.Width, N, ['uint' num2str(info.BitDepth)]);
for I = 1:N
    imbf(:, :, I) = imread([BFC_DIR filesep filebf(I).name]);
end

% crop blockface data
ellipmask90Crop = ellipmask90(crop_row, crop_col);
polymask90Crop = polymask90(crop_row, crop_col);
imbf = imbf(crop_row, crop_col, :);

% save masks
save(BF_MASKS, 'ellipmask90Crop', 'polymask90Crop', ...
    'crop_row', 'crop_col', '-append')

%% pre-processing of blockface images

% correction algorithm parameters
ratio = 1/16;
thr = 45;
radheart = 10;
radpoly = 50;

for I = 1:N
    disp(['I = ' num2str(I)])
    
    % load histology so that we have pixel size values
    scimath = scimat_load([HISTOLO_DIR filesep fileh(I).name]);

    % equalise blockface image
    scimath.data = blockface_equalise_illumination(...
        imbf(:, :, I), polymask90Crop, ellipmask90Crop, ...
        ratio, thr, radheart, radpoly);
    scimath.axis(1).size = size(scimath.data, 1);
    scimath.axis(2).size = size(scimath.data, 2);
    scimath.axis(1).min = -scimath.axis(1).spacing / 2;
    scimath.axis(2).min = -scimath.axis(2).spacing / 2;

    % save blockface
    [~, auxfilename] = fileparts(filebf(I).name);
    scimat_save([BFP_DIR filesep auxfilename '.mha'], scimath);
    scimat_save([BFP_DIR filesep auxfilename '.tif'], scimath);
end

