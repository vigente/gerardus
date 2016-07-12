% validate_reconstruction_noblockface.m
%
% This is a validation script. We assume that
%
%   * a histology reconstruction has been computed with
%     histology_noblockface_mouse.m
%
%   * manual landmarks between histology/histology
%     slices have been traced using script manual_landmarks.m.
%
% Then, this script will apply the reconstruction transforms to the
% landmarks, and compute how far they are from each other at each
% refinement step.
%
% Note that we are validating the reconstruction without blockface, but we
% still need the blockface landmarks so that we can measure how much the
% reconstruction differs from the true anatomy.

% version: 0.1.4
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

PROJDIR = '/home/orie1416/Software/casero2015_3d_histology_diffusion_registration';
DOCDIR = [PROJDIR '/doc'];
SRCDIR = [PROJDIR '/src'];
FIGDIR = [DOCDIR '/figures'];

MOUSE = 'Q53';

HISTO = 'HistologySiriusRed';
HISTOLO = [HISTO 'LoRes'];
HISTDIR = ['/data2/Mouse/' MOUSE filesep HISTO];
HISTLODIR = ['/data2/Mouse/' MOUSE filesep HISTOLO];
IMPROCDIR = ['/data2/Mouse/' MOUSE '/ImageProcessing'];
BFCDATADIR = [IMPROCDIR '/BlockfaceCorrected']; % alignment only
BFPDATADIR = [IMPROCDIR '/BlockfacePreprocessed']; % cropped, intensity equalisation, inversion

% manually traced landmarks
PTS_HIST2BF = ['pts_' HISTOLO '_to_Blockface.mat'];
PTS_HIST2HIST = ['pts_Intra_' HISTOLO '.mat'];

%% transforms to validate

T_NOBF_INTRAHISTLO = ['t_NoBlockfaceIntra_' HISTOLO '.mat'];
T_NOBF_INTRAHISTLO_FILE = [IMPROCDIR filesep T_NOBF_INTRAHISTLO];

T_INTRAHIST = T_NOBF_INTRAHISTLO;
T_INTRAHIST_FILE = T_NOBF_INTRAHISTLO_FILE;

% output error files
ERR_BLOCKHIST = ['err_' MOUSE '_blockface_histology_' T_INTRAHIST];
ERR_INTRAHIST = ['err_' MOUSE '_intra_histology_' T_INTRAHIST];

%% load transforms to validate

% load histology -> blockface transformation
load(T_NOBF_INTRAHISTLO_FILE, 'th2h', 'th2h_bsp')

% compose transforms
ttot = elastix_cat(th2h_bsp, th2h);

% number of levels
L = elastix_length(ttot);

%% blockface and histology images and landmarks

% list of blockface and histology files
[filebf, fileh] = match_file_lists(...
    dir([BFCDATADIR filesep '*.png']), ...
    dir([HISTLODIR filesep '*.tif']), ...
    9:11, ...
    5:7 ...
    );
N = length(filebf);

% load landmarks
load([IMPROCDIR filesep PTS_HIST2BF], 'ptsbf', 'ptsh')
load([IMPROCDIR filesep PTS_HIST2HIST], 'ptsh1', 'ptsh2')

%% reference blockface slice for the next two error computation sections

% crop factors for blockface images
crop_row = 284:1254;
crop_col = 766:1864;

% size scaling = histology / blockface
% this is estimated by measuring the same antomical distance in a blockface
% image and the corresponding histology image
sc = 866.5806 / 594.4855;

% load one of the blockface images (all of them end up with the same size
% and pixel size)
imh = scimat_load([HISTLODIR filesep fileh(1).name], 'HeaderOnly', true);
imbf = scimat_im2scimat(...
    imread([BFCDATADIR filesep filebf(1).name]), ...
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

%% blockface-histology landmark error

% because we are not aligning histology to blockface, to measure the error
% we have to 
% 1) do the histology reconstruction,
% 2) apply transforms to histology landmarks,
% 3) rigid alignment of landmarks to blockface landmarks. This removes the
%    rigid  arbitrary difference between blockface and histology, and
%    leaves the true reconstruction error
%
% When we use the blockface, we can compute the landmark error directly.

% indices of pairs of slices with landmarks
idx = find(~cellfun(@isempty, ptsbf));

err = cell(1, length(fileh));
for I = idx
    
    disp(['I = ' num2str(I) '/' num2str(idx(end))])
    
    % load histology image metadata
    imh = scimat_load([HISTLODIR filesep fileh(I).name], 'HeaderOnly', true);
    
    % allocate image space to hold the landmarks
    imh.data = zeros([imh.axis(1).size, imh.axis(2).size], 'uint8');
    
    % row, col coordinates of the landmarks
    rc = round(scimat_world2index(ptsh{I}, imh));
    
    % transfer landmarks to image. Each landmarks becomes a 3x3 square in
    % the image. The reason why we cannot use a single pixel per landmarks
    % is because otherwise the transformation interpolation sometimes loses
    % them
    for K = 1:size(rc, 1)
        imh.data(rc(K, 1)-1:rc(K, 1)+1, rc(K, 2)-1:rc(K, 2)+1) = K;
    end
    
    % scale the output to the resolution of the histology, to avoid
    % losing too much precising in the landmarks
    xylim = scimat_index2world(size(imbf.data), imbf);
    rcsize = scimat_world2index(xylim, imh, false);
    
    % loop levels of refinement
    for J = L:-1:1
        
        disp(['J = ' num2str(J) '/' num2str(L)])
        
        % keep only the first levels of transform. Because of how elastix
        % defines transforms, the translation is the last level, which
        % appears as the first level in the struct
        taux = elastix_colon(ttot(I), J:L);
        
        % when we transform the landmarks, we want to keep only their
        % integer values, not get fractional values from linear/bspline
        % interpolation
        taux.ResampleInterpolator = 'FinalNearestNeighborInterpolator';
        taux.DefaultPixelValue = 0;
        taux.CompressResultImage = 'false';
        
        % scale the output to the resolution of the histology, to avoid
        % losing too much precising in the landmarks
        taux.Spacing = [imh.axis(2).spacing imh.axis(1).spacing];
        taux.Size = rcsize([2 1]);
        
        % apply transform to landmark image
        opt.verbose = false;
        imh2 = transformix(taux, imh, opt);
        
        % find the location of the transformed landmarks, and sort them in
        % increasing order
        loc = find(imh2.data);
        [r, c] = ind2sub(size(imh2.data), loc);
        xxyy = scimat_index2world([r, c], imh2);
        lab = imh2.data(loc);
        xy{I, J} = [];
        for K = 1:size(rc, 1)
            if (nnz(lab == K)==0)
                warning(['Landmark was lost in transformation, I = ' ...
                    num2str(I) ', J = ' num2str(J) ', K = ' num2str(K)])
            end
            % if a landmark if mapped to more than one pixel, average them
            xy{I, J}(K, :) = mean(xxyy(lab == K, :), 1);
        end
        
    end
    
end

% compute landmark error
for J = L:-1:1

    % rigid alignment of transformed histology points to the blockface. We
    % have to align the landmaks from all slices to get a global rigid
    % alignment
    [~, ~, trig] = procrustes(cat(1, ptsbf{:}), cat(1, xy{:, J}), ...
        'Scaling', false, 'Reflection', false);
   
    for I = idx
        
        % apply the correction to each slice to compute error between
        % blockface and histology landmarks
        aux = trig.b * xy{I, J} * trig.T + trig.c(1:size(xy{I, J}, 1), :);
        
        % compute error for the landmarks in slice I and refinement level J
        err{I}(:, J) = sqrt(sum((aux - ptsbf{I}).^2, 2));
    end
    
end

% save error results
save([SRCDIR filesep ERR_BLOCKHIST], 'err')

% concatenate error by landmarks so that we can see how it changes from
% level to level of refinement. Also, reorder columns front-to-back to that
% we start with the coarsest refinement level
errtot = cat(1, err{:});
errtot = errtot(:, end:-1:1);

subplot(2, 1, 1)
hold off
boxplot(errtot)
subplot(2, 1, 2)
hold off
plot(median(errtot, 1))

%% histology-histology landmark error

% for this one, we don't need to align rigidly to the blockface, as
% intra-histology distances don't depend on the blockface.

% indices of pairs of slices with landmarks
idx = find(~cellfun(@isempty, ptsh1));

err = cell(1, length(fileh));
for I = idx
    
    disp(['I = ' num2str(I) '/' num2str(idx(end))])
    
    % load histology image metadata
    imh1 = scimat_load([HISTLODIR filesep fileh(I).name], 'HeaderOnly', true);
    imh2 = scimat_load([HISTLODIR filesep fileh(I+1).name], 'HeaderOnly', true);
    
    % allocate image space to hold the landmarks
    imh1.data = zeros([imh1.axis(1).size, imh1.axis(2).size], 'uint8');
    imh2.data = zeros([imh2.axis(1).size, imh2.axis(2).size], 'uint8');
    
    % row, col coordinates of the landmarks
    rc1 = round(scimat_world2index(ptsh1{I}, imh1));
    rc2 = round(scimat_world2index(ptsh2{I}, imh2));
    
    % transfer landmarks to image. Each landmarks becomes a 3x3 square in
    % the image. The reason why we cannot use a single pixel per landmarks
    % is because otherwise the transformation interpolation sometimes loses
    % them
    for K = 1:size(rc1, 1)
        imh1.data(rc1(K, 1)-1:rc1(K, 1)+1, rc1(K, 2)-1:rc1(K, 2)+1) = K;
        imh2.data(rc2(K, 1)-1:rc2(K, 1)+1, rc2(K, 2)-1:rc2(K, 2)+1) = K;
    end
    
    % scale the output to the resolution of the histology, to avoid
    % losing too much precising in the landmarks
    xylim = scimat_index2world(size(imbf.data), imbf);
    rcsize = scimat_world2index(xylim, imh1, false);
    
    % initialize error
    err{I} = zeros(size(rc1, 1), L);
    
    for J = L:-1:1
        
        disp(['J = ' num2str(J) '/' num2str(L)])
        
        % keep only the first levels of transform. Because of how elastix
        % defines transforms, the translation is the last level, which
        % appears as the first level in the struct
        taux1 = elastix_colon(ttot(I), J:L);
        taux2 = elastix_colon(ttot(I+1), J:L);
        
        % when we transform the landmarks, we want to keep only their
        % integer values, not get fractional values from linear/bspline
        % interpolation
        taux1.ResampleInterpolator = 'FinalNearestNeighborInterpolator';
        taux2.ResampleInterpolator = 'FinalNearestNeighborInterpolator';
        taux1.DefaultPixelValue = 0;
        taux2.DefaultPixelValue = 0;
        taux1.CompressResultImage = 'false';
        taux2.CompressResultImage = 'false';
        
        % scale the output to the resolution of the histology, to avoid
        % losing too much precising in the landmarks
        taux1.Spacing = [imh1.axis(2).spacing imh1.axis(1).spacing];
        taux2.Spacing = [imh2.axis(2).spacing imh2.axis(1).spacing];
        taux1.Size = rcsize([2 1]);
        taux2.Size = rcsize([2 1]);
        
        % apply transform to landmark image
        opt.verbose = false;
        imh1tr = transformix(taux1, imh1, opt);
        imh2tr = transformix(taux2, imh2, opt);
        
        % find the location of the transformed landmarks, and sort them in
        % increasing order (slice I)
        loc = find(imh1tr.data);
        [r, c] = ind2sub(size(imh1tr.data), loc);
        xxyy = scimat_index2world([r, c], imh1tr);
        lab = imh1tr.data(loc);
        xy1{I,J} = [];
        for K = 1:size(rc1, 1)
            if (nnz(lab == K)==0)
                warning(['Landmark was lost in transformation, I = ' ...
                    num2str(I) ', J = ' num2str(J) ', K = ' num2str(K)])
            end
            % if a landmark if mapped to more than one pixel, average them
            xy1{I,J}(K, :) = mean(xxyy(lab == K, :), 1);
        end
        
        % find the location of the transformed landmarks, and sort them in
        % increasing order (slice I+1)
        loc = find(imh2tr.data);
        [r, c] = ind2sub(size(imh2tr.data), loc);
        xxyy = scimat_index2world([r, c], imh2tr);
        lab = imh2tr.data(loc);
        xy2{I,J} = [];
        for K = 1:size(rc2, 1)
            if (nnz(lab == K)==0)
                warning(['Landmark was lost in transformation, I = ' ...
                    num2str(I) ', J = ' num2str(J) ', K = ' num2str(K)])
            end
            % if a landmark if mapped to more than one pixel, average them
            xy2{I,J}(K, :) = mean(xxyy(lab == K, :), 1);
        end
        
    end
    
    % compute error between histology landmarks
    for J = L:-1:1
        
        err{I}(:, J) = sqrt(sum((xy1{I,J} - xy2{I,J}).^2, 2));
        
    end

%     for J = L:-1:1
%         hold off
%         plot(xy{I,J}(:, 1), xy{I,J}(:, 2), 'kx', 'LineWidth', 1)
%         hold on
%         plot(ptsbf{I}(:, 1), ptsbf{I}(:, 2), 'ko', 'LineWidth', 1)
%         pause
%     end
%     
%     % plot image with landmarks
%     subplot(2, 1, 1)
%     hold off
%     scimat_imagesc(foo);
%     hold on
%     plot(ptsh{I}(:, 1), ptsh{I}(:, 2), 'wo', 'LineWidth', 3)
%     plot(ptsh{I}(:, 1), ptsh{I}(:, 2), 'kx', 'LineWidth', 2)    
%     subplot(2, 1, 2)
%     hold off
%     scimat_imagesc(imh);

end

% save error results
save([SRCDIR filesep ERR_INTRAHIST], 'err')

% concatenate error by landmarks so that we can see how it changes from
% level to level of refinement. Also, reorder columns front-to-back to that
% we start with the coarsest refinement level
errtot = cat(1, err{:});
errtot = errtot(:, end:-1:1);

subplot(2, 1, 1)
hold off
boxplot(errtot)
subplot(2, 1, 2)
hold off
plot(median(errtot, 1))
