% histology_refinement_without_diffusion.m
%
%
% Script to replace the low-res B-spline refinement in
% "histology_blockface_mouse.m":
%
%   - B-spline diffusion refinement of low-res intra-histology.
%     (20 registration sweeps, no neighbor updates).
%
% version: 0.1.3
% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014-2015 University of Oxford
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

DOCDIR = '/home/orie1416/Software/private-gerardus-papers/casero2015_3d_histology_diffusion_registration/doc';
SRCDIR = [DOCDIR '../src'];
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

T_HISTLO2BF = ['t_' HISTOLO '_to_Blockface.mat'];
T_INTRAHISTLO = ['t_Intra_' HISTOLO '.mat'];
T_INTRAHISTLONODIFF = ['t_IntraNoDiff_' HISTOLO '.mat'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Intra-histology registration correction: B-spline without diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
load([IMPROCDIR filesep T_INTRAHISTLO], 'th2h')
ttot = elastix_cat(th2h, th2bf);

% init memory to store transformed histology
[~, nameh] = fileparts(fileh(1).name);
imh = scimat_load([HISTLODIR filesep nameh '.png']);
imh(1:N) = imh;

% apply previous transforms to original histology
for I = 1:N
    
    % load histology
    [~, nameh] = fileparts(fileh(I).name);
    imh(I) = scimat_load([HISTLODIR filesep nameh '.png']);
    
    % apply transforms to histology
    imh(I) = transformix(ttot(I), imh(I), opts);
    
end

% preprocess all images for registration
mask = imh(1);
mask.data = uint8(mask.data(:, :, 1, 1, 1));
mask(1:N) = mask;
for I = 1:N
    
    [~, imh(I), mask(I)] ...
        = histology_preprocessing(imh(end), imh(I));
    
end
for I = 1:N
    
    imh(I).data = rgb2gray(squeeze(imh(I).data));
    
end


%% intra-histology correction. We don't allow diffusion, so we run many 
%% more registration iterations

% B-spline transform diffusion registration of the histology volume
clear optReg optDiff
optReg.verbose = true;
optReg.SpatialSamplesRatio = 0.1;
optReg.MaxIter = 20;
%optReg.MaxVal = 1;
optReg.mask = mask;
optReg.t0 = [];
optReg.RegParam = [SRCDIR '/regParam-bspline-gray-intrahisto.txt'];
optDiff.MaxIter = 1;
tic
[th2h_bsp, infoBspReg, infoBspDiff, imout] = transfdiffreg(...
    'BSplineTransform', imh, optReg, optDiff);
toc


% save result
save([IMPROCDIR filesep T_INTRAHISTLONODIFF], ...
     'th2h', 'th2h_bsp', 'infoBspReg', 'infoBspDiff')
