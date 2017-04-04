% histology_refinement_without_diffusion.m
%
% Script to replace the transformation diffusion low-res B-spline
% refinement in "histology_blockface_mouse.m" by refinement using
% registration sweeps of the stack, without using diffusion. This
% corresponds to the "Reg. only" curve in Fig. 8 in the paper.
%
% version: 0.2.1
% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014-2017 University of Oxford
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


DEBUG = false;

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
HISTOLO_BF_P2_DIR = [HISTOLO_BF_DIR '_Bsp_Preprocessed'];

% transformations
T_HISTOLO_BF = [RESULT_DIR filesep 't_' HISTOLO '_to_Blockface.mat'];

% manual landmarks between blockface and histology
PTS_HIST2BF = [SRCDIR '/hand_tracing/pts_' HISTOLO '_to_Blockface.mat'];
PTS_HIST2HIST = [SRCDIR '/hand_tracing/pts_Intra_' HISTOLO '.mat'];

% landmark errors
ERR = [RESULT_DIR filesep 'err_' MOUSE '_' HISTOLO '.mat'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Intra-histology registration correction: B-spline without diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load histology after rigid refinement

% load histology slices and masks
clear imh mask
for I = 1:N
    
    disp(['I = ' num2str(I) ', Histology = ' fileh(I).name])
    imh(I) = scimat_load([HISTOLO_BF_P2_DIR filesep fileh(I).name]);
    mask(I) = scimat_load([HISTOLO_BF_P2_DIR filesep 'mask-' fileh(I).name]);
    
end

% B-spline transform diffusion registration of the histology volume
clear optReg_bsp optDiff_bsp
optReg_bsp_nodiff.verbose = true;
optReg_bsp_nodiff.SpatialSamplesRatio = 0.1;
optReg_bsp_nodiff.MaxIter = 50;
optReg_bsp_nodiff.mask = mask;
optReg_bsp_nodiff.t0 = [];
optReg_bsp_nodiff.RegParam = [SRCDIR '/regParam-bspline-gray-intrahisto.txt'];
optReg_bsp_nodiff.MaxVal = 0.0; % don't stop before reaching maximum number of iterations
optDiff_bsp_nodiff.Alpha = 0.45;
optDiff_bsp_nodiff.MaxIter = 1;
tic
[th2h_bsp_nodiff, infoReg_bsp_nodiff, infoDiff_bsp_nodiff, imout] = transfdiffreg(...
    'BSplineTransform', imh, optReg_bsp_nodiff, optDiff_bsp_nodiff);
toc

% save result
save(T_HISTOLO_BF, 'th2h_bsp_nodiff', 'infoReg_bsp_nodiff', 'infoDiff_bsp_nodiff', ...
    'optReg_bsp_nodiff', 'optDiff_bsp_nodiff', '-append')
