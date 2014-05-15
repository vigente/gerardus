function [block, mfg, mbg] = amrf_seg(im, blocklen, overlap)
% AMRF_SEG  Adaptive Markov Random Field segmentation.
%
%   Note: This function is in development, and can produce "blocky"
%   segmentations.
%
%   This function splits an image into blocks, estimates typical background
%   and foreground intensities, and uses them to seed a Markov Random Field
%   segmentation algorithm from the Insight Toolbox (itk::MRFImageFilter).
%
% BW = amrf_seg(IM, BLOCKLEN)
%
%   IM is a 3D array with grayscale intesity values, forming a 3D image
%   (e.g. Magnetic Resonance). The background is assumed to be light
%   (higher intensity values), and the foreground dark (lower intensity
%   values).
%
%   IM is split into blocks to estimate the typical foreground and
%   background intensity values.
%
%   BLOCKLEN is a 3-vector with the number of rows, columns and slices in
%   each block. If an integer number of blocks does not fit in the image,
%   the last blocks in each dimension are cropped.
%
%   BW is a binary segmentation of IM. Foreground voxels = 1, and
%   background voxels = 0.
%
% [..., MFG, MBG] = amrf_seg(..., OVERLAP)
%
%   To improve the smoothness of the segmentation result, we don't want
%   sudden changes in the MUFG, MUGB between neighbouring blocks.
%
%   OVERLAP is a 3-vector with the number of rows, columns and slices added
%   on each side of each block so that it overlaps with its neighbours. By
%   default, OVERLAP = [0 0 0] (no overlap is used). Note that large
%   overlaps will increase the memory used by this function substantially.
%
%   MFG, MBG are 3D arrays with the typical (modes) foreground and
%   background intensities for each of the estimation blocks.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.1.0
% $Rev$
% $Date$
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

% check arguments
narginchk(2, 3);
nargoutchk(0, 3);

% defaults
if (nargin < 3 || isempty(overlap))
    overlap = [0 0 0];
end

% image dimensions
NR = size(im, 1);
NC = size(im, 2);
NS = size(im, 3);

%% Split image into overlapping blocks and estimate foreground and
%% background

% "botton left" index of every block in the image, without overlap. If the
% last block starts at the right edge of the image, it has size zero. Thus,
% if we have any block starting there, we remove it with setdiff()
br0 = setdiff(1:blocklen(1):NR, NR);
bc0 = setdiff(1:blocklen(2):NC, NC);
bs0 = setdiff(1:blocklen(3):NS, NS);

% "top right" index of every block in the image, without overlap (taking
% care that the last boxes don't go outside the image)
tr0 = min(NR, br0 + blocklen(1) - 1);
tc0 = min(NC, bc0 + blocklen(2) - 1);
ts0 = min(NS, bs0 + blocklen(3) - 1);

% the last blocks have to cover to the end of the image
tr0(end) = NR;
tc0(end) = NC;
ts0(end) = NS;

% extend the blocks with the overlap margins
br = max(1, br0 - overlap(1));
bc = max(1, bc0 - overlap(2));
bs = max(1, bs0 - overlap(3));
tr = min(NR, tr0 + overlap(1));
tc = min(NC, tc0 + overlap(2));
ts = min(NS, ts0 + overlap(3));

% split the image into blocks. This is necessary if we want to process the
% blocks in parallel. We create a vector array of blocks instead of a 3D
% array of blocks because parfor cannot be nested
block = cell(1, length(br) * length(bc) * length(bs));
for ri = 1:length(br)
    for ci = 1:length(bc)
        for si = 1:length(bs)
            
            block{sub2ind([length(br) length(bc) length(bs)], ri, ci, si)} ...
                = im(br(ri):tr(ri), bc(ci):tc(ci), bs(si):ts(si));
            
        end
    end
end

% parallel processing of the blocks
mfg = zeros(length(br), length(bc), length(bs));
mbg = zeros(length(br), length(bc), length(bs));
parfor I = 1:numel(block)

    [mfg(I), mbg(I)] = im_modes(block{I});
    
end

% clear memory
clear block

% reshape outputs as 3D arrays
mfg = reshape(mfg, length(br), length(bc), length(bs));
mbg = reshape(mbg, length(br), length(bc), length(bs));

%% Interpolate typical values of foreground and background intensity that 
%% couldn't be estimated in the previous phase

% values that need to be interpolated
fgidx = isnan(mfg);
bgidx = isnan(mbg);

if (nnz(fgidx) || nnz(bgidx))

    % centers of each block without overlap
    [gr0, gc0, gs0] ...
        = ndgrid((br0 + tr0)*0.5, (bc0 + tc0)*0.5, (bs0 + ts0)*0.5);
    
    % scale the center coordinates to be within the [0, 1] cube, to avoid
    % ill-conditioned L matrices in the thin-plate spline. Note that
    % because of the scale-invariance of the TPS, this doesn't change the
    % result
    K = max([(br0(end) + tr0(end))*0.5, ...
        (bc0(end) + tc0(end))*0.5, ...
        (bs0(end) + ts0(end))*0.5]);
    gr0 = gr0 / K;
    gc0 = gc0 / K;
    gs0 = gs0 / K;
    
    % interpolate foreground missing values
    mfg(fgidx) = pts_tps_map([gr0(~fgidx), gc0(~fgidx), gs0(~fgidx)], ...
        mfg(~fgidx), [gr0(fgidx), gc0(fgidx), gs0(fgidx)]);
    
    % interpolate background missing values
    mbg(bgidx) = pts_tps_map([gr0(~bgidx), gc0(~bgidx), gs0(~bgidx)], ...
        mbg(~bgidx), [gr0(bgidx), gc0(bgidx), gs0(bgidx)]);
    
end

%% Split the image into non-overlapping blocks, and segment them using a
%% Markov Random Field method

% split the image into blocks, this time without overlap
block = cell(1, length(br) * length(bc) * length(bs));
for ri = 1:length(br0)
    for ci = 1:length(bc0)
        for si = 1:length(bs0)
            
            block{sub2ind([length(br0) length(bc0) length(bs0)], ri, ci, si)} ...
                = im(br0(ri):tr0(ri), bc0(ci):tc0(ci), bs0(si):ts0(si));
            
        end
    end
end

% Markov Random Field segmentation of the blocks
parfor I = 1:numel(block)

    block{I} = itk_imfilter('mrf', block{I}, [mfg(I) mbg(I)]);
    
end

% reshape the vector of blocks as 3D array
block = reshape(block, length(br), length(bc), length(bs));

% format output as segmentation, and invert so that tissue=1, background=0
block = ~cell2mat(block);
