function [imh, imbf] = histology_blockface_preprocessing(imh, imbf, polymask, ellipmask)
% histology_blockface_preprocessing  Prepare histology and blockface for
% registration.
%
% histology_blockface_preprocessing takes a histology and blockface images
% that are going to be registered onto each other, and: 1) inverts them, to
% have a dark background, 2) corrects the scratches of the blockface, 3)
% thresholds the background, to remove noise, 4) extends the histograms to
% cover the intensity range and 5) matches the histology histogram to the
% blockface histogram.
%
% This function should be run after homogenising the illumination in the
% blockface with blockface_equalise_illumination. It is also assumed that
% the blockface has been rotated to make the scratches horizontal.
%
% [IMH2, IMBF2] = histology_blockface_preprocessing(IMH, IMBF, POLYMASK, ELLIPMASK)
%
%   IMH, IMBF are colour or grayscale images with the histology and
%   blockface, respectively.
%
%   POLYMASK, ELLIPMASK are segmentation masks for the blockface image
%   IMBF. POLYMASK is an outline of the wax blockface. ELLIPMASK is an
%   ellipse that contains the heart. The masks can be created with
%   blockface_create_masks.
%
%   IMH2, IMBF2 are the corresponding images after preprocessing.
%
% See also: blockface_create_masks, blockface_equalise_illumination.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.2.0
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

DEBUG = 0;

% check arguments
narginchk(4, 4);
nargoutchk(0, 2);

if (any(size(imbf) ~= size(polymask)))
    error('POLYMASK must be the same size as IMBF')
end
if (any(size(imbf) ~= size(ellipmask)))
    error('ELLIPMASK must be the same size as IMBF')
end

%% preprocessing of histology

% make sure we have grayscale images
if (size(imh, 3) == 3)
    imh = rgb2gray(imh);
end
if (size(imbf, 3) == 3)
    imbf = rgb2gray(imbf);
end

% plot slices
if (DEBUG)
    subplot(2, 1, 1)
    hold off
    imagesc(imh)
    subplot(2, 1, 2)
    hold off
    imagesc(imbf)
end

%% Histology preprocessing

% invert image and make lowest intensity = 0
imh = max(imh(:)) - imh;

% remove background by zeroing anything <= mode (the background forms the
% largest peak). We are quite conservative here in terms of keeping all the
% tissue we can, to avoid losing detail, even if that means that we are
% going to have a bit of background noise
imh(imh <= mode(double(imh(:)))) = 0;

% extend histogram to cover whole dynamic range
minh = double(min(imh(imh > 0)));
maxh = double(max(imh(imh > 0)));
imh = uint8(255 * (double(imh) - minh) / (maxh - minh));

% plot preprocessed histology
if (DEBUG)
    subplot(2, 1, 1)
    imagesc(imh)
end

%% Blockface preprocessing

% invert image
imbf = max(imbf(:)) - imbf;

% correct scratches

% duplicate of the image with intensity values only on the block without
% the heart area, and the rest as 0
imbf2 = imbf;
imbf2(~xor(polymask, ellipmask)) = 0;

% median intensity value of all intensity values on the block
vtarget = double(median(imbf2(imbf2 > 0)));

% iterate rows
for I = 1:size(imbf, 1)
   
    % intensity values in this row that are not 0
    v = double(imbf2(I, :));
    idx = v > 0;
    
    % median intensity value in this row
    vmed = median(v(idx));
    
    % correct intensity in this row to match the target intensity
    imbf(I, idx) = double(imbf(I, idx)) / vmed * vtarget;
end

% the largest peak in the histogram of the wax area gives as the threshold
% background intensity
imbf(imbf <= mode(double(imbf(polymask)))) = 0;

% zero out the area outside the mask
imbf(~polymask) = 0;

% extend histogram to cover whole dynamic range
minbf = double(min(imbf(imbf > 0)));
maxbf = double(max(imbf(imbf > 0)));
imbf = uint8(255 * (double(imbf) - minbf) / (maxbf - minbf));

% match histology histogram to blockface's
idx = imh > 0;
imh(idx) = imhistmatch(imh(idx), imbf(ellipmask & (imbf > 0)));

% plot slices and histograms
if (DEBUG)
    subplot(2, 2, 1)
    imagesc(imh)
    subplot(2, 2, 2)
    hist(double(imh(imh>0)), 50)
    subplot(2, 2, 3)
    imagesc(imbf)
    subplot(2, 2, 4)
    hist(double(imbf(imbf>0)), 50)
end
