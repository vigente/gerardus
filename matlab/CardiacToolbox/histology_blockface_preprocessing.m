function [imh, imbf] = histology_blockface_preprocessing(imh, imbf, maskbf)
% HISTOLOGY_BLOCKFACE_PREPROCESSING  Prepare histology and blockface for
% registration.
%
% histology_blockface_preprocessing takes a histology and blockface images
% that are going to be registered onto each other, and: 1) inverts them, to
% have a dark background, 2) thresholds the background, to remove noise, 3)
% extends the intensity range and 4) matches the histogram of the histology
% image to the blockface.
%
% [IMH2, IMBF2] = histology_blockface_preprocessing(IMH, IMBF)
%
%   IMH, IMBF are colour or grayscale images with the histology and
%   blockface, respectively.
%
%   IMH2, IMBF2 are the corresponding images after preprocessing.
%
% ... = histology_blockface_preprocessing(..., MASKBF)
%
%   MASKBF is a binary segmentation with the same size as IMBF. MASKF
%   contains a binary mask. IMBF will be cropped to the smallest rectangle
%   that contains that binary mask.

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

DEBUG = 0;

% check arguments
narginchk(2, 3);
nargoutchk(0, 2);

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

% invert image 
imh = max(imh(:)) - imh;

% find threshold between background and tissue
thr = gmthr_seg(double(imh));

% remove background
imh(imh <= thr) = 0;

% normalize histogram
imh = uint8(double(imh) / double(max(imh(:))) * 255);

% plot preprocessed histology
if (DEBUG)
    subplot(2, 1, 1)
    imagesc(imh)
end

%% Blockface preprocessing

% if a blockface mask is provided, then we crop the blockface image to the
% smallest rectangle that contains the mask
if (nargin > 2)
    if (any(size(maskbf) ~= size(imbf)))
        error('Blockface mask must be the same size as blockface image')
    end
    
    % top/bottom limits of the mask
    idx = find(sum(maskbf, 2) > 0);
    top = idx(1);
    bot = idx(end);
    
    % left/right limits of the mask
    idx = find(sum(maskbf, 1) > 0);
    lef = idx(1);
    rig = idx(end);
    
    % crop blockface image
    imbf = imbf(top:bot, lef:rig);
    
end

% invert image
imbf = max(imbf(:)) - imbf;

% find threshold between background and tissue
thr = gmthr_seg(double(imbf(:)));

% threshold background
imbf(imbf < thr) = 0;

% plot preprocessed blockface
if (DEBUG)
    subplot(2, 1, 2)
    imagesc(imbf)
end

% fill the range of blockface intensities
idx = imbf>0;
imin = min(imbf(idx));
imax = max(imbf(idx));
imbf(:) = double(imbf(:) - imin(:)) / double(imax - imin) ...
    * double(intmax(class(imbf)));

% match histology histogram to blockface's
idx = imh > 0;
imh(idx) = imhistmatch(imh(idx), imbf(imbf>0));

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
