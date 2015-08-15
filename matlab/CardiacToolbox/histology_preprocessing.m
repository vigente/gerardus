function [imref, im, mask] = histology_preprocessing(imref, im)
% HISTOLOGY_PREPROCESSING  Prepare slices for intra-histology registration.
%
% HISTOLOGY_PREPROCESSING converts two histology images to grayscale,
% inverts and thresholds them (so that the background is black instead of
% white), extends the histograms to cover the dynamic range, and then
% matches the histograms. This prepares them to be registered.
%
% [IMREF2, IM2, MASK] = HISTOLOGY_PREPROCESSING(IMREF, IM)
%
%   IMREF, IM are two input histology images (in RGB colour or grayscale
%   format). The images can be provided as plain arrays (row, col,
%   channel), or as scimat format structs (see "help scimat" for details).
%   This function matches IM to IMREF.
%
%   IMREF2, IM2 are the output images after preprocessing.
%
%   MASK is a binary mask for IM2. This mask can be used to speed up
%   registration (when a mask is provided to elastix, only pixels within
%   the mask are used for the registration metric). The mask includes the
%   tissue and a bit of background around it.
%
% See also: histology_intraframe_reg.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014-2015 University of Oxford
% Version: 0.3.1
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

% check arguments
narginchk(2, 2);
nargoutchk(0, 3);

% save copy if inputs are provided as scimat struct
imrefIsStruct = isstruct(imref);
imIsStruct = isstruct(im);
if (imrefIsStruct)
    imref0 = imref;
    imref = squeeze(imref.data);
end
if (imIsStruct)
    im0 = im;
    im = squeeze(im.data);
end

% invert image intensities
[im, mask] = individual_image_preprocessing(im);
imref = individual_image_preprocessing(imref);

% match histograms to the central slice
for I = 1:size(im, 3)
    
    % select channels
    chref = imref(:, :, I);
    ch = im(:, :, I);
    
    % match histograms ignoring the background
    idxref = chref > 0;
    idx = ch > 0;
    ch(idx) = imhistmatch(ch(idx), chref(idxref));
    
    % replace image channel with processed ones
    imref(:, :, I) = chref;
    im(:, :, I) = ch;
    
end

% if inputs were scimat structs, recover metadata
if (imrefIsStruct)
    imref0.data = reshape(imref, size(imref, 1), size(imref, 2), 1, 1, size(imref, 3));
    imref = imref0;
end
if (imIsStruct)
    im0.data = reshape(im, size(im, 1), size(im, 2), 1, 1, size(im, 3));
    im = im0;
    im0.data = mask;
    mask = im0;
end

end

function [im, mask] = individual_image_preprocessing(im)

% loop image channels
for I = 1:size(im, 3)
    
    % select one channel
    ch = im(:, :, I);
    
    % the image should have two types of blackground pixels: white-ish and
    % pure black. The white-ish ones are from the original image, and the
    % black ones come from the B-spline fill-in. When we invert the image,
    % we don't want to invert the black ones, so we select non-black pixels
    idx = ch ~= 0;
    
    % invert image and make lowest intensity = 0
    ch(idx) = max(ch(idx)) - ch(idx);
    
    % remove background by zeroing anything <= mode (the background forms
    % the largest peak). We are quite conservative here in terms of keeping
    % all the tissue we can, to avoid losing detail, even if that means
    % that we are going to have a bit of background noise
    ch(ch <= mode(double(ch(idx)))) = 0;
    
    % extend histogram to cover whole dynamic range
    minh = double(min(ch(ch > 0)));
    maxh = double(max(ch(ch > 0)));
    ch = uint8(255 * (double(ch) - minh) / (maxh - minh));
    
    % replace image channel with the processed one
    im(:, :, I) = ch;
    
end

% compute a mask of the foreground
mask = uint8(rgb2gray(im) > 2);

% remove some background noise in the mask
se = strel('disk', 1);
mask = imerode(mask, se);

% dilate mask to cover some background
se = strel('disk', 10);
mask = imdilate(mask, se);

end
