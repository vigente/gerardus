function [imref, im, mask] = histology_preprocessing(imref, im, opts)
% HISTOLOGY_PREPROCESSING  Prepare slices for intra-histology registration.
%
% HISTOLOGY_PREPROCESSING converts a list of histology images to grayscale,
% inverts and thresholds them (so that the background is black instead of
% white), extends the histograms to cover the dynamic range, and then
% matches the histograms to a reference slice. This prepares the slices for
% successful intra-histology registration using transfdiffreg().
%
% [IMREF2, IM2, MASK] = HISTOLOGY_PREPROCESSING(IMREF, IM)
%
%   IMREF is the slice whose histogram will be used as reference for the
%   others. IM is the stack of RGB colour or grayscale histology images to
%   preprocess. The images can be provided as:
%
%     * A char with a single path and filename. E.g. 'file.mha'. With this
%       format, an output directory must be provided in OPTS.outdir, or
%       else a temporal one with a random name is created (see below).
%
%     * A cell array of filenames. E.g. {'file1.mha', 'file2.mha'}. IM{I}
%       is the path and filename of the I-th slice. See note about
%       OPTS.outdir in previous item.
%
%     * A plain array IM(R, C, I, CH), R=row, C=column, I=slice,
%       CH=channel.
%
%     * A vector of scimat structs. IM(I) is a scimat struct with the I-th
%       slice. See "help scimat" for details.
%
%   IMREF2, IM2 are the preprocessed images in the same format as the
%   inputs.
%
%   MASK is a binary mask for IM2. This mask can be used to speed up
%   registration (when a mask is provided to elastix, only pixels within
%   the mask are used for the registration metric). The mask includes the
%   tissue and a bit of background around it.
%
% ... = HISTOLOGY_PREPROCESSING(..., OPTS)
%
%   OPTS is a struct with options for the algorithm.
%
%     'outdir': If IM is a list of filenames, this is the path where the
%               preprocessed images are saved to.
%
%     'Grayscale': (def false) true/false whether image will be converted
%              to grayscale.
%
%
% See also: transfdiffreg.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014-2015 University of Oxford
% Version: 0.4.0
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
narginchk(2, 3);
nargoutchk(0, 3);

% defaults
if (nargin < 3 || isempty(opts))
    opts = struct;
end
if (~isfield(opts, 'Grayscale') || isempty(opts.Grayscale))
    opts.Grayscale = false;
end

% check the type and number of the input images
if (ischar(im)) % single filename

    imType = 'filename';
    im = {im};
    N = 1;
    
elseif (iscell(im)) % list of filenames

    imType = 'filename';
    N = length(im);
    
elseif (isstruct(im)) % list of scimat structs
    
    imType = 'scimat';
    N = length(im);
    
else % plain array
    
    imType = 'array';
    N = size(im, 3);
    
end

% convert IMREF to a plain array if it's not already
if (ischar(imref)) % filename

    imref = scimat_load(imref);

end
if (isstruct(imref)) % scimat struct
    
    imref = imref.data;
    
else
    
    % move channels to 5th index, so that it's the same format as IM.data
    imref = reshape(imref, size(imref, 1), size(imref, 2), 1, 1, size(imref, 3));
    
end

% grayscale conversion of IMREF, if requested by user
if (opts.Grayscale)
    
    imref = rgb2gray(squeeze(imref));
    
end

% preprocessing of IMREF
imref = individual_image_preprocessing(imref);

% if images are provided as filenames, either the user must provide an
% output directory, or we generate a temp one with a random name
if (strcmp(imType, 'filename') ...
        && (~isfield(opts, 'outdir') || isempty(opts.outdir)))
    opts.outdir = tempname;
end

% allocate memory for masks
switch (imType)
    
    case  'filename' % filename
        
        mask = cell(1, N);
        
    case 'scimat' % scimat struct
        
        mask = im;
        
    case 'array' % plain array
        
        mask = zeros(size(im, 1), size(im, 2), size(im, 3), 'uint8');
        
end

% loop all images
for I = 1:N

    % convert IM to a scimat struct if it's not already
    switch (imType)
        
        case  'filename' % filename
        
            im0 = scimat_load(im{I});
            
        case 'scimat' % scimat struct
            
            im0 = im(I);
            
        case 'array' % plain array
            
            im0 = scimat_im2scimat(reshape(im(:, :, I, :, :), ...
                [size(im, 1) size(im, 2) 1 1 size(im, 4)]));
        
    end
    
    % grayscale conversion of IM, if requested by user
    if (opts.Grayscale)
        
        im0.data = rgb2gray(squeeze(im0.data));
        
    end

    % preprocessing of IMREF
    mask0 = im0;
    [im0.data, mask0.data] = individual_image_preprocessing(im0.data);

    % match IM histogram to IMREF
    im0.data = match_histograms(imref, im0.data);
    
    % convert processed images to output format
    switch (imType)
        
        case  'filename' % filename
            
            [~, name, ext] = fileparts(im{I});
            im{I} = [opts.outdir filesep name ext];
            scimat_save(im{I}, im0);
            mask{I} = [opts.outdir filesep 'mask-' name ext];
            scimat_save(mask{I}, mask0);
            
        case 'scimat' % scimat struct
            
            im(I) = im0;
            mask(I) = mask0;
            
        case 'array' % plain array
            
            im(:, :, I, :) = reshape(im0.data, ...
                [size(im, 1) size(im, 2) 1 size(im, 5)]);
            mask(:, :, I, :) = mask0.data;
            
    end
    
end

end

%% nested functions

% individual_image_preprocessing:
%
%   Preprocessing of an individual image. Invert, threshold background,
%   extend histogram to whole dynamic range. Compute a rough binary mask of
%   the object.
function [im, mask] = individual_image_preprocessing(im)

% loop image channels
for I = 1:size(im, 5)
    
    % select one channel
    ch = im(:, :, :, :, I);
    
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
    im(:, :, :, :, I) = ch;
    
end

% compute a mask of the foreground
if (size(im, 5) == 1)
    
    mask = uint8(im > 2);
    
else
    
    mask = uint8(rgb2gray(squeeze(im)) > 2);
    
end

% remove some background noise in the mask
se = strel('disk', 1);
mask = imerode(mask, se);

% dilate mask to cover some background
se = strel('disk', 10);
mask = imdilate(mask, se);

end

% match_histograms:
%
%   Match the histogram of each channel in IM to IMREF, ignoring background
%   voxels.
%
%   IMREF: reference image, array (R, C, CH)
function im = match_histograms(imref, im)

% match histograms to the central slice
for I = 1:size(im, 5)
    
    % select channels
    chref = imref(:, :, :, :, I);
    ch = im(:, :, :, :, I);
    
    % match histograms ignoring the background
    idxref = chref > 0;
    idx = ch > 0;
    ch(idx) = imhistmatch(ch(idx), chref(idxref));
    
    % replace image channel with processed ones
    %imref(:, :, I) = chref;
    im(:, :, :, :, I) = ch;
    
end

end
