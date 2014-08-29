function [ellipmask, polymask, immean] = blockface_create_masks(indir, files)
% blockface_create_masks  User interface to create the elliptical and
% polygonal segmentation masks used to correct blockface illumination.
%
% blockface_create_mask shows an average of the whole blockface volume, and
% prompts the user to draw an ellipse tightly around the heart. Then, it
% prompts the user to draw a polygon around the block's edge. Both drawings
% are converted to binary segmentation masks. Those masks can then be used
% in blockface_equalise_illumination() to correct the blockface
% illumination.
%
% [ELLIPMASK, POLYMASK] = blockface_create_masks(INDIR, FILES)
%
%   INDIR is a string with the directory where the input files are kept.
%
%   FILES is the result of a dir() command with a list of the files to
%   correct.
%
%   ELLIPMASK, POLYMASK are binary segmentations of the elliptical and
%   polygonal masks drawn by the user.
%
% [..., IMMEAN] = blockface_create_masks(FILE_EXPR)
%
%   IMMEAN is the average blockface intensity image shown to the user. This
%   can be useful for debugging purposes.

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

% check arguments
narginchk(2, 2);
nargoutchk(0, 3);

% number of files
N = length(files);

if (N == 0)
    ellipmask = [];
    polymask = [];
    return
end

% load the first image, so that we know the number of rows and columns in
% the blockface volume
im = imread([indir filesep files(1).name]);

% initialize matrix to keep the average image
immean = zeros(size(im, 1), size(im, 2));

% compute an average of the slice intensities, that the user can employ as
% a guide to draw the masks
for I = 1:N
    
    % load image from file
    im = imread([indir filesep files(I).name]);
    
    % convert to grayscale if image is in colour
    if (size(im, 3) == 3)
        im = rgb2gray(im);
    end
    
    % add to the average
    immean = immean + double(im) / N;
    
end

% plot the average intensity image
him = imagesc(immean);
title('Draw an ellipse tightly around the heart. Double-click on ellipse to finish.')

% let the user draw an elliptical mask around the heart
hmask = imellipse;
wait(hmask);

% convert the drawing to a binary segmentation
ellipmask = createMask(hmask, him);

% let the user draw a polygonal mask around the block
title('Draw a polygon on the block''s edge. Double-click on polygon to finish.')
hmask = impoly;
wait(hmask);

% convert the drawing to a binary segmentation
polymask = createMask(hmask, him);
