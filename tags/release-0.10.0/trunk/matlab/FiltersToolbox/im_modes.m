function [mfg, mbg] = im_modes(im)
% IM_MODES  Estimate typical values (modes) of foreground/background
% intensities in a greyscale image.
%
% [MFG, MBG] = im_modes(IM)
%
%   IM is an N-array of grayscale intensity values.
%
%   MFG, MBG are scalars with the typical intensity values of foreground
%   (darker) and foreground (lighter) voxels. "Typical" means the mode of
%   the corresponding distributions.
%
%   The modes are estimated by computing histograms with up to 500 bins for
%   the intensity values, and smoothing the distribution until only 1 or 2
%   peaks are visible. In case of 1 peak, it is assume that the image
%   contains only background voxels, and MFG=NaN.

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
narginchk(1, 1);
nargoutchk(0, 2);

% ignore intensity values = 0. Those are considered to be masked out
im2 = im(im > 0);

% if the input image is empty, or too small, we assume that we don't have
% enough information to estimate the background and foreground
if (nnz(im2) < 100)
    mbg = nan;
    mfg = nan;
    return
end

% number of bins to use. We assume that we need at least 10 samples per
% bin, but we don't need more than 500 bins in total. Histograms with more
% bins are slower to process and smooth
nbin = min(500, ceil(nnz(im2)/10));

% compute histogram of the intensity values
[fhist, xhist] = hist(im2, nbin);
fhist = fhist / sum(fhist);

% initial estimation of peaks in the histogram
[pks, loc] = findpeaks(fhist, 'minpeakheight', 0.5e-3);

% number of peaks found
npks = length(pks);

% smooth the histogram until we find just 1 or 2 peaks
tol = 0.5e-10;
while ~((npks == 1) ...
        || ((npks == 2) && (abs(diff(loc)) >= 50)))
    
    % increase the smoothing parameter
    tol = tol * 2;
    
    % smooth the histogram
    [~, fhist2] = spaps(xhist, fhist, tol);
    
    % find the peaks
    [pks, loc] = findpeaks(fhist2, 'minpeakheight', 0.5e-3);
    
    % number of peaks found
    npks = length(pks);
    
end

% deal with the number of peaks
if (npks == 2)
    
    % we have background and tissue. The background is lighter. Sort the
    % peaks in darker to lighter order
    loc = sort(loc, 'ascend');

    % extract typical intensities for background and foreground
    mfg = xhist(loc(1));
    mbg = xhist(loc(2));
    
    % DEBUG
    if (DEBUG)
        subplot(2, 1, 1)
        hold off
        plot(xhist, fhist, 'b')
        hold on
        plot(xhist, fhist2, 'r')
        plot(mfg*[1 1], [0 max(fhist)], 'g')
        plot(mbg*[1 1], [0 max(fhist)], 'k')
        
        subplot(2, 1, 2)
        hold off
        imagesc(im(:, :, round((size(im, 3)+1)/2)))
    end

elseif (npks == 1)
    
    % we assume there's only background (although note that this could be a
    % case with only tissue)
    mfg = nan;
    mbg = xhist(loc);
    
    % DEBUG
    if (DEBUG)
        hold off
        plot(xhist, fhist, 'b')
        hold on
        plot(xhist, fhist2, 'r')
        plot(mbg*[1 1], [0 max(fhist)], 'r')
    end
    
else

    error('Assertion fail: The histogram has no peaks')
    
end

end
