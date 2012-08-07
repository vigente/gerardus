function ispeak = findpeaksn(im, usedim, minheight, mindist)
% FINDPEAKSN  Find local peaks (local maxima) in n-dimensional data
%
% ISPEAK = findpeaks(IM)
%
%   IM is an n-dimensional array, e.g. a vector, 2D image or 3D volume.
%
%   ISPEAK is a boolean array of the same size. True elements correspond to
%   local peaks (local maxima) in IM. A local peak must be a peak in each
%   of the n dimensions.
%
%   Note: to mirror Matlab's function findpeaks(), we consider as peaks
%   only the end of increasing slopes. So if we have a raising slope, a
%   flat region, and then a falling slope, we only consider one peak (the
%   end of the raising slope).
%
% ISPEAK = findpeaks(IM, USEDIM, MINHEIGHT, MINDIST)
%
%   USEDIM is a boolean vector with length n. Each element says whether the
%   corresponding dimension should be taken into account to define the
%   peaks. By default, USEDIM=true(1, ndims(IM)) and all dimensions are
%   taken into account.
%
%   MINHEIGHT is a scalar with the minimum value an element in IM must have
%   to be considered a peak. This can be used to reduce false positives
%   given by background noise. By default, MINHEIGHT=0 and all peaks are
%   accepted.
%
%   MINDIST is a scalar with the minimum distance (in pixel/voxel units)
%   between any two peaks. Smaller peaks too close to larger peaks can be
%   removed using this paramenter.
%
% See also: findpeaks.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012 University of Oxford
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
error(nargchk(1, 4, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% defaults
if (nargin < 2 || isempty(usedim))
    usedim = true(1, ndims(im));
end
if (nargin < 3) || isempty(minheight)
    minheight = 0;
end
if (nargin < 4) || isempty(mindist)
    mindist = 0;
end

% init boolean volume with the same size as the image. True elements
% indicate where a peak is
ispeak = true(size(im));

% a local peak in n-dimensions has to be a peak in each of the n dimensions
% individually
for D = 1:ndims(im)

    % skip current dimension if the user doesn't want it to count towards
    % the peak
    if ~usedim(D)
        continue
    end
    
    %% find all peaks
    
    % shift the dimensions so that the current dimension is along the rows.
    % This makes then writing the code easier
    [im, perm, nshifts] = shiftdata(im, D);
    
    % assume that the background beyond the image is always 0, so add a row
    % of zeros along the current dimension
    sz = size(im);
    sz(1) = 1;
    aux = cat(1, im, zeros(sz));
    
    % first we sweep the dimension in the ascending order
    aux = sign(diff(aux, 1));
    
    % remove peak candidates that are not peaks in this dimension
    %
    % to mirror Matlab's findpeaks(), we consider peaks only the end of
    % increasing slopes. So if we have a raising slope, a flat, and then a
    % falling slope, we only consider one peak (the end of the raising
    % slope)
    %
    % in the first sweep, peaks must have aux <= 0
    ispeak = ispeak & unshiftdata(aux<=0, perm, nshifts);

    % second we sweep rows in the descending order, with the zero-padding
    % as the first row
    aux = cat(1, zeros(sz), im);
    aux = flipdim(sign(diff(flipdim(aux, 1), 1)), 1);
    
    % remove peak candidates that are not peaks in this dimension
    %
    % to mirror Matlab's findpeaks(), we consider peaks only the end of
    % increasing slopes. So if we have a raising slope, a flat, and then a
    % falling slope, we only consider one peak (the end of the raising
    % slope)
    %
    % in the second sweep, peaks must have aux == -1
    ispeak = ispeak & unshiftdata(aux==-1, perm, nshifts);
    
    % undo the shift of dimensions
    im = unshiftdata(im, perm, nshifts);
    
end

%% remove peaks that are too small

ispeak(im < minheight) = false;
    
%% remove smaller peaks that are too close to other bigger peaks

if (mindist > 0) % we want to remove smaller peaks that are too close to bigger peaks
    
    % get list of peaks in the image, and their values
    idx = find(ispeak);
    val = im(idx);
    
    % sort peaks according to value in descending order
    [val, subidx] = sort(val, 'descend');
    idx = idx(subidx);
    
    % coordinates of peaks
    [sub{1:ndims(im)}] = ind2sub(size(im), idx); % each cell has 1 coordinate of all points
    sub = cell2mat(sub);
    sub = mat2cell(sub, ones(1, length(val)), ndims(im)); % each cell has all coordinates of 1 point

    % initialise vector to note which peaks need to be deleted
    todelete = false(size(val));

    % loop every peak
    for I = 1:length(val)
        if ~todelete(I) % we have found a peak that doesn't have to be removed
            
            % compute the distance of the other peaks to this one
            d = cellfun(@(x) norm(x - sub{I}), sub);
            
            % peaks that are too close are tagged to be removed, except the
            % current peak itself
            todelete = todelete | (d < mindist);
            todelete(I) = false;
            
        end
    end
    
    % remove smaller peaks that are too close
    ispeak(idx(todelete)) = 0;

end

% % return peaks as list of indices to save memory
% ispeak = find(ispeak);

end
