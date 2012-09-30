function [im, thr] = gmthr_seg(im, nobj, nsubs)
% GMTHR_SEG  Segment an image estimating threshold as intersection of two
% Gaussians from Gaussian mixture model
%
% [BW, THR] = gmthr_seg(IM)
%
%   IM is an input n-dim image. This function assumes that IM contains a
%   darker object over a brighter background.
%
%   BW is an output segmentation mask, where voxels == true correspond to
%   the darker object.
%
%   THR is a scalar with the estimated threshold value.
%
%   A Gaussian mixture model is fitted to the image intensities, and the
%   intersection point between the Gaussian maxima is computed. The object
%   in the image is segmented using this intersection value as the
%   segmentation threshold. Finally, all segmented objects are removed
%   except for the largest one. This last step is useful to remove
%   segmentation noise.
%
% ... = gmthr_seg(..., NOBJ, NSUBS)
%
%   NOBJ is a scalar. Only the largest NOBJ are kept in the segmentation.
%   By default, NOBJ = 1.
%
%   NSUBS is a scalar. For large images, the variance estimate will be too
%   small for the Gaussian mixture fitting function, which will return an
%   error. This problem can be solved doing a random subsampling of the
%   image to estimate the Gaussian mixture model. NSUBS is the subsampling
%   factor. E.g. NSUBS=100 will randomly sample numel(IM)/NSUBS voxels in
%   the image. By default, NSUBS=1 and no subsampling is performed.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012 University of Oxford
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
narginchk(1, 3);
nargoutchk(0, 1);

% defaults
if (nargin < 2 || isempty(nobj))
    nobj = 1;
end
if (nargin < 3 || isempty(nsubs))
    nsubs = 1;
end

% % DEBUG: approximate pdf of whole image
% [ftot, xout] = hist(im(:), 100);
% inc = xout(2) - xout(1);
% ftot = ftot / numel(im) / inc;

% compute gaussian mixture model. Note that we need to randomly subsample
% the image so that the variance of the Gaussians is not too small
% (otherwise, gmdistribution.fit() gives an error)
if (nsubs > 1)
    idx = randi(numel(im), round(numel(im)/nsubs), 1);
    obj = gmdistribution.fit(im(idx), 2);
else
    obj = gmdistribution.fit(im(:), 2);
end

% get mixture of Gaussians parameters
[mutis, idx] = min(obj.mu);
vartis = obj.Sigma(idx);
[mubak, idx] = max(obj.mu);
varbak = obj.Sigma(idx);

% compute intersection points between two gaussians
thr = intersect_gaussians(mutis, mubak, sqrt(vartis), sqrt(varbak));

% keep the one that is between both mean values
thr = thr(thr > mutis & thr < mubak);

% % DEBUG: create Gaussian curves for display purposes
% ftis = normpdf(xout, mutis, sqrt(vartis));
% fbak = normpdf(xout, mubak, sqrt(varbak));

% % DEBUG: plot histogram curves
% hold off
% plot(xout, ftot)
% hold on
% plot(xout, ftis, 'r')
% plot(xout, fbak, 'g')
% plot([thr, thr], [0 max([ftis(:); fbak(:)])], 'k')
% legend('all', 'tissue', 'background')
% xlabel('intensity')

% threshold segmentation
im = (im <= thr);

% get connected components
cc = bwconncomp(im);

% keep only the largest component to remove background noise
%
% note: it's better to clear the whole image, and then add the largest
% components, than trying to delete the smaller components. The latter
% doesn't remove all the noise, for some reason.
len = cellfun(@length, cc.PixelIdxList);
[~, idx] = sort(len, 2, 'descend');
im = zeros(size(im), 'uint8');
idx = cc.PixelIdxList{idx(1:nobj)};
im(idx) = 1;
