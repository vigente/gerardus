function [thr, q, obj, im] = gmthr_seg(im, nobj, nsubs)
% GMTHR_SEG  Segment an image estimating threshold as intersection of two
% Gaussians from Gaussian mixture model
%
% THR = gmthr_seg(IM)
%
%   IM is an input n-dim image. This function assumes that IM contains a
%   darker object over a brighter background.
%
%   THR is a scalar with the estimated threshold value between dark
%   voxels (object) and lighter voxels (background). A Gaussian mixture
%   model is fitted to the image intensities, and the intersection point
%   between the Gaussian maxima is computed. The object in the image is
%   segmented using this intersection value as the segmentation threshold.
%
%   If the object and the background are too similar compared to the number
%   of samples in the image (i.e. the Gaussians intersect outside of the
%   interval between the Gaussian maxima), then this method cannot provide
%   a threshold to separate object and background. In that case, THR is
%   returned as NaN. This is the case, for example, if the image only
%   contains background, or only object voxels.
%
% [THR, Q, OBJ, BW] = gmthr_seg(IM, NOBJ, NSUBS)
%
%   NOBJ is a scalar. Only the largest NOBJ are kept in the segmentation.
%   This last step is useful to remove segmentation noise. By default,
%   NOBJ=1.
%
%   NSUBS is a scalar. For large images, the variance estimate will be too
%   small for the Gaussian mixture fitting function, which will return an
%   error. This problem can be solved doing a random subsampling of the
%   image to estimate the Gaussian mixture model. NSUBS is the subsampling
%   factor. E.g. NSUBS=100 will randomly sample numel(IM)/NSUBS voxels in
%   the image. By default, NSUBS=1 and no subsampling is performed.
%
%   Q is a quality measure of the threshold. Q takes values in [0, 1].
%   Values close to 0 mean that both Gaussians have a lot of overlap, so
%   the threshold between object and background cannot be trusted very
%   much. Values close to 1 mean that both Gaussians are well separated,
%   and the threshold value can be trusted to provide a good segmentation.
%
%   OBJ is the Gaussian mixture object. See help('gmdistribution.fit') for
%   details. The mean and variance of the Gaussians can be extracted as
%   obj.mu and obj.Sigma, respectively.
%
%   BW is an output segmentation mask, where voxels == true correspond to
%   the darker object.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012 University of Oxford
% Version: 0.5.1
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
nargoutchk(0, 4);

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

% % DEBUG: create Gaussian curves for display purposes
% ftis = normpdf(xout, mutis, sqrt(vartis));
% fbak = normpdf(xout, mubak, sqrt(varbak));

% compute intersection points between two gaussians
thr = intersect_gaussians(mutis, mubak, sqrt(vartis), sqrt(varbak));

% keep the one that is between both mean values
thr = thr(thr > mutis & thr < mubak);

% if there's no intersection point between the maxima, then returning a
% threshold is meaningless, and instead we return NaN. This is the case,
% e.g. if there's only background and no object
if isempty(thr)
    thr = nan;
end

% quality of the clustering measure. Integral under the tissue Gaussian
% in [thr, Inf] and integral under the background Gaussian in [-Inf, thr]:
% the sum represents the Gaussian overlap area. This overlap has a value in
% [0, 1], with 0 for a lot of overlap, and 1 for no overlap. The quality
% measure is then 1-overlap
if (nargout < 2)
    return
end
if (isnan(thr))
    q = nan;
else
    q = 1 - normcdf(2*mutis-thr, mutis, sqrt(vartis))...
        - normcdf(thr, mubak, sqrt(varbak));
end

% % DEBUG: plot histogram curves
% hold off
% plot(xout, ftot)
% hold on
% plot(xout, ftis, 'r')
% plot(xout, fbak, 'g')
% plot([thr, thr], [0 max([ftis(:); fbak(:)])], 'k')
% legend('all', 'tissue', 'background')
% xlabel('intensity')

% no need to waste time segmenting the image if the user doesn't ask for
% the output segmentation
if (nargout < 4)
    return
end

% threshold segmentation
im = (im <= thr);

%% remove segmentation noise
% we could use function bwrmsmallcomp() here, but we don't want having to
% replicate the segmentation data too many times. That could create memory
% problems for very large volumes. Instead, we have copied the code in
% bwrmsmallcomp() directly here:

% get connected components
cc = bwconncomp(im);

% keep only the largest components to remove background noise
%
% note: it's better to clear the whole image, and then add the largest
% components, than trying to delete the smaller components. The latter
% doesn't remove all the noise, for some reason.
len = cellfun(@length, cc.PixelIdxList);
[~, idx] = sort(len, 2, 'descend');
im = zeros(size(im), 'uint8');
if ~isempty(idx)
    idx = cc.PixelIdxList{idx(1:nobj)};
    im(idx) = 1;
end
