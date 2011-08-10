function y = crest_envelope(x, c, ITER, P)
% CREST_ENVELOPE  Compute the envelope for the Right Ventricle's crest
%
% Y = CREST_ENVELOPE(X, C)
%
%   X is a 3-column matrix with the coordinates of the Right Ventricle's
%   crest points.
%
%   C is a 3-column matrix with the coordinates of the Left Ventricle's
%   centroid curve.
%
%   Y is a 3-column matrix with the coordinates of the crest envelope.
%
% Y = CREST_ENVELOPE(X, C, ITER, P)
%
%   ITER is a scalar with the number of times the distance curve is
%   iterated to remove local maxima. By default, ITER=4.
%
%   P is a scalar with the interpolation factor. P=0 means perfect
%   smoothing (least squares straight line fitting), P=1 means perfect
%   interpolation (no smoothing at all). That is, the smaller the P, the
%   bigger the smoothing. By default, P=0.99.

% Author: Ramón Casero <ramon.casero@comlab.ox.ac.uk>
% Copyright © 2010 University of Oxford
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
error(nargchk(2, 4, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% defaults
if (nargin < 3 || isempty(ITER))
    ITER=4;
end
if (nargin < 4 || isempty(P))
    P=.99;
end


% remove NaNs
x = x(~isnan(sum(x, 2)), :);

% number of points
N = size(x, 1);

% compute distance from each crest point to each centroid curve point
d = dmatrix(x', c');

% keep only the shortest distance from each crest point
% this also provides a correspondence from each crest point to a centroid
[d, cidx] = min(d, [], 2);

% number of elements to repeat at the beginning and end of the vector, so
% that it is considered a cyclic vector
PAD = 200;

% pad at the end with the beginning of the vector to make it cyclic
d = d([(end-PAD+1):end 1:end 1:PAD]);

% find local minima
idx = localminima(d);

% interpolate and resample to have a vector from the local minima of the
% same length as the distance vector. We use linear interpolation instead
% of spline because the latter originates ringing
dmin = interp1(idx,d(idx),1:length(d),'linear', 'extrap')';

% even linear interpolation can in some cases produce interpolated distance
% values larger than the original distance values. We make sure that in
% those cases we use the original curve instead of the 
idxmin = find((d - dmin) < 0);
dmin(idxmin) = d(idxmin);

% find local minima in the equalized signal
idx = localminima(dmin);

% for 3 consecutive minima, if they have a "^" shape, remove the minimum in
% the middle
for I = 1:ITER
    idx = removewaves(idx, dmin);
end

% interpolate
dmin = interp1(idx,dmin(idx),1:length(dmin),'cubic', 'extrap')';

% chop off the padding
d = d(PAD+1:end-PAD);
dmin = dmin(PAD+1:end-PAD);

% direction vector from centroid to crest
v = x - c(cidx, :);

% compute points of the envelope. To do this, we find an intermediate point
% at dmin along the straigt line connecting the centroid and the crest
% point
k = dmin./d;
y = c(cidx, :) + k(:, ones(1, size(x, 2))) .* v;

% make envelope cyclic
y = y([1:end 1], :);

% compute Lee's centripetal knot points
t = cumsum([0;((diff(y).^2)...
    *ones(size(y, 2),1)).^(1/4)]).';

% compute smoothing cubic spline for each coordinate
ppx = csaps(t, y(:, 1), P);
ppy = csaps(t, y(:, 2), P);
ppz = csaps(t, y(:, 3), P);

% sample splines
y = [ppval(ppx, t)' ppval(ppy, t)' ppval(ppz, t)'];

end

function idx = removewaves(idx, x)
% for 3 consecutive minima, if they have a "^" shape, remove the minimum in
% the middle

flag = true(size(idx));
for I = 2:length(idx)-1
    flag(I) = (x(idx(I-1)) > x(idx(I)) || ...
        x(idx(I+1)) > x(idx(I)));
end
idx = idx(flag);
end

function idx = localminima(x)

% sign of the slope of the crest heights
signslope = sign(diff(x));

% find points that have a negative slope followed by a positive slope
idx = find([signslope(1)>0; diff(signslope)>0; signslope(end)<0]);

end
