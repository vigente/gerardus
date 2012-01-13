function tips = correct_rv_crescent_tips(tips, valve, dthr)
% CORRECT_RV_CRESCENT_TIPS  Correct the segmentation of the Right
% Ventricle's crescent tips using the segmentation of the tricuspid and
% pulmonary annula
%
% TIPS = CORRECT_RV_CRESCENT_TIPS(TIPS, VALVE)
%
%   TIPS is a 3-column matrix with the coordinates of all the points in the
%   Right Ventricle crest, that corresponds to the crescent shape tips when
%   you cut the ventricle by an axial plane.
%
%   VALVE is a 3-colum matrix with the coordinates of either the tricuspid
%   valve or pulmonary valve segmetations.
%
%   The points in the segmentations are assumed to be consecutive.
%
% TIPS = CORRECT_RV_CRESCENT_TIPS(TIPS, VALVE, DTHR)
%
%   DTHR is the threshold distance value between the crest line and the
%   valve.

% Author: Ramon Casero <rcasero@gmail.com>
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
error( nargchk( 2, 3, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

% number of crest points
P = size(tips, 1);

% compute distance from each crest point to the valve
d = inf(P, 1);
for I = find(~isnan(sum(tips, 2)));
    d(I) = min(dmatrix(valve', tips(I,:)'));
end

% default value of dthr
if (nargin < 3 || isempty(dthr)) 
    dthr = min(d) * 5;
end

% find the first and last times the crest gets within a small distance of
% the valve
idx = find(d <= dthr);
if isempty(idx)
    error('The crest does not get close enough to the valve')
end
t1 = idx(1);
t2 = idx(end);

% % find one of the intersection points of the valve with the crest (crest
% % point index)
% [aux, t1] = min(d);
% 
% % we are going to ignore a bunch of points either side of the intersection
% % point
% L = round(size(valve, 1)/4);
% d2 = d;
% d2(t1-L:t1+L) = Inf;
% 
% % find second intersection point (crest point index)
% [aux, t2] = min(d2);
% aux = sort([t1 t2]);
% t1 = aux(1);
% t2 = aux(2);

% find intersection points (valve point index)
[aux, v1] = min(dmatrix(valve', tips(t1,:)'));
[aux, v2] = min(dmatrix(valve', tips(t2,:)'));

% find which half of the valve is closer to the centroid
aux = sort([v1 v2]);
v1 = aux(1);
v2 = aux(2);
if (mean(dmatrix(valve(v1:v2, :)', [0;0;0])) < ...
        mean(dmatrix(valve([v2:end 1:v1], :)', [0;0;0])))
    valveseg = valve(v1:v2, :);
else
    valveseg = valve([v2:end 1:v1], :);
end

% if the first point we find in the valve segment is closer to the _last_
% point we find in the crest gap, then we have to invert the segment's
% order
if (dmatrix(valveseg(1, :)', tips(t1,:)') > ...
        dmatrix(valveseg(end, :)', tips(t1,:)'))
    valveseg = valveseg(end:-1:1, :);
end

% replace the incorrect crest segment by the valve segment
tips = [tips(1:t1, :); valveseg; tips(t2:end, :)];
