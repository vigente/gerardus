function [x, y, xidx, yidx] = match_file_lists(x, y, xnameidx, ynameidx)
% match_file_lists  Remove elements from a file list that are not in
% another.
%
% [X2, Y2, XIDX, YIDX] = match_file_lists(X, Y, XNAMEIDX, YNAMEIDX)
%
%   X, Y are file lists obtained with dir(). X, Y are arrays of structs
%   with a field X.name, Y.name for the filenames.
%
%   XNAMEIDX, YNAMEIDX are indices that determine which part of the
%   filenames are matched.
%
%   X2, Y2 are the returned lists. X2 is X minus the elements not found in
%   Y. Y2 is Y minus the elements not found in X.
%
%   XIDX, YIDX are the indices of the matched elements. X2 = X(XIDX),
%   Y2 = Y(YIDX).
%
%   For example,
%
%   X(1).name = 'Q53_55_0004.png';
%
%   Y(1).name = 'Q53_002 - 2014-06-24 16.43.14.png';
%   Y(2).name = 'Q53_004 - 2014-06-24 16.44.51.png';
%   Y(3).name = 'Q53_006 - 2014-06-24 16.46.35.png';
%
%   XNAMEIDX = 9:11;
%   YNAMEIDX = 5:7;
%
%   This function matches X(1).name(9:11) = '004' to the other three
%   elements, '002', '004' and '006'. Because there's a match with the
%   second element in Y, the function returns
%
%   X = X(1).name;
%   Y = Y(2).name;

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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
narginchk(2, 4);
nargoutchk(0, 4);

% defaults
if (nargin < 3)
    xnameidx = [];
end
if (nargin < 4)
    ynameidx = [];
end

% convert list of filenames to char array
xs = cat(1, x.name);
ys = cat(1, y.name);

% crop to the part of the filename that we have to match
if (~isempty(xs))
    xs = xs(:, xnameidx);
end
if (~isempty(ys))
    ys = ys(:, ynameidx);
end

% find elements of first list that are also in the second list
xidx = ismember(xs, ys, 'rows');

% find elements of second list that are also in the first list
yidx = ismember(ys, xs, 'rows');

% remove elements that are not in both lists
x = x(xidx);
y = y(yidx);
