function [r, c] = kurtosis_ijkl_to_rc(i, j, k, l)
%KURTOSIS_IJKL_TO_RC converts indices in a kurtosis tensor from 3x3x3x3
%format to 6x6 format

% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright © 2014-2016 University of Oxford
% Version: 0.1.0
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

switch i
    case 1
        r1 = [1 2 3];
    case 2
        r1 = [2 4 5];
    case 3
        r1 = [3 5 6];
end
switch j
    case 1
        r2 = [1 2 3];
    case 2
        r2 = [2 4 5];
    case 3
        r2 = [3 5 6];
end

if i == j
    switch i
        case 1
            r =1;
        case 2
            r = 4;
        case 3
            r = 6;
    end
else
    r = intersect(r1, r2);
end


switch k
    case 1
        c1 = [1 2 3];
    case 2
        c1 = [2 4 5];
    case 3
        c1 = [3 5 6];
end
switch l
    case 1
        c2 = [1 2 3];
    case 2
        c2 = [2 4 5];
    case 3
        c2 = [3 5 6];
end

if k == l
    switch k
        case 1
            c =1;
        case 2
            c = 4;
        case 3
            c = 6;
    end
else
    c = intersect(c1, c2);
end