function [ T ] = vector2tensor( x )
% VECTOR2TENSOR converts a tensor from 6x1 or 51x1 format to 3x3 or 6x6
% format respectively.

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

if length(x) == 6
    T = [x(1), x(2), x(3);
         x(2), x(4), x(5);
         x(3), x(5), x(6)];
elseif length(x) == 15

    T = [   x( 1), x( 2), x( 3), x( 4), x( 5), x( 6);
            x( 2), x( 4), x( 5), x( 7), x( 8), x( 9);
            x( 3), x( 5), x( 6), x( 8), x( 9), x(10);
            x( 4), x( 7), x( 8), x(11), x(12), x(13);
            x( 5), x( 8), x( 9), x(12), x(13), x(14);
            x( 6), x( 9), x(10), x(13), x(14), x(15)];
end

end

