function [ x ] = tensor2vector( T, scale )
% TENSOR2VECTOR converts a tensor from 3x3 or 6x6 format to 6x1 or 15x1
% format respectively.
% The scale parameter tells you whether to account for the fact that
% certain elements are repeated - you might use this if you want to
% multiply tensors together

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

if nargin == 1
    scale = 0;
end

if size(T,1) == 3
    x = [T(1),T(2),T(3),T(5),T(6),T(9)];
    
    if scale
        x = x .* [1 2 2 1 2 1];
    end
    
elseif size(T,1) == 6
    x = [T(1), T( 2), T( 3), T( 4), T( 5), T( 6), T(10), T(11), T(12), ...
        T(18), T(22), T(23), T(24), T(30), T(36)];
    
    if scale
        x = x .* [1 2 2 3 4 3 2 4 4 2 1 2 3 2 1];
    end
    
end


end

