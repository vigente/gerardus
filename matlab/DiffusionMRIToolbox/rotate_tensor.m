function [ OUT ] = rotate_tensor( IN, vec )
%ROTATE_TENSOR Rotates a tensor to a new reference frame
%   IN is the input tensor (can be 3x3 for a diffusion tensor, or 6x6 for a
%   kurtosis tensor)
%   vec are the eigenvectors defining a new reference frame (must be 3x3)

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

if size(IN,1) == 3
    OUT = pinv(vec) * IN * vec;
elseif size(IN,1) == 6
    
    % this comes from dipy: https://github.com/nipy/dipy/blob/330f95b1af1c1d120394c5c2b074310ad278167d/dipy/reconst/dki.py#L1536
    
    
    W = kurtosis_66_to_3333(IN);
    
    W_rotated = W * 0;
    for W_i = 1:3 % for each element in the kurtosis tensor
        for W_j = 1:3
            for W_k = 1:3
                for W_l = 1:3
                    for i = 1:3 % get the sum of the contributions from the 4 dimensions
                        for j = 1:3
                            for k = 1:3
                                for l = 1:3
                                    W_rotated(W_i, W_j, W_k, W_l) = W_rotated(W_i, W_j, W_k, W_l) + vec(i, W_i) * vec(j, W_j) * vec(k, W_k) * vec(l, W_l) * W(i,j,k,l);
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    OUT = kurtosis_3333_to_66(W_rotated);
    
end 
    
end

