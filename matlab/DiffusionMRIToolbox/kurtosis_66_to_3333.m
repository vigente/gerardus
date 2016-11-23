function W = kurtosis_66_to_3333(K)    
%KURTOSIS_66_TO_3333 Converts a 6x6 kurtosis tensor to 3x3x3x3 format
%   K is the 6x6 kurtosis tensor
%   W is the 3x3x3x3 tensor

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


% K_index =[01, 02, 03, 04, 05, 06; % xx
%           02, 04, 05, 07, 08, 09; % xy
%           03, 05, 06, 08, 09, 10; % xz
%           04, 07, 08, 11, 12, 13; % yy
%           05, 08, 09, 12, 13, 14; % yz
%           06, 09, 10, 13, 14, 15];% zz
%         % xx  xy  xz  yy  yz  zz              
                                
W = zeros([3 3 3 3]);

for i = 1:3
    for j = 1:3
         for k = 1:3
             for l = 1:3
                [r, c] = kurtosis_ijkl_to_rc(i, j, k, l);

                % if you have it as a 15x1 vector, use this
                %K2_index = K_index(r,c);

                W(i,j,k,l) = K(r,c);
             end
         end
    end
end