function [ DT, FA, ADC, eVal, eVect, angle ] = DT_from_Smolydn_pt3( MRI_signal, gradient_directions, bvalue )
%DT_FROM_SMOLDYN_PT3 Calculates DT from MRI signal
%  
% PART 3 OF THE CALCULATION OF THE DIFFUSION TENSOR FROM THE SMOLDYN MODEL DATA
% This part calculates the diffusion tensor from the MRI signals.
%
% Inputs are:
%   MRI_signal: vector of no of directions, as calculated in part 2.
%   gradient_directions: matrix of size no of directions x 3 which gives
%       the x,y,z components of the gradient directions
%   bvalue: scalar. (um^2/us)
%
% Outputs are:
%   DT: diffusion tensor
%   FA: fractional anisotropy
%   ADC: apparant diffusion coefficient
%   eVal: eigenvalues of the diffusion tensor
%   eVect: eigenvectors of the diffusion tensor
%   angle: angle of primary eigenvector with the x? direction.

% Author: Jo Bates <jobates81@gmail.com>
% Copyright © 2014 University of Oxford
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

no_of_directions = size(gradient_directions,1);
H = zeros(no_of_directions,6);
for j = 1:no_of_directions
        h = [gradient_directions(j,1)^2, gradient_directions(j,2)^2, gradient_directions(j,3)^2, 2*gradient_directions(j,1)*gradient_directions(j,2), 2*gradient_directions(j,1)*gradient_directions(j,3), 2*gradient_directions(j,2)*gradient_directions(j,3)];
        H(j,:) = h;
end

% CONVERT DIFFUSION SIGNAL INTO TENSOR
Y = (-log(MRI_signal))/bvalue;
d = H\Y;

DT = [d(1), d(4), d(5); d(4), d(2), d(6); d(5), d(6), d(3)];

[eVect, eVal] = eig (DT);
eVal = diag(eVal); 
[eVal, order] = sort(eVal, 'descend');
eVect = eVect(:,order);
eVal = abs(eVal);
ADC = mean(eVal);
FA = sqrt(1.5)*sqrt(((eVal(1)-ADC)^2+(eVal(2)-ADC)^2+(eVal(3)-ADC)^2)/(eVal(1)^2+eVal(2)^2+eVal(3)^2));
angle = atand(eVect(2,1)/eVect(1,1));
       
end

