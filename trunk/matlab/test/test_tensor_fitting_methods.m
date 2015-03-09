% test_tensor_fitting_nethods.m

% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright © 2015 University of Oxford
% Version: 0.1.1
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.


% This script generates a signal with a known FA and ADC, and compares two
% different methods for tensor fitting. The first is fit_DT, which fits the
% 6 tensor elements, as well as a 7th parameter to handle S0. The second
% method is fit_DT_YHd_method, which requires an input of S0 from
% un-attenuated images, which in practice are never exactly zero. 

% The other difference between the scripts is the b value input arguments.
% fit_DT requires a 3x3 b matrix, whereas fit_DT_YHd_method requires only a
% b value and unit vector. Cross-terms aren't handled in the latter case -
% it is assumed that the diffusion gradient is applied only in the
% direction of the unit vector.

% fit_DT also has some extra functionality, such as being able to handle
% any size input images, doing a weighted linear fit, and doing a
% non-linear fit on voxels within a mask. It can also fit a positve
% constrained linear or non-linear model.



close all
clear all
clc


% start with 2 b=15 (approximates b=0) and 6x b = 600

dirdata = [0 0 1;
           0 0 1;
          1/sqrt(2), 0, -1/sqrt(2);
         -1/sqrt(2), 0, -1/sqrt(2);
         0, -1/sqrt(2), -1/sqrt(2);
         0, -1/sqrt(2), 1/sqrt(2);
         1/sqrt(2), -1/sqrt(2), 0;
        -1/sqrt(2), -1/sqrt(2), 0];

% b values - when the b=0 values are actually zero, and there is no 
% noise, the two methods are the same
b = [15 15 600 600 600 600 600 600]';


% convert the b value and unit vectors into b matrices
bval = zeros(3,3,8);

for n = 1:size(dirdata,1)
    bval(:,:,n)=b(n)*dirdata(n,:)'*dirdata(n,:);
end


% generate the signal:

FA = 0.25;
ADC = 1.2E-3;

disp('Aim for:')
disp(['ADC = ' num2str(ADC)])
disp(['FA = ' num2str(FA)])


v1 = 1; % this can be anything

% get v2 (and v3) based on the FA
v2 = v1 * (sqrt(1-(2*FA^2-1)*(FA^2-1)) - 1) / (2*FA^2 - 1);

D = [v1, 0 0; 0  v2 0; 0 0  v2];

% scale to get the correct ADC
D = D ./ sum(D(:)) * ADC * 3;

% rotate from this orientation:
v1 = [1 0 0];

% to this orientation: (arbitrarily chosen)
v2 = [1/sqrt(2) 1/sqrt(2) 0];  

R = vrrotvec2mat(vrrotvec(v1, v2));

% rotate the tensor
D_rot = R * D * R';

% convert to a vector (diagonally symmetric)
M = [D_rot(1,1), D_rot(1,2), D_rot(1,3), D_rot(2,2), D_rot(2,3), D_rot(3,3)]';


Bv=squeeze([bval(1,1,:),2*bval(1,2,:),2*bval(1,3,:),bval(2,2,:),2*bval(2,3,:),bval(3,3,:)])';

% get the signal for this tensor
S = zeros(1, size(bval, 3));
for i = 1:length(S)

    Bv_i = Bv(i,:);

    S(i) = exp(-Bv_i * M);

end


% add noise (Rician)

noise_level = 0.01;
S = abs(S + randn(size(S))*noise_level + 1i * randn(size(S))*noise_level);


% fit_DT method:
[ DT, FA, ADC, VectorField, EigVals] = fit_DT( S, bval);

disp('Method 1:')
disp(['ADC = ' num2str(ADC)])
disp(['FA = ' num2str(FA)])



% fit_DT_YHd method:

S0 = mean(S(1:2));

[DT, FA, ADC, EigVects, EigVals] = fit_DT_YHd_method(permute(S(3:end), [1 3 4 2]), ...
    permute(S0, [1 3 4 2]), dirdata(3:end, :) , b(3:end) );

disp('Method 2:')
disp(['ADC = ' num2str(ADC)])
disp(['FA = ' num2str(FA)])

% In Tunnicliffe et al 2014, they compensate for the non-zero b=0 images by
% subtracting this value from the upper b values
[DT, FA, ADC, EigVects, EigVals] = fit_DT_YHd_method(permute(S(3:end), [1 3 4 2]), ...
    permute(S0, [1 3 4 2]), dirdata(3:end, :) , b(3:end) - b(1) );

disp('Method 2 - compensating for non-zero b=0 images:')
disp(['ADC = ' num2str(ADC)])
disp(['FA = ' num2str(FA)])




