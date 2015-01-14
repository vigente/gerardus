function [COU, FA_STD, ADC_STD] = wild_bootstrapping_DTI( I, bval, Mask, Nreps )
% WILD_BOOTSTRAPPING_DTI Returns uncertainty measures from diffusion tensor
% images via wild bootstrapping. See Whitcher 2008 - Using the Wild 
% Bootstrap to Quantify Uncertainty in Diffusion Tensor Imaging
%
%
% Inputs:
%
%   I is the input image, of any dimensionality, with the diffusion 
%   scans in the last dimension
%     
%   BVAL is the b-matrix, of size [3 3 N]
%
%   MASK is the mask of voxels for analysis, same size as the first N-1
%   dimensions of I
%
%   NREPS is the number of Monte Carlo repetitions (default 1000)
%
% Outputs:
%
%   COU is the cone of uncertainty (in degrees) (See Jones 2003 - 
%   Determining and Visualizing Uncertainty in Estimates of Fiber 
%   Orientation From Diffusion Tensor MRI). Here, the 95th percentile is
%   returned.
%
%   FA_STD is the standard deviation of the fractional anisotropy over 
%   the repetitions
%
%   ADC_STD is the standard deviation of the ADC over the repetitions
    
% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright � 2014 University of Oxford
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
narginchk(3,4);
nargoutchk(0, 3);

if nargin == 3
    Nreps = 1000;
end


sz = size(I);

% save memory and time by fitting in vector form
I = reshape(I, [prod(sz(1:3)), sz(4)]);

Ivector = I(Mask, :);

% fit the tensor
[DT, FA, ADC, VectorField, EigVals] = fit_DT(Ivector, bval);


% Model fitted image
Ifit = dt2image(DT, bval);
% Residuals
Resids = Ivector - Ifit;

DT_reps = zeros([size(DT), Nreps]);
FA_reps = zeros([size(FA), Nreps]);
ADC_reps = zeros([size(ADC), Nreps]);
VectorField_reps = zeros([size(VectorField), Nreps]);
EigVals_reps = zeros([size(EigVals), Nreps]);



% for N iterations, randomly multiply the residuals by 1 or -1 and compute
% parameters
for n = 1:Nreps
    disp(n)
    
    % either 1 or -1, with 50% probability each 
    F = rand(size(Resids)) > 0.5;
    F = F * 2 - 1;

    Resids_to_add = Resids .* F;

    Inew = Ifit + Resids_to_add;

    [DT2, FA2, ADC2, VectorField2, EigVals2] = fit_DT(Inew, bval);

    DT_reps(:,:,n) = DT2;
    FA_reps(:,1,n) = FA2;
    ADC_reps(:,1,n) = ADC2;
    VectorField_reps(:,:,:,n) = VectorField2;
    EigVals_reps(:,:,:,n) = EigVals2;
end

FA_STD = zeros(size(Mask));
FA_STD(Mask) = std(FA_reps, [], 3);

ADC_STD = zeros(size(Mask));
ADC_STD(Mask) = std(ADC_reps, [], 3);

% angle between original data and bootstrapped data
Angle_deviation = zeros(size(FA_reps));
for n = 1:Nreps
    A = real(squeeze(VectorField(:,1,1:3)));

    B = real(squeeze(VectorField_reps(:,1,1:3,n)));
    
    % ensure unit magnitude
    
    A = bsxfun(@rdivide, A, sqrt(sum(A.^2, 2))+eps);
    B = bsxfun(@rdivide, B, sqrt(sum(B.^2, 2))+eps);
    
    C = dot(A, B, 2);
    
    Angle_deviation(:,1,n) = acos(C) / pi * 180;
end

% bigger than 90 degrees? flip it over
Angle_deviation(Angle_deviation > 90) = 180 - Angle_deviation(Angle_deviation > 90);

% 95th percentile
COU = zeros(size(Mask));
COU(Mask) = prctile(Angle_deviation, 95, 3);




