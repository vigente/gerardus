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
%   BVAL is the b-matrix, of size [3 3 N] (same as in fit_DT)
%
%   MASK is the mask of voxels for analysis, same size as the first N-1
%   dimensions of I. 
%
%   NREPS is the number of Monte Carlo repetitions (default 1000)
%
% Outputs:
%
%   COU is the cone of uncertainty (in degrees) (See Jones 2003 - 
%   Determining and Visualizing Uncertainty in Estimates of Fiber 
%   Orientation From Diffusion Tensor MRI). Here, the 95th percentile is
%   returned. The COU of the second and third eigenvectors is concatenated
%   in the last dimension.
%
%   FA_STD is the standard deviation of the fractional anisotropy over 
%   the repetitions
%
%   ADC_STD is the standard deviation of the ADC over the repetitions
    
% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright � 2014 University of Oxford
% Version: 0.1.2
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
narginchk(2,4);
nargoutchk(0, 3);

sz = size(I);

if nargin < 3
    Mask = zeros(sz(1:end-1));
end

if nargin < 4
    Nreps = 1000;
end


% save memory and time by fitting in vector form
I = reshape(I, [prod(sz(1:3)), sz(4)]);

Ivector = I(Mask(:), :);

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
    
    VectorField2 = real(VectorField2);
    
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
Angle_deviation_primary = zeros(size(FA_reps));
Angle_deviation_secondary = zeros(size(FA_reps));
Angle_deviation_tertiary = zeros(size(FA_reps));

for n = 1:Nreps
    % primary eigenvectors
    v1 = real(squeeze(VectorField(:,1,1:3)));
    v2 = real(squeeze(VectorField_reps(:,1,1:3,n)));
    
    % ensure unit magnitude
    v1 = bsxfun(@rdivide, v1, sqrt(sum(v1.^2, 2))+eps);
    v2 = bsxfun(@rdivide, v2, sqrt(sum(v2.^2, 2))+eps);
    
    theta = dot(v1, v2, 2);
    Angle_deviation_primary(:,1,n) = acos(theta) / pi * 180;
    
    % secondary eigenvectors
    v1 = real(squeeze(VectorField(:,1,4:6)));
    v2 = real(squeeze(VectorField_reps(:,1,4:6,n)));
    
    % ensure unit magnitude
    v1 = bsxfun(@rdivide, v1, sqrt(sum(v1.^2, 2))+eps);
    v2 = bsxfun(@rdivide, v2, sqrt(sum(v2.^2, 2))+eps);
    
    theta = dot(v1, v2, 2);
    Angle_deviation_secondary(:,1,n) = acos(theta) / pi * 180;
    
    % tertiary eigenvectors
    v1 = real(squeeze(VectorField(:,1,7:9)));
    v2 = real(squeeze(VectorField_reps(:,1,7:9,n)));
    
    % ensure unit magnitude
    v1 = bsxfun(@rdivide, v1, sqrt(sum(v1.^2, 2))+eps);
    v2 = bsxfun(@rdivide, v2, sqrt(sum(v2.^2, 2))+eps);
    
    theta = dot(v1, v2, 2);
    Angle_deviation_tertiary(:,1,n) = acos(theta) / pi * 180;
    
    
end

% bigger than 90 degrees? flip it over
Angle_deviation_primary(Angle_deviation_primary > 90) = 180 - Angle_deviation_primary(Angle_deviation_primary > 90);
Angle_deviation_secondary(Angle_deviation_secondary > 90) = 180 - Angle_deviation_secondary(Angle_deviation_secondary > 90);
Angle_deviation_tertiary(Angle_deviation_tertiary > 90) = 180 - Angle_deviation_tertiary(Angle_deviation_tertiary > 90);

% 95th percentile
COU_primary = zeros(size(Mask));
COU_primary(Mask) = prctile(Angle_deviation_primary, 95, 3);
COU_secondary = zeros(size(Mask));
COU_secondary(Mask) = prctile(Angle_deviation_secondary, 95, 3);
COU_tertiary = zeros(size(Mask));
COU_tertiary(Mask) = prctile(Angle_deviation_tertiary, 95, 3);

n = ndims(COU_primary);
COU = cat(n+1, COU_primary, COU_secondary, COU_tertiary);


