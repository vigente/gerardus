function [AA, FD, gamma_all, A_all, Npoints] = fit_stretched_exponential(I0, Ihigh, bval0, bvalhigh, mask, noise_thresh)
% FIT_STRETCHED_EXPONENTIAL    Fit the stretched exponential model of 
% diffusion to an image
%
%   Traditional tensor model: S = S0 exp( - b*D )
%   Stretched exponential model has two definitions:
%           Bennet: S = S0 exp(-(b*DDC)^alpha)
%           Hall:   S = S0 exp(- b^gamma * A)
%                   where gamma = alpha, A = DDC^gamma
%   Bennet 2003 - Characterization of Continuously Distributed Cortical
%   Water Diffusion Rates With a Stretched-Exponential Model
%   and 
%   Hall 2008 - From Diffusion-Weighted MRI to Anomalous Diffusion Imaging
%
%   Here, we use the Hall 2008 definition because it is easily linearisable
%
%   Inputs:
%       I0 is the un-attenuated image [row col slice repetition]
%           If there was only a single I0 acquired, just pass a 3D array
%       IHIGH is the attenuated images [row col slice direction shell]
%           E.g. if you have 21 different directions and 3 shells, the data
%           size would be [r c s 21 3]
%           The function assumes the shells are in ascending order of
%           b-value.
%       BVAL0 is the b-matrix of I0, of size [3 3 rep]
%       BVALHIGH is the b-matrix of Ihigh, of size [3 3 direction shell]
%       MASK is the mask of voxels that you want to be fit with
%           a nonlinear regression. All voxels get fit with a linear fit 
%           (very fast), whereas the voxels inside the mask are repeated 
%           using a nonlinear regression with bounds on parameters (very 
%           slow). I have found that the nonlinear fit changes the AA and 
%           FD by 2-3% in 3-shell data. (optional)
%       NOISE_THRESH is the noise threshold for fitting the model. Only
%           points in IHIGH that are above this value will be fit - the
%           others will be ignored. For example, in fluid only the 
%           innermost shells will be above the noise floor, so we only fit 
%           them. But in tissue, more shells will be above the noise floor. 
%           See the output variable NPOINTS. (optional, but recommended)
%     
%   Outputs:
%       AA is Anomalous Anisotropy [r c s]
%       FD is Fractal Dimension [r c s]
%       GAMMA_ALL are the coefficients of gamma [r c s dir]
%       A_ALL are the coefficients of A [r c s dir]
%       NPOINTS is the number of fit shells at each voxel (see Hall 2008 
%           Figure 3e)       
%
%   See also fit_DT


% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.1.3
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
narginchk(4, 6);
nargoutchk(0, 5);

if (nargin < 5) || isempty(mask)
    mask = 0;
else
    mask = mask == 1; % convert to logical
end

sz = size(Ihigh);

if (nargin < 6) || isempty(noise_thresh)
   % assume there is just air in the top left corner
    noise_thresh = mean(mean(mean(Ihigh(1:min(sz(1),10), 1:min(sz(2),10), 1:min(sz(3),10), 1,1))));
end

I0_mean = mean(I0, 4);


gamma_all = zeros(sz(1:4));
A_all = zeros(sz(1:4));

%mask = (I0_mean > 2E4) & (I0_mean ./ Ihigh(:,:,:,1,1) < 3);


% as in Hall 2008, we only want to fit the linear regression where we
% are above the noise floor
Imin = min(Ihigh, [], 4);
Npoints = sum(Imin > noise_thresh, 5);

for dir = 1:sz(4) % for each direction

    b0 = bval0(:,:,1);
    b1 = bvalhigh(:,:,dir,1);
    b2 = bvalhigh(:,:,dir,2);
    b3 = bvalhigh(:,:,dir,3);

    % magnitude of b
    b = [sqrt(sum(b1(:).^2)), sqrt(sum(b2(:).^2)), sqrt(sum(b3(:).^2))]';

    % normalise by I0
    I_dir = permute(Ihigh(:,:,:,dir,:), [1 2 3 5 4]) ./ repmat(I0_mean, [1 1 1 size(Ihigh,5)]);
    
    Ivector = reshape(I_dir, [prod(sz(1:3)), sz(5)]);
    
    coeff = zeros(prod(sz(1:3)), 2);
    
    for n = 2:max(Npoints(:)) % n = number of points for linear fitting
        idx = Npoints(:) == n;
        
        % take the log for linear fitting
        loglogsig = log(abs(log(Ivector(idx, 1:n))));
        logb = log(b(1:n));

        b_array = [logb(:), ones(size(logb(:)))];

        % alternatively, fit the linear model using a weighted fit
    %     w = zeros(3);
    %     w(1,1) = 1;
    %     w(2,2) = 0.5;
    %     w(3,3) = 0.1;
    %     coeff(idx,:) = (pinv(b_array' * (w.^2) * b_array) * b_array' * (w.^2) * (loglogsig'))';
        coeff(idx, :) = (b_array \ loglogsig')';
    end
    
    if max(mask(:)) == 1 % if the mask is not empty
    
        % use x as a best guess for the non-linear fitting
        I_dir = squeeze(Ihigh(:,:,:,dir,:)); % not normalised this time
        b = [sqrt(sum(b0(:).^2)), sqrt(sum(b1(:).^2)), sqrt(sum(b2(:).^2)), sqrt(sum(b3(:).^2))];

        % do the nonlinear fitting only on the voxels in the mask
        Ivector = reshape(I_dir, [prod(sz(1:3)), sz(5)]);
        Ivector = double([I0_mean(mask), Ivector(mask(:), :)]);

        options = optimset('display','off');
        lb = zeros(1,3);
        ub = [1, 10E-3, inf]; % See Bennet et al. for bounds

        new_coeffs = zeros(size(Ivector,1), 3);
        old_coeffs = coeff(mask(:),:);

        s = matlabpool('size');
        if (s == 0) && (size(new_coeffs,1) > 100) % if there is a decent number of voxels
            disp('Open matlab pool to use parallel processing')
        end
        
        parfor i = 1:size(new_coeffs,1)
            
            x_guess = double([old_coeffs(i,1), exp(old_coeffs(i,2)), Ivector(i,1)]);

            modelfun = @(x) x(3) .* exp(- x(2) * (b.^x(1))) - Ivector(i,:);

            new_coeffs(i,:) = lsqnonlin(modelfun, abs(x_guess), lb, ub, options);
        end

        coeff(mask(:), :) = new_coeffs(:,1:2);
    end
    
    
    gamma = reshape(coeff(:,1), sz(1:3));
    
    A = reshape(exp(coeff(:,2)), sz(1:3));

    gamma_all(:,:,:,dir) = gamma;
    A_all(:,:,:,dir) = A;
end

% Anomalous anisotropy
mean_gamma = mean(gamma_all, 4);
N = sz(4);
AA = sqrt(N / (N - 1) * sum((gamma_all - repmat(mean_gamma, [1 1 1 N])).^2, 4)...
    ./ sum(gamma_all.^2, 4));

% Fractal dimension
FD = 2 / mean_gamma;

