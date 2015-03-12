function [AA, FD, S0_all, gamma_all, A_all, Npoints] = fit_stretched_exponential(I0, Ihigh, bval0, bvalhigh, mask, noise_thresh)
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
%       S0_ALL are the coefficients of S0 [r c s dir]
%       GAMMA_ALL are the coefficients of gamma [r c s dir]
%       A_ALL are the coefficients of A [r c s dir]
%       NPOINTS is the number of fit shells at each voxel (see Hall 2008 
%           Figure 3e)       
%
%   See also fit_DT, weighted_linear_fit


% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.1.4
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
nargoutchk(0, 6);

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
S0_all = zeros(sz(1:4));
%mask = (I0_mean > 2E4) & (I0_mean ./ Ihigh(:,:,:,1,1) < 3);


% as in Hall 2008, we only want to fit the linear regression where we
% are above the noise floor
Imin = min(Ihigh, [], 4);
Npoints = sum(Imin > noise_thresh, 5)+1; % (includes b = 0)

nvoxels = prod(sz(1:3));

for dir = 1:sz(4) % for each direction
    disp(['Dir = ' num2str(dir)])
    
    bval_all = cat(3, bval0, squeeze(bvalhigh(:,:,dir,:)));
    
    b = zeros(1, size(bval_all,3));
    for nshell = 1:size(bval_all,3)
        b(nshell) = bval_all(1,1,nshell) + bval_all(2,2,nshell) + bval_all(3,3,nshell);
    end

    %disp(b)
    
    I_all = cat(4, I0_mean, permute(Ihigh(:,:,:,dir,:), [1 2 3 5 4]));
    
    Ivector = reshape(I_all, [prod(sz(1:3)), size(I_all,4)]);
    
    S0 = Ivector(:,1);
    A = zeros(nvoxels,1);
    gamma = zeros(nvoxels,1);
    
    
    for n = 3:max(Npoints(:)) % n = number of points for linear fitting
        %disp(['Num shells = ' num2str(n)])
        
        idx = Npoints(:) == n;
        S0_idx = S0(idx);
        
        % eqn 1 - use S0 and solve for gamma, A
        
        
        I_idx = Ivector(idx,1:n);
        
        % take the log for linear fitting
        loglogsig = log(abs(log(eps + bsxfun(@rdivide, I_idx, S0_idx))));
        logb = log(b(1:n));

        b_array = [logb(:), ones(size(logb(:)))];


        %coeff = (b_array(2:end,:) \ loglogsig(:,2:end)')';
        func = @(x) exp(-abs(exp(x)));
        coeff = weighted_linear_fit(loglogsig(:,2:end), b_array(2:end,:)', func,0);

        gamma_idx = coeff(:,1);
        A_idx = exp(coeff(:,2));

%             sig_recon = S0_idx(1) .* exp(-A_idx(1) .* b(1:n).^gamma_idx(1));
%             err = I_idx(1,:) - sig_recon;
% 
%             %rmse = S0_idx .* exp(A_idx .* b(1:n).^gamma_idx) - Ivector
%             disp(['gamma = ' num2str(gamma_idx(1)), ...
%                 ', A = ' num2str(A_idx(1)) ', S0 = ' num2str(S0_idx(1)), ', rmse = ' num2str(sqrt(mean(err.^2)))])
% 
%             % eqn 2 - use gamma and solve for S0, A
%             gamma_idx = constrain(gamma_idx, 0, 1.5);
%             logsignal = log(Ivector(idx,1:n));
% 
%             % the equation can be linearised in another way, as 
%             % log(S) = log(S0) - A*b^gamma
%
%             for vox = 1:size(logsignal,1)
%                 b_pow = b(1:n).^gamma_idx(vox);
%                 b_array = [b_pow(:), ones(size(b_pow(:)))];
% 
%                 % only using n-1 here because I had issues from the log
%                 % operation
%                 %coeff = b_array(1:end-1,:) \ logsignal(vox,1:end-1)';
%                 coeff = weighted_linear_fit(logsignal(vox,1:end), b_array(1:end,:)', @(x)exp(x), 0); 
%                 
%                 A_idx(vox) = -coeff(1);
%                 S0_idx(vox) = exp(coeff(2));
% 
%             end
%             
%             sig_recon = S0_idx(1) .* exp(-A_idx(1) .* b(1:n).^gamma_idx(1));
%             err = I_idx(1,:) - sig_recon;
% 
%             %rmse = S0_idx .* exp(A_idx .* b(1:n).^gamma_idx) - Ivector
%             disp(['gamma = ' num2str(gamma_idx(1)), ...
%                 ', A = ' num2str(A_idx(1)) ', S0 = ' num2str(S0_idx(1)), ', rmse = ' num2str(sqrt(mean(err.^2)))])

        
        gamma(idx) = gamma_idx;
        A(idx) = A_idx;
        S0(idx) = S0_idx;
        
    end
    
    if max(mask(:)) == 1 % if the mask is not empty
    
        
        options = optimoptions('lsqcurvefit','Jacobian','on', 'DerivativeCheck', 'off', 'display', 'off');
        lb = zeros(1,3);
        ub = [1, 10E-3, inf]; % See Bennet et al. for bounds

        x_guess = abs(double([gamma(mask(:)), A(mask(:)), Ivector(mask(:),1)]));
        new_coeffs = zeros(size(x_guess));
        
        I_nlfit = double(Ivector(mask(:), :));
        
        
        if (matlabpool('size') == 0) && (size(new_coeffs,1) > 100) % if there is a decent number of voxels
            disp('Open matlab pool if you want to use parallel processing')
        end
        
        parfor i = 1:size(new_coeffs,1)

            % alternative (slower) method:
            %modelfun = @(x) x(3) .* exp(- x(2) * (b.^x(1))) - I_nlfit(i,:);
            %new_coeffs(i,:) = lsqnonlin(modelfun, x_guess(i,:), lb, ub, options);
            
            new_coeffs(i,:) = lsqcurvefit(@stretched_exp_model, x_guess(i,:), b, I_nlfit(i,:), lb, ub, options);
            
        end

        gamma(mask(:)) = new_coeffs(:,1);
        A(mask(:)) = new_coeffs(:,2);
        S0(mask(:)) = new_coeffs(:,3);

    end
    

    gamma_all(:,:,:,dir) = reshape(gamma, sz(1:3));
    A_all(:,:,:,dir) =  reshape(A, sz(1:3));
    S0_all(:,:,:,dir) = reshape(S0, sz(1:3));
end

% Anomalous anisotropy
mean_gamma = mean(gamma_all, 4);
N = sz(4);
AA = sqrt(N / (N - 1) * sum((gamma_all - repmat(mean_gamma, [1 1 1 N])).^2, 4)...
    ./ sum(gamma_all.^2, 4));

% Fractal dimension
FD = 2 / mean_gamma;



% the function to fit inside the parfor loop
function [F, J] = stretched_exp_model(x, b)

F = x(3) .* exp(- x(2) * (b.^x(1)));

if nargout > 1 % Jacobian - makes fitting about 20% faster
    J = zeros(length(b), 3);
    J(:,1) = -x(2) * x(3) * b.^(x(1)) .* log(b) .* exp(- x(2) * (b.^x(1)));
    J(:,2) = - b.^(x(1)) * x(3) .* exp(- x(2) * (b.^x(1)));
    J(:,3) = exp(- x(2) * (b.^x(1)));
end
