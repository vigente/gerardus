function [ Model_coeffs ] = fit_stretched_exponential( im, b, thresh_val, bounded )
% FIT_STRETCHED_EXPONENTIAL    Fit the stretched exponential model of 
% diffusion to an image
%
%   Traditional tensor model: S = S0 exp( - b*D )
%   Stretched exponential model: S = S0 exp(-(b*DDC)^alpha)
%   See Bennet 2003 - Characterization of Continuously Distributed Cortical
%   Water Diffusion Rates With a Stretched-Exponential Model
%
%   IM is the input image, of any dimensionality, with the diffusion 
%   scans in the last dimension
%     
%   B is the b-matrix, of size [3 3 N]. If the measured b matrix isn't 
%   available, it can be approximated by:
%   b(:,:,n)=b_values(n)*unit_vectors(n,:)'*unit_vectors(n,:);
% 
%   THRESH_VAL is a threshold for skipping the FA, ADC and vector field 
%   computation in voxels x where im(x) < thresh_val.
%
%   BOUNDED forces all coefficients to be within certain bounds (recommended)
%
%   MODEL_COEFFS is the same size as the original image, with the last
%   dimension as the coefficients of the stretched exponential model as
%   follows:  1-6 - DDC (xx, xy, xz, yy, yz, zz)
%               7 - S0
%               8 - alpha
%
%   A note on usage: This model shouldn't really be used unless you have
%   multiple shells. Given a single shell, it can't accurately estimate 
%   the alpha parameter. 
%   While the function *can* fit a single alpha for all directions, the
%   approach in Bennet (2003) and Hall (2008) seems to be fitting
%   parameters in one direction at a time, with many samples with
%   differing |b| values in a particular direction.
%
%   See also fit_DT


% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright © 2014 University of Oxford
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
narginchk(2, 4);
nargoutchk(0, 1);

if nargin < 3
    thresh_val = -inf;
end
if nargin < 4
    bounded = true;
end

sz = size(im);

% take the log of the image to linearise the equation
imlog = log(im);

% Sort all b matrices in to a vector Bv=[Bxx,2*Bxy,2*Bxz,Byy,2*Byz,Bzz];
Bv=squeeze([b(1,1,:),2*b(1,2,:),2*b(1,3,:),b(2,2,:),2*b(2,3,:),b(3,3,:)])';


% Add another column to Bv to handle the constant term:
% Slog = Bv * M + log(S0)
% becomes:
% Slog = [Bv, -1] * [M; -log(S0)]
Bv = [Bv, -ones(size(Bv,1), 1)];

% reshape Ilog so we can use the efficient \ operator
imlog = reshape(imlog, [prod(sz(1:end-1)), sz(end)]);
imvector = reshape(im, [prod(sz(1:end-1)), sz(end)]);

% Fit the traditional tensor model (alpha = 1)
M = (Bv \ -imlog')';

% Turn log(S0) into S0
M(:,7) = exp(M(:,7));

% Set up the anomalous diffusion model
if bounded
    % we have to use lsqnonlin, because it lets you pick upper and lower
    % bounds
    options = optimset('display','off');
    lb = zeros(1,8);
    ub = [ones(1,6)*10E-3, inf, 1]; % See Bennet et al.
else
    % we can use nlinfit (slightly faster than lsqnonneg)
    modelfun = @(x, Bv) x(7) .* exp(-1 * ...
                 (x(1) * Bv(:,1) + x(2) * Bv(:,2) + x(3) * Bv(:,3) + ...
                  x(4) * Bv(:,4) + x(5) * Bv(:,5) + x(6) * Bv(:,6)).^x(8));
end


Model_coeffs = zeros(size(M,1), 8);
for i = 1:size(M,1)
    disp(['Completed ' num2str(i) ' of ' num2str(size(M,1))])
    if im(i) < thresh_val % skip if signal is low
        continue
    end
    
    % use a nonlinear fit, because I can't figure out a way to linearise
    % this equation, so it is very slow.
    x0 = [M(i, :), 1]; % initial fit has alpha = 1
    if bounded    
        modelfun = @(x) x(7) .* exp(-1 * (x(1) * Bv(:,1) + ...
                        x(2) * Bv(:,2) + x(3) * Bv(:,3) + x(4) * Bv(:,4) + ...
                        x(5) * Bv(:,5) + x(6) * Bv(:,6)).^x(8)) ...
                        - imvector(i,:)';
        Model_coeffs(i,:) = lsqnonlin(modelfun, abs(x0), lb, ub, options);
    else
        Model_coeffs(i,:) = nlinfit(Bv, imvector(i, :)', modelfun, x0);
    end
end

% reshape to original size
Model_coeffs = reshape(Model_coeffs, [sz(1:end-1), 8]);


