function [ DKI ] = fit_DKI_model( im, B_mat, thresh_val, method)
% FIT_DT    Fits the diffusion kurtosis model voxelwise to an image
%           S = S0 exp(-(bD + 1/6 b^2 * K))
%
%
% Inputs:
%
%   IM is the input image, of any dimensionality, with the diffusion 
%   scans in the last dimension
%     
%   B is the b-matrix, of size [3 3 N]. If the measured b matrix isn't 
%   available, it can be approximated by:
%   b(:,:,n)=b_values(n)*unit_vectors(n,:)'*unit_vectors(n,:);
% 
%   THRESH_VAL is a threshold for skipping the FA, ADC and vector field 
%   computation in voxels x where im(x) < thresh_val. THRESH_VAL can
%   alternatively be a mask of voxels to process.
%
%   METHOD is a string, and can take the following values:
%       'linear' fits the tensor in the log domain (default)
%       'nonlinear' fits the tensor in the signal domain (recommended for
%       high b-values or when SNR is low)
%
% Outputs:
%
%   DKI is the array of coefficients
% 
%   FA is the fractional anisotropy
% 
%   ADC is the apparent diffusion coefficient
% 
%   VECTORFIELD is the primary, secondary and tertiary unit vectors, concatenated
%	in the last dimension
%
%	EIGVALS are the eigenvalues corresponding to the vector field, concatenated
%	in the last dimension
%
%

    
% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright © 2014-2015 University of Oxford
% Version: 0.1.0
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
%narginchk(2, 5);
%nargoutchk(0, 5);

sz = size(im);

if nargin < 3
    thresh_val = -inf;
end

if (nargin < 4) || isempty(method)% default method is linear least squares
    method = 'linear';
end


% check that a valid method has been given
if ~max([strcmp(method, 'linear'), strcmp(method, 'nonlinear')])
    disp('Unrecognised method, performing fast linear fitting')
end



% handle a vector with the wrong orientation
if (length(sz) == 2) && (sz(2) == 1)
    im = im';
    sz = size(im);
end

% reshape image
I_vector = reshape(im, [prod(sz(1:end-1)), sz(end)]);


% normalise and take the log of the image to linearise the equation
imlog = log(bsxfun(@rdivide, I_vector, I_vector(:,1)));


% get the bvalue
bvalue = squeeze(B_mat(1,1,:) + B_mat(2,2,:) + B_mat(3,3,:))';

B_mat = B_mat ./ repmat(permute(bvalue, [3,1,2]), [3 3 1]);



% define A matrix

A_D = squeeze([B_mat(1,1,:), B_mat(2,2,:), B_mat(3,3,:), 2*B_mat(1,2,:), 2*B_mat(1,3,:), 2*B_mat(2,3,:)])';
% A_D in the paper is not scaled by the b value

A_K = squeeze([B_mat(1,1,:).^2, B_mat(2,2,:).^2, B_mat(3,3,:).^2, ...
    4*B_mat(1,1,:).*B_mat(1,2,:), 4*B_mat(1,1,:).*B_mat(1,3,:), 4*B_mat(1,2,:).*B_mat(2,2,:), ...
    4*B_mat(1,3,:).*B_mat(3,3,:), 4*B_mat(2,2,:).*B_mat(2,3,:), 4*B_mat(2,3,:).*B_mat(3,3,:), ...
    2*B_mat(1,1,:).*B_mat(2,2,:) + 4*B_mat(1,2,:).*B_mat(1,2,:), ...
    2*B_mat(1,1,:).*B_mat(3,3,:) + 4*B_mat(1,3,:).*B_mat(1,3,:), ...
    2*B_mat(2,2,:).*B_mat(3,3,:) + 4*B_mat(2,3,:).*B_mat(2,3,:), ...
    8*B_mat(1,2,:).*B_mat(1,3,:) + 4*B_mat(1,1,:).*B_mat(2,3,:), ...
    8*B_mat(1,2,:).*B_mat(2,3,:) + 4*B_mat(2,2,:).*B_mat(1,3,:), ...
    8*B_mat(1,3,:).*B_mat(2,3,:) + 4*B_mat(3,3,:).*B_mat(1,2,:)])';
% this takes care of cross terms by saying that (for example) 6 x y^2
% = 2 (x^2) (y^2) + 4 (xy) (xy)
% so it uses these elements of the B matrix.


A = [-bsxfun(@times, bvalue', A_D), 1/6 * bsxfun(@times, bvalue'.^2, A_K)];

% weighted fit
W = diag(mean(I_vector));

%M = (pinv(A) * imlog')';
M = (pinv(A' * (W.^2) * A) * A' * (W.^2) * (imlog'))';
% weighted version of M = (pinv(A) * imlog')';

% add the S0 column
M = [I_vector(:,1), M];


% perform nonlinear fitting (if required)
if strcmp(method, 'nonlinear')

    if isscalar(thresh_val)
        thresh_val = I_vector(:,1) > thresh_val;
    else
        thresh_val = logical(thresh_val);
    end
    
    % declare variables for nonlinear fitting
    I_nlfit = double(I_vector(thresh_val(:), :));
    M_nl = double(M(thresh_val(:), :));

    % get error
    I_recon = bsxfun(@times, M_nl(:,1) , exp(A * M_nl(:,2:end)')');
    rmse = sqrt(mean((I_nlfit(:) - I_recon(:)).^2));
    disp(['RMSE from quick fit = ' num2str(rmse)])
    
    options = optimoptions('lsqcurvefit','Jacobian','on', 'DerivativeCheck', 'off', ...
        'display', 'off', 'MaxIter', 40);

    
    lb = [0 0 0 0 -1 -1 -1 ... % [S0 DTI
            0 0 0 -1 -1 -1 -1 -1 -1 0 0 0 -1 -1 -1] * 3E-3; % DKI ]
        
    ub = [inf, ones(1,21) * 3E-3]; % [S0 DTI DKI]
    
    % constrain initial guess between lower and upper bounds
    M_nl = bsxfun(@max, M_nl, lb);
    M_nl = bsxfun(@min, M_nl, ub);
    
    
    if (size(I_nlfit,1) > 1000) && (matlabpool('size') == 0)
        disp('Open matlabpool for parallel processing')
    end
    
    % change A for nonlinear fitting
    A = -A';
    
    parfor i = 1:size(I_nlfit,1)
    
        M_nl(i,:) = lsqcurvefit(@DKI_model, M_nl(i,:), A, I_nlfit(i,:), lb, ub, options);
        
    end
    
    I_recon = bsxfun(@times, M_nl(:,1) , exp(-A' * M_nl(:,2:end)')');
    rmse = sqrt(mean((I_nlfit(:) - I_recon(:)).^2));
    disp(['RMSE from slow fit = ' num2str(rmse)])


    M(thresh_val(:),:) = M_nl;
    
end
    
    
% reshape back again
DKI = reshape(M, [sz(1:end-1), 22]);

end


function [ F, J ] = DKI_model( x, Bv )

F = x(1) * exp(-(x(2:end) * Bv));
    

if nargout > 1 % Jacobian
    J = zeros(length(F), length(x));
    
    J(:,1) = exp(-(x(2:end) * Bv));

    J(:,2:end) = -bsxfun(@times, Bv, F)';
    
end
    


end

