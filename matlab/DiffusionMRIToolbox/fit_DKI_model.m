function [ DKI, K_ADC, K_AKC, Kv1, Kv2, Kv3, I_DKI ] = fit_DKI_model( im, B_mat, thresh_val, method)
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


    
% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright © 2014-2015 University of Oxford
% Version: 0.1.1
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

if (nargin < 3) || isempty(thresh_val)
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

% handle mask or threshold
if isscalar(thresh_val)
    thresh_val = I_vector(:,1) > thresh_val;
else
    thresh_val = logical(thresh_val);
end

% keep only the masked voxels
I_vector = double(I_vector(thresh_val(:), :));

[~, b0] = max(mean(I_vector, 1));

% normalise and take the log of the image to linearise the equation
imlog = log(abs(bsxfun(@rdivide, I_vector, I_vector(:,b0))));


% get the bvalue
bvalue = squeeze(B_mat(1,1,:) + B_mat(2,2,:) + B_mat(3,3,:))';

% normalise the B matrix
B_mat_norm = B_mat ./ repmat(permute(bvalue, [3,1,2]), [3 3 1]);


% define A matrix

A_D = squeeze([B_mat_norm(1,1,:), B_mat_norm(2,2,:), B_mat_norm(3,3,:), 2*B_mat_norm(1,2,:), 2*B_mat_norm(1,3,:), 2*B_mat_norm(2,3,:)])';
% A_D in the paper is not scaled by the b value


A_K = squeeze([B_mat_norm(1,1,:).^2, B_mat_norm(2,2,:).^2, B_mat_norm(3,3,:).^2, ...
    4*B_mat_norm(1,1,:).*B_mat_norm(1,2,:), ...
    4*B_mat_norm(1,1,:).*B_mat_norm(1,3,:), ...
    4*B_mat_norm(1,2,:).*B_mat_norm(2,2,:), ...
    4*B_mat_norm(1,3,:).*B_mat_norm(3,3,:), ...
    4*B_mat_norm(2,2,:).*B_mat_norm(2,3,:), ...
    4*B_mat_norm(2,3,:).*B_mat_norm(3,3,:), ...
    2*B_mat_norm(1,1,:).*B_mat_norm(2,2,:) + 4*B_mat_norm(1,2,:).*B_mat_norm(1,2,:), ...
    2*B_mat_norm(1,1,:).*B_mat_norm(3,3,:) + 4*B_mat_norm(1,3,:).*B_mat_norm(1,3,:), ...
    2*B_mat_norm(2,2,:).*B_mat_norm(3,3,:) + 4*B_mat_norm(2,3,:).*B_mat_norm(2,3,:), ...
    8*B_mat_norm(1,2,:).*B_mat_norm(1,3,:) + 4*B_mat_norm(1,1,:).*B_mat_norm(2,3,:), ...
    8*B_mat_norm(1,2,:).*B_mat_norm(2,3,:) + 4*B_mat_norm(2,2,:).*B_mat_norm(1,3,:), ...
    8*B_mat_norm(1,3,:).*B_mat_norm(2,3,:) + 4*B_mat_norm(3,3,:).*B_mat_norm(1,2,:)])';
% this takes care of cross terms by saying that (for example) 6 x^2 y^2
% = 2 (x^2) (y^2) + 4 (xy) (xy)
% so it uses these elements of the B matrix.



A = [-bsxfun(@times, bvalue', A_D), 1/6 * bsxfun(@times, bvalue'.^2, A_K)];

% only keep b values less than 3000 in the linear fit
lin_bvals = bvalue < 5000;

% weighted fit
param.rician = 0;
param.unique_weights = 0;
param.verbose = 1;
M = weighted_linear_fit(imlog(:,lin_bvals), A(lin_bvals,:)', @(z)exp(z), param);
%M = weighted_linear_fit(imlog, A', @(z)exp(z), param);

%M = (pinv(A) * imlog')';
%M = (pinv(A' * (W.^2) * A) * A' * (W.^2) * (imlog'))';
% weighted version of M = (pinv(A) * imlog')';

% add the S0 column
M = [I_vector(:,1), M];

M_lin = M;

% perform nonlinear fitting (if required)
if strcmp(method, 'nonlinear')

    % get error
    I_recon = bsxfun(@times, M(:,1) , exp(M(:,2:end) * A'));
    rmse_lin = sqrt(mean((I_vector(:) - I_recon(:)).^2));
    disp(['RMSE from quick fit = ' num2str(rmse_lin)])
    
    options = optimoptions('lsqcurvefit','Jacobian','on', 'DerivativeCheck', 'off', ...
        'display', 'off', 'MaxIter', 40);

    
    lb = [0 0 0 0 -1 -1 -1 ... % [S0 DTI
            0 0 0 -1 -1 -1 -1 -1 -1 0 0 0 -1 -1 -1] * 3E-3; % DKI ]
        
    ub = [inf, ones(1,21) * 3E-3]; % [S0 DTI DKI]
    
    % constrain initial guess between lower and upper bounds
    M = bsxfun(@max, M, lb);
    M = bsxfun(@min, M, ub);
    M(isnan(M(:))) = 0;
    M(isinf(M(:))) = 0;
    
    if (size(I_vector,1) > 1000) && (matlabpool('size') == 0)
        disp('Open matlabpool for parallel processing')
    end
    
    A_hat = -A';
    
    parfor i = 1:size(I_vector,1)
    
        M(i,:) = lsqcurvefit(@DKI_model, M(i,:), A_hat, I_vector(i,:), lb, ub, options);
        
    end
    
    
    I_recon = bsxfun(@times, M(:,1) , exp(M(:,2:end) * A'));
    rmse = sqrt(mean((I_vector(:) - I_recon(:)).^2));
    
    bad_fits = rmse > rmse_lin;
    disp(['Number of bad fits: ' num2str(sum(bad_fits(:)))])
    M(bad_fits,:) = M_lin(bad_fits,:);
    
    disp(['RMSE from slow fit = ' num2str(rmse)])
    
end

% compute the returned signal
I_DKI = bsxfun(@times, M(:,1), exp(M(:, 2:end) * A'));


K_v1 = zeros(1, size(M,1));
K_v2 = zeros(1, size(M,1));
K_v3 = zeros(1, size(M,1));
ADC = zeros(1, size(M,1));
AKC = zeros(1, size(M,1));

[x, y, z] = psphere(200);
x = x(1:end/2); y = y(1:end/2); z = z(1:end/2);


% Compute the kurtosis (mean and directional)
for v = 1:size(M,1)
    
    % reshape the coefficients
    
    % S0 = M(i,1);
    
    D_hat_i = M(v,2:7);
    
    % A_D = xx, yy, zz, xy, xz, yz

    % rearrange for convenience
    D_hat = D_hat_i([1 4 5 2 6 3]);
    
    D33 = [D_hat(1), D_hat(2), D_hat(3); 
           D_hat(2), D_hat(4), D_hat(5); 
           D_hat(3), D_hat(5), D_hat(6)]; 


    [vec, val] = eig(D33);
    % sort
    [~, ix] = sort(sum(val,1), 'descend');
    vec = vec(:, ix);
    val = sum(val(:, ix));
    
    % eigenvalues should be positive
    % taking the magnitude is a very bad way of handling them
    val = abs(val);
    
    ADC(v) = sum(val) / 3;
    
    D2i = M(v, 8:end);

    % A_K = xxxx, yyyy, zzzz, 1 2 3 
    %       xxxy, xxxz, xyyy, 4 5 6
    %       xzzz, yyyz, yzzz, 7 8 9
    %       xxyy, xxzz, yyzz, 10 11 12
    %       xxyz, xyyz, xyzz ,13 14 15

    
    % re-arrange for convenience
    D2 = D2i([1, 4, 5, 10, 13, 11, 6, 14, 15, 7, 2, 8, 12, 9, 3]);


    % kurtosis tensor (includes ADC^2 internally)
    K = [ D2(01), D2(02), D2(03), D2(04), D2(05), D2(06); % xx
          D2(02), D2(04), D2(05), D2(07), D2(08), D2(09); % xy
          D2(03), D2(05), D2(06), D2(08), D2(09), D2(10); % xz
          D2(04), D2(07), D2(08), D2(11), D2(12), D2(13); % yy
          D2(05), D2(08), D2(09), D2(12), D2(13), D2(14); % yz
          D2(06), D2(09), D2(10), D2(13), D2(14), D2(15)];% zz
        %   xx      xy      xz      yy      yz      zz                      
    
    % Parallel kurtosis (kurtosis in the direction of highest diffusivity)
    v1 = vec(:, 1);
    
    % Diffusion in this direction
    D_v1 = v1' * D33 * v1; % = eigenvalue 1
    
    % kurtosis in this direction
    v1_sq = v1 * v1'; 
    v1_sq_vect = [v1_sq(1), 2*v1_sq(2), 2*v1_sq(3), v1_sq(5), 2*v1_sq(6), v1_sq(9)];
    K_v1(v) = (v1_sq_vect * K * v1_sq_vect') / (D_v1^2);

    
    % Parallel kurtosis (v2)
    v2 = vec(:, 2);
    
    % Diffusion in this direction
    D_v2 = v2' * D33 * v2; % = eigenvalue 2
    
    % kurtosis in this direction
    v2_sq = v2 * v2'; 
    v2_sq_vect = [v2_sq(1), 2*v2_sq(2), 2*v2_sq(3), v2_sq(5), 2*v2_sq(6), v2_sq(9)];
    K_v2(v) = (v2_sq_vect * K * v2_sq_vect') / (D_v2^2);
    
    
    % Parallel kurtosis (v3)
    v3 = vec(:, 3);
    
    % Diffusion in this direction
    D_v3 = v3' * D33 * v3; % = eigenvalue 3
    
    % kurtosis in this direction
    v3_sq = v3 * v3'; 
    v3_sq_vect = [v3_sq(1), 2*v3_sq(2), 2*v3_sq(3), v3_sq(5), 2*v3_sq(6), v3_sq(9)];
    K_v3(v) = (v3_sq_vect * K * v3_sq_vect') / (D_v3^2);
    
    % mean kurtosis
    % there are clever ways of getting this, but the easiest is to just
    % average the kurtosis over a bunch of uniformly spaced directions
    K_vi = zeros(1,length(x));
    for i = 1:length(x)
        vi = [x(i), y(i), z(i)]';
        D_vi = vi' * D33 * vi;
        vi_sq = vi * vi'; 
        vi_sq_vect = [vi_sq(1), 2*vi_sq(2), 2*vi_sq(3), vi_sq(5), 2*vi_sq(6), vi_sq(9)];
        K_vi(i) = (vi_sq_vect * K * vi_sq_vect') / (D_vi^2);
    end
    
    AKC(v) = mean(K_vi);

end

% fill regions outside mask with zeros, and reshape
DKI = zeros(prod(sz(1:end-1)), 22);
DKI(thresh_val(:), :) = M;
DKI = reshape(DKI, [sz(1:end-1), 22]);

IDKI = zeros(prod(sz(1:end-1)), sz(end));
IDKI(thresh_val(:), :) = I_DKI;
I_DKI = reshape(IDKI, sz);

K_ADC = zeros(size(thresh_val));
K_ADC(thresh_val) = ADC;

K_AKC = zeros(size(thresh_val));
K_AKC(thresh_val) = AKC;

Kv1 = zeros(size(thresh_val));
Kv1(thresh_val) = K_v1;

Kv2 = zeros(size(thresh_val));
Kv2(thresh_val) = K_v2;

Kv3 = zeros(size(thresh_val));
Kv3(thresh_val) = K_v3;



end


function [ F, J ] = DKI_model( x, Bv )

    F = x(1) * exp(-(x(2:end) * Bv));

    if nargout > 1 % Jacobian
        J = zeros(length(F), length(x));

        J(:,1) = exp(-(x(2:end) * Bv));

        J(:,2:end) = -bsxfun(@times, Bv, F)';

    end

end


%% This converts the 3x3 kurtosis tensor into a 3x3x3x3 kurtosis tensor
% It is not needed right now, but good to have for the future.

%     K_index =[01, 02, 03, 04, 05, 06; % xx
%               02, 04, 05, 07, 08, 09; % xy
%               03, 05, 06, 08, 09, 10; % xz
%               04, 07, 08, 11, 12, 13; % yy
%               05, 08, 09, 12, 13, 14; % yz
%               06, 09, 10, 13, 14, 15];% zz
%             % xx  xy  xz  yy  yz  zz              
%                                 
%     W = zeros([3 3 3 3]);
%     
%     for i = 1:3
%         for j = 1:3
%              for k = 1:3
%                  for l = 1:3
%                     switch i
%                         case 1
%                             r1 = [1 2 3];
%                         case 2
%                             r1 = [2 4 5];
%                         case 3
%                             r1 = [3 5 6];
%                     end
%                     switch j
%                         case 1
%                             r2 = [1 2 3];
%                         case 2
%                             r2 = [2 4 5];
%                         case 3
%                             r2 = [3 5 6];
%                     end
%                     
%                     if i == j
%                         switch i
%                             case 1
%                                 r =1;
%                             case 2
%                                 r = 4;
%                             case 3
%                                 r = 6;
%                         end
%                     else
%                         r = intersect(r1, r2);
%                     end
% 
% 
%                     switch k
%                         case 1
%                             c1 = [1 2 3];
%                         case 2
%                             c1 = [2 4 5];
%                         case 3
%                             c1 = [3 5 6];
%                     end
%                     switch l
%                         case 1
%                             c2 = [1 2 3];
%                         case 2
%                             c2 = [2 4 5];
%                         case 3
%                             c2 = [3 5 6];
%                     end
%                     
%                     if k == l
%                         switch k
%                             case 1
%                                 c =1;
%                             case 2
%                                 c = 4;
%                             case 3
%                                 c = 6;
%                         end
%                     else
%                         c = intersect(c1, c2);
%                     end
%                     
%                     D2_index = K_index(r,c);
%                     
%                     W(i,j,k,l) = D2(D2_index);
%                  end
%              end
%         end
%     end
%     v1 = vec(:, 1);
%     % Diffusion in this direction
%     D_v1 = v1' * D33 * v1;
%     K_v1_vox = 0;
%     for i = 1:3
%         for j = 1:3
%             for k = 1:3
%                 for l = 1:3
%                     K_v1_vox = K_v1_vox + v1(i) * v1(j) * v1(k) * v1(l) * W(i,j,k,l);
%                 end
%             end
%         end
%     end
%     K_v1(v) = K_v1_vox / D_v1^2;
