function [strexp_D_all, strexp_gamma_all, biexp_f_all, biexp_Dfast_all, ...
    biexp_Dslow_all, dki_D_all, dki_K_all, monoexp_I, strexp_I, biexp_I, dki_I, ...
    stat_mu_all, stat_sigma_all, stat_I, gamma_kappa_all, gamma_theta_all, gamma_I, ...
    trunc_mu_all, trunc_sigma_all, trunc_D_upper_all, trunc_I, ...
    beta_alpha_all, beta_beta_all, beta_D0_all, beta_I] = ...
    fit_local_models( im, b, thresh_val, method, weighting)


    
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


sz = size(im);

if nargin < 3
    thresh_val = -inf;
end

if (nargin < 4) || isempty(method)% default method is linear least squares
    method = 'linear';
end

perform_weighting = false;
if nargin > 4
    perform_weighting = true;
    if length(weighting) ~= sz(end)
        disp('Weighting vector is wrong length')
        perform_weighting = false;
    end
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


% reshape to a matrix
im_matrix = reshape(im, [prod(sz(1:end-1)), sz(end)]);

% get the mask
if isscalar(thresh_val)
    mask = im_matrix(:,1) > thresh_val;
else
    mask = logical(thresh_val);
end

if ~any(mask(:))
    error('The mask is empty')
end

% take only the voxels inside the mask
I = im_matrix(mask(:), :);

if ~isreal(I)
    warning('Some voxels are complex. Taking magnitude.')
    I = abs(I);
end

b_value = b(1,1,:) + b(2,2,:) + b(3,3,:);
if min(b_value) > 50
    warning('Did you forget to include a b=0 image? The function doesn''t assume you divided it out')
end


% take the log of the image to linearise the equation
imlog = log(abs(I));
imlog(isinf(imlog)) = 0;

% Sort all b matrices in to a vector Bv=[Bxx,2*Bxy,2*Bxz,Byy,2*Byz,Bzz];
Bv=squeeze([b(1,1,:),2*b(1,2,:),2*b(1,3,:),b(2,2,:),2*b(2,3,:),b(3,3,:)])';



% If you have unit vectors and b values instead of b matrices, do this:
% b=zeros([3 3 size(unit_vectors,1)]);
% for i=1:size(unit_vectors,1),
%     b(:,:,i)=b_values(i)*unit_vectors(i,:)'*unit_vectors(i,:);
% end


% Add another column to Bv to handle the constant term:
% Slog = Bv * M + log(S0)
% becomes:
% Slog = [Bv, -1] * [M; -log(S0)]
Bv = [Bv, -ones(size(Bv,1), 1)];



% Fit the tensor model (linear, unconstrained)
%M = (Bv \ -imlog')';
% This is the same as (pinv(Bv) * (-imlog'))';
% also same as (pinv(Bv' * Bv) * Bv' * (-imlog'))';
% weighting (NxN) goes   ^ here      ^ and here
if perform_weighting

    W = eye(length(Bv));
    W(W == 1) = weighting(:);
    M = (pinv(Bv' * (W.^2) * Bv) * Bv' * (W.^2) * (-imlog'))';
else
    param.rician = 0;
    param.unique_weights = 0;
    param.verbose = 0;
    M = weighted_linear_fit(-imlog, Bv', @(z)exp(-z), param);
end

% convert log(S0) to S0
M(:,7) = exp(M(:,7));
M(isnan(M)) = 0;
M(isinf(M)) = 0;

Bv = Bv(:,1:6)';


lb = [0 -3E-3 -3E-3 0 -3E-3 0 0]; % cross terms are allowed to be -ve
ub = [zeros(1,6)+3E-3, inf];

M = bsxfun(@max, M, lb);
M = bsxfun(@min, M, ub);

% perform nonlinear fitting of DTI
if strcmp(method, 'nonlinear')

    % make sure you have doubles
    I = double(I);
    M = double(M);
    
    

    options = optimoptions('lsqcurvefit','Jacobian','on', 'DerivativeCheck', 'off', ...
        'display', 'off', 'TypicalX', mean(M,1));
    % M = [xx,xy,xz,yy,yz,zz,S0]
    
    
%     if (size(I,1) > 1000) && (matlabpool('size') == 0)
%         disp('Open matlabpool for parallel processing')
%     end
    
    parfor j = 1:size(I,1)

        M(j,:) = lsqcurvefit(@mono_exp_model, M(j,:), Bv, I(j,:), lb, ub, options);

    end
    

end
    
monoexp_I = zeros([prod(sz(1:end-1)), sz(end)]);
monoexp_I(mask,:) = bsxfun(@times, M(:,7), exp(-(M(:,1:6) * Bv)));
monoexp_I = reshape(monoexp_I, sz);



% initialise variables
strexp_D = zeros(size(M,1), 3);
strexp_gamma = zeros(size(M,1), 3);
biexp_f = zeros(size(M,1), 1);
biexp_Dfast = zeros(size(M,1), 3);
biexp_Dslow = zeros(size(M,1), 3);
dki_D = zeros(size(M,1),3);
dki_K = zeros(size(M,1),3);
stat_mu = zeros(size(M,1),3);
stat_sigma = zeros(size(M,1),3);
gamma_kappa = zeros(size(M,1),3);
gamma_theta = zeros(size(M,1),3);
trunc_mu = zeros(size(M,1),3);
trunc_sigma = zeros(size(M,1),3);
trunc_D_upper = zeros(size(M,1),1);
beta_alpha = zeros(size(M,1),3);
beta_beta = zeros(size(M,1),3);
beta_D0 = zeros(size(M,1),1);

strexp_I = zeros(size(M,1), sz(end));
biexp_I = zeros(size(M,1), sz(end));
dki_I = zeros(size(M,1), sz(end));
stat_I = zeros(size(M,1), sz(end));
gamma_I = zeros(size(M,1), sz(end));
trunc_I = zeros(size(M,1), sz(end));
beta_I = zeros(size(M,1), sz(end));

options = optimoptions('lsqcurvefit', 'display', 'off', 'MaxIter', 400);

strexp_lb = [0 0 0 0 0 0 0];
strexp_ub = [inf 10E-3 10E-3 10E-3 1 1 1];
biexp_lb = [0 0 0 0 0 0 0 0];
biexp_ub = [inf 10E-3 10E-3 10E-3 10E-3 10E-3 10E-3 1];
dki_lb = [0 0 0 0 0 0 0];
dki_ub = [inf 10E-3 10E-3 10E-3 3 3 3];
stat_lb = zeros(1,7) - inf;
stat_ub = zeros(1,7) + inf;
gamma_lb = zeros(1,7);
gamma_ub = zeros(1,7) + inf;

%parfor_progress(size(M,1));

disp('Beginning non-Gaussian models')


loop_interval = max(12, ceil(size(M,1) / 100));

for outer_loop = 1:100
    
    disp(outer_loop)
    
    loop_lower = loop_interval * (outer_loop - 1) + 1;
    loop_upper = min(size(M,1), loop_interval * outer_loop);
    
    M_loop = M(loop_lower:loop_upper, :);
    I_loop = I(loop_lower:loop_upper, :);

    
    loop_length = length(loop_lower:loop_upper);
    
    strexp_D_loop = zeros(loop_length, 3);
    strexp_gamma_loop = zeros(loop_length, 3);
    biexp_f_loop = zeros(loop_length, 1);
    biexp_Dfast_loop = zeros(loop_length, 3);
    biexp_Dslow_loop = zeros(loop_length, 3);
    dki_D_loop = zeros(loop_length,3);
    dki_K_loop = zeros(loop_length,3);
    stat_mu_loop = zeros(loop_length,3);
    stat_sigma_loop = zeros(loop_length,3);
    gamma_kappa_loop = zeros(loop_length,3);
    gamma_theta_loop = zeros(loop_length,3);
    trunc_mu_loop = zeros(loop_length,3);
    trunc_sigma_loop = zeros(loop_length,3);
    trunc_D_upper_loop = zeros(loop_length,1);
    beta_alpha_loop = zeros(loop_length,3);
    beta_beta_loop = zeros(loop_length,3);
    beta_D0_loop = zeros(loop_length,1);
    
    strexp_I_loop = zeros(loop_length, sz(end));
    biexp_I_loop = zeros(loop_length, sz(end));
    dki_I_loop = zeros(loop_length, sz(end));
    stat_I_loop = zeros(loop_length, sz(end));
    gamma_I_loop = zeros(loop_length, sz(end));
    trunc_I_loop = zeros(loop_length, sz(end));
    beta_I_loop = zeros(loop_length, sz(end));
    
    
    parfor i = 1:loop_length

        Mi = M_loop(i,:);
        Ii = double(I_loop(i,:));

        % The DiffusionTensor (Remember it is a symetric matrix,
        % thus for instance Dxy == Dyx)
        DiffusionTensor=[Mi(1) Mi(2) Mi(3); Mi(2) Mi(4) Mi(5); Mi(3) Mi(5) Mi(6)];

        % Calculate the eigenvalues and vectors, and sort the 
        % eigenvalues from small to large
        [EigenVectors,D]=eig(DiffusionTensor); 
        EigenValues=diag(D);
        [~,index]=sort(EigenValues, 'descend'); 
        EigenValues=EigenValues(index); 
        EigenVectors=EigenVectors(:,index);

        EigenValues = abs(EigenValues);


        % rotate the b-matrix so that it matches the reference frame
        % defined by the diffusion tensor (see De Santis 2011)

        b1 = zeros(size(b,3),1);
        b2 = zeros(size(b,3),1);
        b3 = zeros(size(b,3),1);

        for b_counter = 1:length(b1)
            b1(b_counter) = EigenVectors(:,1)' * b(:,:,b_counter) * EigenVectors(:,1);
            b2(b_counter) = EigenVectors(:,2)' * b(:,:,b_counter) * EigenVectors(:,2);
            b3(b_counter) = EigenVectors(:,3)' * b(:,:,b_counter) * EigenVectors(:,3);
        end
        b123 = [b1, b2, b3];


        %% fit the stretched exponential model
        x = double([Mi(7), EigenValues', 0.9, 0.9 0.9]); % [S0 D1 D2 D3 g1 g2 g3]

        % S0 exp( - D1 b1^g1 - D2 b2^g2 - D3 b3^g3 )

        mdl = @(x,b123) x(1) * exp(-x(2)*(b123(:,1).^x(5)) - x(3)*(b123(:,2).^x(6)) - x(4)*(b123(:,3).^x(7)));

        strexp_fit = lsqcurvefit(mdl, x, b123, Ii', strexp_lb, strexp_ub, options);


        strexp_D_loop(i,:) = strexp_fit(2:4);
        strexp_gamma_loop(i,:) = strexp_fit(5:7);

        strexp_I_loop(i,:) = mdl(strexp_fit, b123);

        %% fit the kurtosis model

        % this can be linearised!

        b_mag = sum(b123,2);
        b_lin = [ones(size(b123(b_mag < 5001,:),1), 1), -b123(b_mag < 5001,:), 1/6*b123(b_mag < 5001,:).^2];

        x = weighted_linear_fit(log(abs(Ii(b_mag < 5001)')), b_lin, @(x) exp(x));
        x(1) = exp(x(1)); % log(S0) -> S0
        %       S0      D1 D2 D3    K1 K2 K3
        %       1       2  3  4     5  6  7



        mdl = @(x,b_lin) x(1) * exp(x(2) * b_lin(:,2) + x(5) * b_lin(:,5) + ...
                                    x(3) * b_lin(:,3) + x(6) * b_lin(:,6) + ...
                                    x(4) * b_lin(:,4) + x(7) * b_lin(:,7));

        % check that the linear fit worked...
        SSE_dti = sum((Ii(b_mag < 5001) -  Mi(7) * exp(-(Mi(1:6) * Bv(:,b_mag < 5001)))).^2);
        SSE_dki = sum((Ii(b_mag < 5001)' - mdl(x, b_lin)).^2);

        if SSE_dti < SSE_dki
            % the linear fit failed, use the DTI fit
            x = [Mi(7), EigenValues', [0.3, 0.5, 0.5] .* EigenValues'.^2 ];

        end


        x = max(x, dki_lb);
        x = min(x, dki_ub);


        dki_fit = lsqcurvefit(mdl, x, b_lin, Ii(b_mag < 5001)', dki_lb, dki_ub, options);


        dki_D_loop(i,:) = dki_fit(2:4);
        dki_K_loop(i,:) = dki_fit(5:7) ./ (dki_fit(2:4).^2);

        % use all of the b values to reconstruct signal
        b_lin = [ones(size(b123,1), 1), -b123, 1/6*b123.^2];

        dki_I_loop(i,:) = mdl(dki_fit, b_lin);


        %% fit the bi-exponential model
        x = double([Mi(7), EigenValues', EigenValues'/10, 0.9]);
        %           1       2 3 4          5 6 7          8
        % ^ this works better than the clever way below 

    %     b_lin = [ones(size(b123,1), 1), -b123];
    %     
    %     % do a least squares fit of the b values < 2100
    %     b_mag = sum(b123,2);
    %     low_idx = b_mag <= 1001;
    %     b123_low = b123(low_idx, :);
    %     b_lin_low = [ones(size(b123_low,1), 1), -b123_low];
    % 
    %     low_fit = b_lin_low \ log(abs(double(Ii(low_idx))'));
    %     I_low = exp(low_fit' * b_lin');
    %     %low_fit(1) = exp(low_fit(1)); % log(S0) -> S0
    % 
    %     high_idx = b_mag > 6000;
    %     b123_high = b123(high_idx, :);
    %     b_lin_high = [ones(size(b123_high,1), 1), -b123_high];
    % 
    %     high_fit = b_lin_high \ log(abs(double(Ii(high_idx))'));
    %     I_high = exp(high_fit' * b_lin');
    %     %high_fit(1) = exp(high_fit(1)); % log(S0) -> S0
    %   
    %     % now that I have two tensors describing the signal, I just need to
    %     % find their ratios
    %     
    %     highlow = [I_low; I_high];
    %     W = Ii' \ highlow';
    %    
    %     S0 = exp(low_fit(1));
    %     f = W(1) / sum(W);
    %     x = [S0 low_fit(2:4)' high_fit(2:4)' f];
    %     
    %     % this should really be S0 = w1 S01 + w2 S02 but the low b values
    %     % estimate S0 much better, so I'm using that.


        mdl = @(x,b123) x(1) * (x(8) * exp(-x(2) * b123(:,1)) + (1 - x(8)) * exp(-x(5) * b123(:,1))) .* ...
                               (x(8) * exp(-x(3) * b123(:,2)) + (1 - x(8)) * exp(-x(6) * b123(:,2))) .* ...
                               (x(8) * exp(-x(4) * b123(:,3)) + (1 - x(8)) * exp(-x(7) * b123(:,3)));

        biexp_fit = lsqcurvefit(mdl, x, b123, Ii', biexp_lb, biexp_ub, options);                   


        biexp_f_loop(i,:) = biexp_fit(8);
        biexp_Dfast_loop(i,:) = biexp_fit(2:4);
        biexp_Dslow_loop(i,:) = biexp_fit(5:7);

        biexp_I_loop(i,:) = mdl(biexp_fit, b123);

        %% try a tri-exponential model
        %x = [biexp_fit(1:7)
    %     S0 = biexp_fit(1) ^ (1/3);
    %     x = double([EigenValues'*1.5, EigenValues', EigenValues'/10, 0.25*S0, 0.5*S0, 0.25*S0]);
    %     %               1 2 3            4 5 6          7 8 9          10      11       12
    % 
    %     mdl = @(x,b123) (x(10) * exp(-x(1) * b123(:,1)) + x(11) * exp(-x(4) * b123(:,1)) + x(12) * exp(-x(7)*b123(:,1)) ) .* ...
    %                     (x(10) * exp(-x(2) * b123(:,2)) + x(11) * exp(-x(5) * b123(:,2)) + x(12) * exp(-x(8)*b123(:,2)) ) .* ...
    %                     (x(10) * exp(-x(3) * b123(:,3)) + x(11) * exp(-x(6) * b123(:,3)) + x(12) * exp(-x(9)*b123(:,3)) );
    %                        
    %                        triexp_lb = zeros(1,12);
    %                        triexp_ub = zeros(1,12) + inf;
    %     triexp_fit = lsqcurvefit(mdl, x, b123, double(Ii)', triexp_lb, triexp_ub, options);                   

    % note: this was extremely sensitive to initial conditions! 





        %% statistical model (truncated Gaussian)

        % S0 D1 D2 D3 sigma1 sigma2 sigma3
        % 1  2  3  4    5      6      7

        mdl = @(x,b123) x(1) * ...
            (1 + erf(x(2)/(x(5)*sqrt(2)) - b123(:,1)*x(5)/sqrt(2))) ./ (1 + erf(x(2)/(x(5) * sqrt(2)))) .*  ...
            (1 + erf(x(3)/(x(6)*sqrt(2)) - b123(:,2)*x(6)/sqrt(2))) ./ (1 + erf(x(3)/(x(6) * sqrt(2)))) .*  ...
            (1 + erf(x(4)/(x(7)*sqrt(2)) - b123(:,3)*x(7)/sqrt(2))) ./ (1 + erf(x(4)/(x(7) * sqrt(2)))) .*  ...
            exp(-b123(:,1) * x(2) + 1/2 * b123(:,1).^2 * x(5).^2 + ...
                -b123(:,2) * x(3) + 1/2 * b123(:,2).^2 * x(6).^2 + ...
                -b123(:,3) * x(4) + 1/2 * b123(:,3).^2 * x(7).^2);



        % initialise with dki fit
        x = dki_fit;
        x(5:7) = sqrt(1/3 * dki_D_loop(i,:).^2 .* dki_K_loop(i,:));

        stat_fit = lsqcurvefit(mdl, x, b123, Ii', stat_lb, stat_ub, options);

        stat_mu_loop(i,:) = stat_fit(2:4);
        stat_sigma_loop(i,:) = stat_fit(5:7);

        stat_I_loop(i,:) = mdl(stat_fit, b123);

        %% statistical model (gamma)
        % I(b) = I0 (1 + b sigma^2 / D) ^ (- D^2 / sigma^2)
        % I0 k1 k2 k3 theta1 theta2 theta3
        % 1  2  3  4  5      6      7


        mdl = @(x, b123) x(1) * ...
            (1 + b123(:,1) * x(5)) .^ (- x(2)) .* ...
            (1 + b123(:,2) * x(6)) .^ (- x(3)) .* ...
            (1 + b123(:,3) * x(7)) .^ (- x(4));

        x = [dki_fit(1), 3 ./ (2*dki_K_loop(i,:)), dki_D_loop(i,:) .* (2*dki_K_loop(i,:)) / 3];

        gamma_fit =  lsqcurvefit(mdl, x, b123, Ii', gamma_lb, gamma_ub, options);

        gamma_kappa_loop(i,:) = gamma_fit(2:4);
        gamma_theta_loop(i,:) = gamma_fit(5:7);
        gamma_I_loop(i,:) = mdl(gamma_fit, b123);


        %% statistical model (doubly truncated gaussian)
         % this one is very good, but I don't want to be criticised for
         % proposing a new model when many others already exist, and this
         % isn't the point of the paper

        % S0 D1 D2 D3 sigma1 sigma2 sigma3 D_upper
        % 1  2  3  4    5      6      7     8


            % x1
        mdl = @(x,b123) x(1) * ...
            (normcdf(x(8), x(2) - b123(:,1)*x(5)^2, x(5)) - normcdf(0, x(2) - b123(:,1)*x(5)^2, x(5))) ...
          / (normcdf(x(8), x(2), x(5)) - normcdf(0, x(2), x(5)) + eps) .* ...
            (normcdf(x(8), x(3) - b123(:,2)*x(6)^2, x(6)) - normcdf(0, x(3) - b123(:,2)*x(6)^2, x(6))) ...
          / (normcdf(x(8), x(3), x(6)) - normcdf(0, x(3), x(6)) + eps) .* ...
            (normcdf(x(8), x(4) - b123(:,3)*x(7)^2, x(7)) - normcdf(0, x(4) - b123(:,3)*x(7)^2, x(7))) ...
          / (normcdf(x(8), x(4), x(7)) - normcdf(0, x(4), x(7)) + eps) .* ...
            exp(-b123(:,1) * x(2) + 1/2 * b123(:,1).^2 * x(5).^2 + ...
                -b123(:,2) * x(3) + 1/2 * b123(:,2).^2 * x(6).^2 + ...
                -b123(:,3) * x(4) + 1/2 * b123(:,3).^2 * x(7).^2);


        % initialise with statistical model fit
        x = [stat_fit, 2.3E-3];
        
        tst = mdl(x, b123);
        if any(isnan(tst))
            disp(i)
            disp('Trunc fit gives nans')
            x = [dki_fit, 2.3E-3];
            x(5:7) = sqrt(1/3 * dki_D_loop(i,:).^2 .* dki_K_loop(i,:));
        end

        trunc_fit = lsqcurvefit(mdl, x, b123, Ii', [stat_lb 0], [stat_ub 5E-3], options);


        trunc_mu_loop(i,:) = trunc_fit(2:4);
        trunc_sigma_loop(i,:) = trunc_fit(5:7);

        trunc_D_upper_loop(i,:) = trunc_fit(8);

        trunc_I_loop(i,:) = mdl(trunc_fit, b123);


        %% beta distribution

        % S = S0 Fch( alpha, alpha + beta, -bD0)

        % S0 alpha1 alpha2 alpha3 beta1 beta2 beta3 D0
        % 1   2      3      4      5     6     7    8

        mdl = @(x, b123) beta_function(x, b123);

    %     mdl = @(x, b123) x(1) * ...
    %         real(arrayfun(@KummerComplex, repmat(x(2), [size(b123,1), 1]), ...
    %             repmat(x(2) + x(5), [size(b123,1), 1]), -b123(:,1)*x(8))) .* ...
    %         real(arrayfun(@KummerComplex, repmat(x(3), [size(b123,1), 1]), ...
    %             repmat(x(3) + x(6), [size(b123,1), 1]), -b123(:,2)*x(8))) .* ...
    %         real(arrayfun(@KummerComplex, repmat(x(4), [size(b123,1), 1]), ...
    %             repmat(x(4) + x(7), [size(b123,1), 1]), -b123(:,3)*x(8)));


        % alpha = beta = 1.5 
        % x = [gamma_fit(1), zeros(1,6) + 1.5, 2.3E-3];

        Dmax = 2.3E-3;

        x = [dki_fit(1), 3./(2*dki_K_loop(i,:))-dki_D_loop(i,:)/Dmax - 3 ./ (2*dki_K_loop(i,:)) .* dki_D_loop(i,:) / Dmax, ...
             1/2*(3 ./ (2*dki_K_loop(i,:)).*dki_D_loop(i,:)/Dmax + 3./(2*dki_K_loop(i,:))*Dmax./dki_D_loop(i,:) + dki_D_loop(i,:)/Dmax - 6 ./ (2*dki_K_loop(i,:)) - 1), Dmax];



        beta_fit = lsqcurvefit(mdl, x, b123, Ii', [stat_lb 0], [stat_ub 5E-3], options);

        beta_alpha_loop(i,:) = beta_fit(2:4);
        beta_beta_loop(i,:) = beta_fit(5:7);

        beta_D0_loop(i,:) = beta_fit(8);

        beta_I_loop(i,:) = mdl(beta_fit, b123);


        %parfor_progress;

    end
    
    
    strexp_D(loop_lower:loop_upper,:) = strexp_D_loop;
    strexp_gamma(loop_lower:loop_upper,:) = strexp_gamma_loop;
    biexp_f(loop_lower:loop_upper,:) = biexp_f_loop;
    biexp_Dfast(loop_lower:loop_upper,:) = biexp_Dfast_loop;
    biexp_Dslow(loop_lower:loop_upper,:) = biexp_Dslow_loop;
    dki_D(loop_lower:loop_upper,:) = dki_D_loop;
    dki_K(loop_lower:loop_upper,:) = dki_K_loop;
    stat_mu(loop_lower:loop_upper,:) = stat_mu_loop;
    stat_sigma(loop_lower:loop_upper,:) = stat_sigma_loop;
    gamma_kappa(loop_lower:loop_upper,:) = gamma_kappa_loop;
    gamma_theta(loop_lower:loop_upper,:) = gamma_theta_loop;
    trunc_mu(loop_lower:loop_upper,:) = trunc_mu_loop;
    trunc_sigma(loop_lower:loop_upper,:) = trunc_sigma_loop;
    trunc_D_upper(loop_lower:loop_upper,:) = trunc_D_upper_loop;
    beta_alpha(loop_lower:loop_upper,:) = beta_alpha_loop;
    beta_beta(loop_lower:loop_upper,:) = beta_beta_loop;
    beta_D0(loop_lower:loop_upper,:) = beta_D0_loop;
    
    strexp_I(loop_lower:loop_upper,:) = strexp_I_loop;
    biexp_I(loop_lower:loop_upper,:) = biexp_I_loop;
    dki_I(loop_lower:loop_upper,:) = dki_I_loop;
    stat_I(loop_lower:loop_upper,:) = stat_I_loop;
    gamma_I(loop_lower:loop_upper,:) = gamma_I_loop;
    trunc_I(loop_lower:loop_upper,:) = trunc_I_loop;
    beta_I(loop_lower:loop_upper,:) = beta_I_loop;

end

%parfor_progress(0);


% quick hack because reshape sz needs to be at least length 2
if length(sz) == 2
    sz = [sz(1), 1, sz(2)];
end

strexp_D_all = zeros([prod(sz(1:end-1)), 3]);
strexp_D_all(mask(:),:) = strexp_D;
strexp_D_all = reshape(strexp_D_all, [sz(1:end-1), 3]);
strexp_gamma_all = zeros([prod(sz(1:end-1)), 3]);
strexp_gamma_all(mask(:),:) = strexp_gamma;
strexp_gamma_all = reshape(strexp_gamma_all, [sz(1:end-1), 3]);
biexp_f_all = zeros([prod(sz(1:end-1)), 1]);
biexp_f_all(mask(:)) = biexp_f;
biexp_f_all = reshape(biexp_f_all, sz(1:end-1));
biexp_Dfast_all = zeros([prod(sz(1:end-1)), 3]);
biexp_Dfast_all(mask(:),:) = biexp_Dfast;
biexp_Dfast_all = reshape(biexp_Dfast_all, [sz(1:end-1), 3]);
biexp_Dslow_all = zeros([prod(sz(1:end-1)), 3]);
biexp_Dslow_all(mask(:),:) = biexp_Dslow;
biexp_Dslow_all = reshape(biexp_Dslow_all, [sz(1:end-1), 3]);
dki_D_all = zeros([prod(sz(1:end-1)), 3]);
dki_D_all(mask(:),:) = dki_D;
dki_D_all = reshape(dki_D_all, [sz(1:end-1) 3]);
dki_K_all = zeros([prod(sz(1:end-1)), 3]);
dki_K_all(mask(:),:) = dki_K;
dki_K_all = reshape(dki_K_all, [sz(1:end-1), 3]);
stat_mu_all = zeros([prod(sz(1:end-1)), 3]);
stat_mu_all(mask(:),:) = stat_mu;
stat_mu_all = reshape(stat_mu_all, [sz(1:end-1), 3]);
stat_sigma_all = zeros([prod(sz(1:end-1)), 3]);
stat_sigma_all(mask(:),:) = stat_sigma;
stat_sigma_all = reshape(stat_sigma_all, [sz(1:end-1), 3]);
gamma_kappa_all = zeros([prod(sz(1:end-1)), 3]);
gamma_kappa_all(mask(:),:) = gamma_kappa;
gamma_kappa_all = reshape(gamma_kappa_all, [sz(1:end-1), 3]);
gamma_theta_all = zeros([prod(sz(1:end-1)), 3]);
gamma_theta_all(mask(:),:) = gamma_theta;
gamma_theta_all = reshape(gamma_theta_all, [sz(1:end-1), 3]);    
trunc_mu_all = zeros([prod(sz(1:end-1)), 3]);
trunc_mu_all(mask(:),:) = trunc_mu;
trunc_mu_all = reshape(trunc_mu_all, [sz(1:end-1), 3]);
trunc_sigma_all = zeros([prod(sz(1:end-1)), 3]);
trunc_sigma_all(mask(:),:) = trunc_sigma;
trunc_sigma_all = reshape(trunc_sigma_all, [sz(1:end-1), 3]);
trunc_D_upper_all = zeros([prod(sz(1:end-1)), 1]);
trunc_D_upper_all(mask(:),:) = trunc_D_upper;
trunc_D_upper_all = reshape(trunc_D_upper_all, [sz(1:end-1), 1]);
beta_alpha_all = zeros([prod(sz(1:end-1)), 3]);
beta_alpha_all(mask(:),:) = beta_alpha;
beta_alpha_all = reshape(beta_alpha_all, [sz(1:end-1), 3]);
beta_beta_all = zeros([prod(sz(1:end-1)), 3]);
beta_beta_all(mask(:),:) = beta_beta;
beta_beta_all = reshape(beta_beta_all, [sz(1:end-1), 3]);
beta_D0_all = zeros([prod(sz(1:end-1)), 1]);
beta_D0_all(mask(:),:) = beta_D0;
beta_D0_all = reshape(beta_D0_all, [sz(1:end-1), 1]);


strexp_I_vect = zeros([prod(sz(1:end-1)), sz(end)]);
strexp_I_vect(mask(:),:) = strexp_I;
strexp_I = reshape(strexp_I_vect, [sz(1:end-1), sz(end)]);
biexp_I_vect = zeros([prod(sz(1:end-1)), sz(end)]);
biexp_I_vect(mask(:),:) = biexp_I;
biexp_I = reshape(biexp_I_vect, [sz(1:end-1), sz(end)]);
dki_I_vect = zeros([prod(sz(1:end-1)), sz(end)]);
dki_I_vect(mask(:),:) = dki_I;
dki_I = reshape(dki_I_vect, [sz(1:end-1), sz(end)]);
stat_I_vect = zeros([prod(sz(1:end-1)), sz(end)]);
stat_I_vect(mask(:),:) = stat_I;
stat_I = reshape(stat_I_vect, [sz(1:end-1), sz(end)]);
gamma_I_vect = zeros([prod(sz(1:end-1)), sz(end)]);
gamma_I_vect(mask(:),:) = gamma_I;
gamma_I = reshape(gamma_I_vect, [sz(1:end-1), sz(end)]);
trunc_I_vect = zeros([prod(sz(1:end-1)), sz(end)]);
trunc_I_vect(mask(:),:) = trunc_I;
trunc_I = reshape(trunc_I_vect, [sz(1:end-1), sz(end)]);
beta_I_vect = zeros([prod(sz(1:end-1)), sz(end)]);
beta_I_vect(mask(:),:) = beta_I;
beta_I = reshape(beta_I_vect, [sz(1:end-1), sz(end)]);

end

% the function to fit inside the parfor loop
function [F, J] = mono_exp_model(x, Bv)

F = x(7) * exp(-(x(1:6) * Bv));
    
if nargout > 1 % Jacobian
    J = zeros(length(F), length(x));
    
    J(:,1) =  x(7) * -Bv(1,:) .* exp(-(x(1:6) * Bv));
    J(:,2) =  x(7) * -Bv(2,:) .* exp(-(x(1:6) * Bv));
    J(:,3) =  x(7) * -Bv(3,:) .* exp(-(x(1:6) * Bv));
    J(:,4) =  x(7) * -Bv(4,:) .* exp(-(x(1:6) * Bv));
    J(:,5) =  x(7) * -Bv(5,:) .* exp(-(x(1:6) * Bv));
    J(:,6) =  x(7) * -Bv(6,:) .* exp(-(x(1:6) * Bv));        
    J(:,7) = exp(-(x(1:6) * Bv));
end

end

function y = beta_function(x, b123)
% I found that the KummerComplex function on the mathworks website worked, but was too slow.

y1 = myKummerComplex(x(2), x(2) + x(5), -b123(:,1)*x(8));
y2 = myKummerComplex(x(3), x(3) + x(6), -b123(:,2)*x(8));
y3 = myKummerComplex(x(4), x(4) + x(7), -b123(:,3)*x(8));

    
y =  real(x(1) * y1 .* y2 .* y3);

end
    
    
    

