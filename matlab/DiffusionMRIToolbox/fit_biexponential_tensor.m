function [M, ADC1, ADC2, FA1, FA2, VectorF1, VectorF2] = fit_biexponential_tensor(I, bval, mask)
%FIT_BIEXPONENTIAL_TENSOR Fit the bi-exponential tensor model
% fits the function S = S0_fast * e^(b*D_fast) + S0_slow e^(b*D_slow)
%
%
% Inputs:
%   I is the input image, of any dimensionality, with the diffusion 
%   scans in the last dimension
%     
%   BVAL is the b-matrix, of size [3 3 N]. If the measured b matrix isn't 
%   available, it can be approximated by:
%   bval(:,:,n)=b_values(n)*unit_vectors(n,:)'*unit_vectors(n,:);
%
%   MASK is a mask of the voxels for non-linear fitting, ADC, FA, etc
%   computation. The voxels outside the mask are fit with a very crude
%   estimate of the model.
%
% Outputs:
%   M is the model coefficients, of size [row col slice 14]. The last
%   dimension is as follows:
%       [S0_fast  D1_fast ... D6_fast S0_slow D1_slow ... D6_slow]
%       [1        2       ... 7       8       9       ... 14]
%
%   ADC1, ADC2 are the ADC maps for the fast and slow component
%   respectively
%
%   FA1, FA2 are the FA maps (bear in mind that FA2 might be very noisy)
%
%   VECTORF1, VECTORF2 are eigenvector fields for fast and slow components
%
%
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
narginchk(2,3);
nargoutchk(0, 7);

sz = size(I);

if nargin < 3
    mask = ones(sz(1:end-1));
end

mask = logical(mask);


% first of all, fit a normal diffusion tensor to anything below 1000
% (Ingliss 2001) for an estimate of the fast component.

bval_trace = squeeze(bval(1,1,:) + bval(2,2,:) + bval(3,3,:));
fast_bvals = bval_trace < 1000;

DT_fast = reshape(fit_DT(I(:,:,:,fast_bvals), bval(:,:,fast_bvals)), [prod(sz(1:3)), 7]);

% likewise, fit another tensor to anything about 1000 for the slow
% component
slow_bvals = bval_trace > 1000;

DT_slow = reshape(fit_DT(I(:,:,:,slow_bvals), bval(:,:,slow_bvals)), [prod(sz(1:3)), 7]);


I_fast = reshape(dt2image(DT_fast, bval), [prod(sz(1:3)), sz(4)]);
I_slow = reshape(dt2image(DT_slow, bval), [prod(sz(1:3)), sz(4)]);

Ivector = reshape(I, [prod(sz(1:3)), sz(4)]);

DT_fast = DT_fast(:, [7, 1:6]); % put S0 first
DT_slow = DT_slow(:, [7, 1:6]);



M = zeros(size(Ivector,1), 14);

Bv=squeeze([bval(1,1,:),2*bval(1,2,:),2*bval(1,3,:),bval(2,2,:),2*bval(2,3,:),bval(3,3,:)]);

if (matlabpool('size') == 0) && (size(M,1) > 100) % if there is a decent number of voxels
    disp('Open matlab pool if you want to use parallel processing')
end

if min(mask(:)) == 0

    % in unimportant voxels, just take a linear combination of fast and slow
    I_linfit = Ivector(~mask(:),:);
    DT_fast_linfit = DT_fast(~mask(:), :);
    DT_slow_linfit = DT_slow(~mask(:), :);
    I_fast_linfit = I_fast(~mask(:), :);
    I_slow_linfit = I_slow(~mask(:), :);
    coeffs_linfit = M(~mask(:), :);

    tic
    parfor i = 1:size(I_linfit,1)

        fastslow = [I_fast_linfit(i,:); I_slow_linfit(i,:)]';
        w = fastslow \I_linfit(i,:)';
        coeffs_linfit(i,:) = [w(1) * DT_fast_linfit(i,1), DT_fast_linfit(i, 2:7), ...
                              w(2) * DT_slow_linfit(i,1), DT_slow_linfit(i, 2:7)];

    end
    t = toc;

    disp(['Linear fit complete: finished ' num2str(size(I_linfit,1)) ' voxels in ' num2str(t) ' seconds'])
    disp(['That is ' num2str(size(I_linfit,1)/t) ' voxels per second'])
    disp(' ')
    M(~mask(:), :) = coeffs_linfit;
end

if max(mask(:)) == 1

    % in the important voxels, do a nonlinear fit

    I_nlfit = double(Ivector(mask(:),:));
    DT_fast_nlfit = double(DT_fast(mask(:), :));
    DT_slow_nlfit = double(DT_slow(mask(:), :));
    I_fast_nlfit = double(I_fast(mask(:), :));
    I_slow_nlfit = double(I_slow(mask(:), :));
    coeffs_nlfit = double(M(mask(:), :));


    options = optimoptions('lsqcurvefit','Jacobian','on', 'DerivativeCheck', 'off', 'display', 'off');
    lb = [0 0 -3E-3 -3E-3 0 -3E-3 0 0 0 -3E-3 -3E-3 0 -3E-3 0]; % cross terms are allowed to be -ve
    ub = [inf, zeros(1,6)+3E-3, inf, zeros(1,6)+3E-3];

    tic
    parfor i = 1:size(I_nlfit,1)

        fastslow = [I_fast_nlfit(i,:); I_slow_nlfit(i,:)]';
        w = fastslow \I_nlfit(i,:)';

        x_guess = [w(1) * DT_fast_nlfit(i,1), DT_fast_nlfit(i, 2:7), ...
                   w(2) * DT_slow_nlfit(i,1), DT_slow_nlfit(i, 2:7)];

        coeffs_nlfit(i,:) = lsqcurvefit(@bi_exp_model, x_guess, Bv, I_nlfit(i,:), lb, ub, options);

    end

    t = toc;

    disp(['Non-linear fit complete: finished ' num2str(size(I_nlfit,1)) ' voxels in ' num2str(t) ' seconds'])
    disp(['That is ' num2str(size(I_nlfit,1)/t) ' voxels per second'])

    M(mask(:), :) = coeffs_nlfit;

    % fix the voxels where the fast and slow component have swapped (if any)
    fs_swap = M(:,8) > M(:,1);
    M(fs_swap, :) = M(fs_swap, [8:14, 1:7]);



    % compute FA, ADC, etc for both tensors in the nonlinear fitted region only
    M_nl = M(mask,:);


    % initialise variables
    FA = zeros(size(M_nl,1), 1);
    ADC = zeros(size(M_nl,1),1);
    VectorF = zeros(size(M_nl,1), 3);
    VectorF2 = zeros(size(M_nl,1), 3);
    VectorF3 = zeros(size(M_nl,1), 3);
    EigVals = zeros(size(M_nl,1), 3);

    for i = 1:size(M_nl,1)


        Mi = M_nl(i,2:7);

        % The DiffusionTensor (Remember it is a symetric matrix,
        % thus for instance Dxy == Dyx)
        DiffusionTensor=[Mi(1) Mi(2) Mi(3); Mi(2) Mi(4) Mi(5); Mi(3) Mi(5) Mi(6)];

        % Calculate the eigenvalues and vectors, and sort the 
        % eigenvalues from small to large
        [EigenVectors,D]=eig(DiffusionTensor); 
        EigenValues=diag(D);
        [~,index]=sort(EigenValues); 
        EigenValues=EigenValues(index); 
        EigenVectors=EigenVectors(:,index);
        EigenValues_old=EigenValues;

        % Regulating of the eigen values (negative eigenvalues are
        % due to noise and other non-idealities of MRI)
        EigenValues=abs(EigenValues);

        % Apparent Diffuse Coefficient
        ADCv=(EigenValues(1)+EigenValues(2)+EigenValues(3))/3;

        % Fractional Anistropy (2 different definitions exist)
        % First FA definition:
        %FAv=(1/sqrt(2))*( sqrt((EigenValues(1)-EigenValues(2)).^2+(EigenValues(2)-EigenValues(3)).^2+(EigenValues(1)-EigenValues(3)).^2)./sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2) );
        % Second FA definition:
        FA(i)=sqrt(1.5)*( sqrt((EigenValues(1)-ADCv).^2+(EigenValues(2)-ADCv).^2+(EigenValues(3)-ADCv).^2)./sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2) );
        ADC(i)=ADCv;
        VectorF(i,:)=EigenVectors(:,3)*EigenValues_old(3);
        VectorF2(i,:)=EigenVectors(:,2)*EigenValues_old(2);
        VectorF3(i,:)=EigenVectors(:,1)*EigenValues_old(1);

        EigVals(i,:) = EigenValues_old;

    end

    M = zeros([prod(sz(1:end-1)), 14]);
    M(mask(:), :) = M_nl;
    M = reshape(M, [sz(1:end-1), 14]);

    FA1 = zeros(sz(1:end-1));
    FA1(mask(:)) = FA;
    ADC1 = zeros(sz(1:end-1));
    ADC1(mask(:)) = ADC;
    VectorF1 = zeros([prod(sz(1:end-1)), 3]);
    VectorF1(mask(:),:) = VectorF;
    VectorF1 = reshape(VectorF1, [sz(1:end-1), 3]);


    for i = 1:size(M_nl,1)


        Mi = M_nl(i,9:14);

        % The DiffusionTensor (Remember it is a symetric matrix,
        % thus for instance Dxy == Dyx)
        DiffusionTensor=[Mi(1) Mi(2) Mi(3); Mi(2) Mi(4) Mi(5); Mi(3) Mi(5) Mi(6)];

        % Calculate the eigenvalues and vectors, and sort the 
        % eigenvalues from small to large
        [EigenVectors,D]=eig(DiffusionTensor); 
        EigenValues=diag(D);
        [~,index]=sort(EigenValues); 
        EigenValues=EigenValues(index); 
        EigenVectors=EigenVectors(:,index);
        EigenValues_old=EigenValues;

        % Regulating of the eigen values (negative eigenvalues are
        % due to noise and other non-idealities of MRI)
        EigenValues=abs(EigenValues);

        % Apparent Diffuse Coefficient
        ADCv=(EigenValues(1)+EigenValues(2)+EigenValues(3))/3;

        % Fractional Anistropy (2 different definitions exist)
        % First FA definition:
        %FAv=(1/sqrt(2))*( sqrt((EigenValues(1)-EigenValues(2)).^2+(EigenValues(2)-EigenValues(3)).^2+(EigenValues(1)-EigenValues(3)).^2)./sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2) );
        % Second FA definition:
        FA(i)=sqrt(1.5)*( sqrt((EigenValues(1)-ADCv).^2+(EigenValues(2)-ADCv).^2+(EigenValues(3)-ADCv).^2)./sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2) );
        ADC(i)=ADCv;
        VectorF(i,:)=EigenVectors(:,3)*EigenValues_old(3);
        VectorF2(i,:)=EigenVectors(:,2)*EigenValues_old(2);
        VectorF3(i,:)=EigenVectors(:,1)*EigenValues_old(1);

        EigVals(i,:) = EigenValues_old;

    end

    FA2 = zeros(sz(1:end-1));
    FA2(mask(:)) = FA;
    ADC2 = zeros(sz(1:end-1));
    ADC2(mask(:)) = ADC;
    VectorF2 = zeros([prod(sz(1:end-1)), 3]);
    VectorF2(mask(:),:) = VectorF;
    VectorF2 = reshape(VectorF2, [sz(1:end-1), 3]);



end


M = reshape(M, [sz(1:3), 14]);

end




% the function to fit inside the parfor loop
function [F, J] = bi_exp_model(x, Bv)


F = x(1) * exp(-(x(2:7) * Bv)) + ...
    x(8) * exp(-(x(9:14) * Bv));
    
if nargout > 1 % Jacobian (speeds up by 3x)
    J = zeros(length(F), length(x));
    J(:,1) = exp(-(x(2:7) * Bv));
    J(:,2) =  x(1) * -Bv(1,:) .* exp(-(x(2:7) * Bv));
    J(:,3) =  x(1) * -Bv(2,:) .* exp(-(x(2:7) * Bv));
    J(:,4) =  x(1) * -Bv(3,:) .* exp(-(x(2:7) * Bv));
    J(:,5) =  x(1) * -Bv(4,:) .* exp(-(x(2:7) * Bv));
    J(:,6) =  x(1) * -Bv(5,:) .* exp(-(x(2:7) * Bv));
    J(:,7) =  x(1) * -Bv(6,:) .* exp(-(x(2:7) * Bv));
        
    J(:,8) = exp(-(x(9:14) * Bv));
    J(:,9) =  x(8) * -Bv(1,:) .* exp(-(x(9:14) * Bv));
    J(:,10) =  x(8) * -Bv(2,:) .* exp(-(x(9:14) * Bv));
    J(:,11) =  x(8) * -Bv(3,:) .* exp(-(x(9:14) * Bv));
    J(:,12) =  x(8) * -Bv(4,:) .* exp(-(x(9:14) * Bv));
    J(:,13) =  x(8) * -Bv(5,:) .* exp(-(x(9:14) * Bv));
    J(:,14) =  x(8) * -Bv(6,:) .* exp(-(x(9:14) * Bv));
end
    
    
end




