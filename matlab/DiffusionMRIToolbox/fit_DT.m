function [ DT, FA, ADC, VectorField, EigVals] = fit_DT( im, b, thresh_val, method)
% FIT_DT    Fits the diffusion tensor model voxelwise to an image
%           S = S0 exp(-bD)
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
%       'weighted linear' fits the tensor with a weighted linear regression
%       (also quite fast, weights are automatically assigned)
%       'nonlinear' fits the tensor in the signal domain (recommended for
%       high b-values or when SNR is low)
%
%
% Outputs:
%
%   DT is the diffusion tensor coefficients [xx,xy,xz,yy,yz,zz,S0]
% 
%   FA is the fractional anisotropy
% 
%   ADC is the apparent diffusion coefficient
% 
%   VECTORFIELD is the primary, secondary and tertiary unit vectors, concatenated
%	in the last dimension ([r c s xyz v1v2v3])
%
%	EIGVALS are the eigenvalues corresponding to the vector field, concatenated
%	in the last dimension
%
%
% An alternative method is implemented in fit_DT_YHd_method.
%
% See also dt2image.

    
% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright © 2014-2015 University of Oxford
% Version: 0.1.10
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
nargoutchk(0, 5);

sz = size(im);

if nargin < 3
    thresh_val = -inf;
end

if (nargin < 4) || isempty(method)% default method is weighted linear least squares
    method = 'weighted linear';
end


% check that a valid method has been given
if ~max([strcmp(method, 'linear'), strcmp(method, 'weighted linear'), strcmp(method, 'nonlinear')])
    disp('Unrecognised method, performing weighted linear fitting')
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
if strcmp(method, 'linear')
    M = (Bv \ -imlog')';%
else
    % this will always yield a lower solution than the linear fitting
    % because there is a check in there that replaces any voxels that fail
    % with the linear fit version!
    
    param.rician = 1;
    param.unique_weights = 1;
    param.verbose = 0;
    M = weighted_linear_fit(-imlog, Bv', @(z)exp(-z), param);
end

% convert log(S0) to S0
M(:,7) = exp(M(:,7));
M(isnan(M)) = 0;
M(isinf(M)) = 0;

% perform nonlinear fitting (if required)
if strcmp(method, 'nonlinear')

    % make sure you have doubles
    I = double(I);
    M = double(M);
    
    Bv = Bv(:,1:6)';

    options = optimoptions('lsqcurvefit','Jacobian', 'on', 'DerivativeCheck', 'off', ...
        'display', 'off', 'TypicalX', mean(M,1));
    % M = [xx,xy,xz,yy,yz,zz,S0]
    lb = [0 -3E-3 -3E-3 0 -3E-3 0 0]; % cross terms are allowed to be -ve
    ub = [zeros(1,6)+3E-3, inf];

    v = ver('MATLAB');
    
    if (size(I,1) > 1000) && (str2double(v.Version) <= 8.1)
        if (matlabpool('size') == 0)
            disp('You should open matlabpool if you want parallel processing')
        end
    end
    
    
    
    parfor i = 1:size(I,1)

        M(i,:) = lsqcurvefit(@mono_exp_model, M(i,:), Bv, I(i,:), lb, ub, options);

    end

    
end


% return the diffusion tensor (fill with zeros in parts outside mask)
DT = zeros([prod(sz(1:end-1)), 7]);
DT(mask, :) = M;
DT = reshape(DT, [sz(1:end-1), 7]);
    

if nargout > 1 % if you want the FA, ADC, etc.

    % initialise variables
    FA_mask = zeros(size(M,1), 1);
    ADC_mask = zeros(size(M,1),1);
    VectorF_mask = zeros(size(M,1), 3);
    VectorF2_mask = zeros(size(M,1), 3);
    VectorF3_mask = zeros(size(M,1), 3);
    EigVals_mask = zeros(size(M,1), 3);
   
    parfor i = 1:size(M,1)
     
        Mi = M(i,:);

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
        ADCv = mean(EigenValues);%(EigenValues(1)+EigenValues(2)+EigenValues(3))/3;

        % Fractional Anistropy (2 different definitions exist)
        % First FA definition:
        %FAv=(1/sqrt(2))*( sqrt((EigenValues(1)-EigenValues(2)).^2+(EigenValues(2)-EigenValues(3)).^2+(EigenValues(1)-EigenValues(3)).^2)./sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2) );
        % Second FA definition:
        FA_mask(i)=sqrt(1.5)*( sqrt((EigenValues(1)-ADCv).^2+(EigenValues(2)-ADCv).^2+(EigenValues(3)-ADCv).^2)./sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2) );
        ADC_mask(i)=ADCv;
%         VectorF_mask(i,:)=EigenVectors(:,3)*EigenValues_old(3);
%         VectorF2_mask(i,:)=EigenVectors(:,2)*EigenValues_old(2);
%         VectorF3_mask(i,:)=EigenVectors(:,1)*EigenValues_old(1);
        VectorF_mask(i,:)=EigenVectors(:,3);
        VectorF2_mask(i,:)=EigenVectors(:,2);
        VectorF3_mask(i,:)=EigenVectors(:,1);

        EigVals_mask(i,:) = EigenValues_old([3 2 1]);

        
    end
    
    
    
    % reshape to match input dimensions, with zero filling outside the mask
    VectorF = zeros([prod(sz(1:end-1)), 3]);
    VectorF(mask(:), :) = VectorF_mask;
    VectorF = reshape(VectorF, [sz(1:end-1), 3]);
    VectorF2 = zeros([prod(sz(1:end-1)), 3]);
    VectorF2(mask(:), :) = VectorF2_mask;
    VectorF2 = reshape(VectorF2, [sz(1:end-1), 3]);
    VectorF3 = zeros([prod(sz(1:end-1)), 3]);
    VectorF3(mask(:), :) = VectorF3_mask;
    VectorF3 = reshape(VectorF3, [sz(1:end-1), 3]);
    VectorField = cat(ndims(VectorF)+1, VectorF, VectorF2, VectorF3);
    
    EigVals = zeros([prod(sz(1:end-1)), 3]);
    EigVals(mask(:), :) = EigVals_mask;
    EigVals = reshape(EigVals, [sz(1:end-1), 3]);
    
    % quick hack because reshape sz needs to be at least length 2
    if length(sz) == 2
        sz = [sz(1), 1, sz(2)];
    end
    
    FA = zeros([prod(sz(1:end-1)), 1]);
    FA(mask(:)) = FA_mask;
    FA = reshape(FA, sz(1:end-1));
    ADC = zeros([prod(sz(1:end-1)), 1]);
    ADC(mask(:)) = ADC_mask;
    ADC = reshape(ADC, sz(1:end-1));
    
    
    
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