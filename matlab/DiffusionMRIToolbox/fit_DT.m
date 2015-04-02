function [ DT, FA, ADC, VectorField, EigVals] = fit_DT( im, b, thresh_val, method, weighting)
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
%       'nonlinear' fits the tensor in the signal domain (recommended for
%       high b-values or when SNR is low)
%
%   WEIGHTING can be used to do a weighted linear tensor fit. It must be
%   a vector of length N. This can help to minimise artefacts from linear
%   fitting. (optional, not needed if SNR is good)
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
%	in the last dimension
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
% Version: 0.1.8
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
narginchk(2, 5);
nargoutchk(0, 5);

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

% take the log of the image to linearise the equation
imlog = log(im);

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

% reshape Ilog so we can use the efficient \ operator
imlog = reshape(imlog, [prod(sz(1:end-1)), sz(end)]);

% Fit the tensor model (linear, unconstrained)
M = (Bv \ -imlog')';
% IDX = kmeans(imlog, 3);
% M = zeros(size(imlog,1), 7);
% for i = 1:3
%     disp(num2str(i))
%     M(IDX == i, :) = weighted_linear_fit(-imlog(IDX == i,:), Bv', @(z)exp(-z));
% end
%M = weighted_linear_fit(-imlog, Bv', @(z)exp(-z));

% This is the same as (pinv(Bv) * (-imlog'))';
% also same as (pinv(Bv' * Bv) * Bv' * (-imlog'))';
% weighting (NxN) goes   ^ here      ^ and here
if perform_weighting

    W = eye(length(Bv));
    W(W == 1) = weighting(:);
    % weighted fit

    M = (pinv(Bv' * (W.^2) * Bv) * Bv' * (W.^2) * (-imlog'))';
end

% convert log(S0) to S0
M(:,7) = exp(M(:,7));
M(isnan(M)) = 0;
M(isinf(M)) = 0;

% perform nonlinear fitting (if required)
if strcmp(method, 'nonlinear')

    I_vector = reshape(im, [prod(sz(1:end-1)), sz(end)]);
    
    if isscalar(thresh_val)
        thresh_val = I_vector(:,1) > thresh_val;
    else
        thresh_val = logical(thresh_val);
    end
    
    % declare variables for nonlinear fitting
    I_nlfit = double(I_vector(thresh_val(:), :));
    M_nl = double(M(thresh_val(:), :));
    Bv = Bv(:,1:6)';

    options = optimoptions('lsqcurvefit','Jacobian','on', 'DerivativeCheck', 'off', ...
        'display', 'off', 'TypicalX', mean(M_nl));
    % M = [xx,xy,xz,yy,yz,zz,S0]
    lb = [0 -3E-3 -3E-3 0 -3E-3 0 0]; % cross terms are allowed to be -ve
    ub = [zeros(1,6)+3E-3, inf];

    
    if (size(I_nlfit,1) > 1000) && (matlabpool('size') == 0)
        disp('Open matlabpool for parallel processing')
    end
    
    parfor i = 1:size(I_nlfit,1)

        M_nl(i,:) = lsqcurvefit(@mono_exp_model, M_nl(i,:), Bv, I_nlfit(i,:), lb, ub, options);

    end


    M(thresh_val(:),:) = M_nl;
    
end
    
    
% return the diffusion tensor
DT = reshape(M, [sz(1:end-1), 7]);



if nargout > 1 % if you want the FA, ADC, etc.

    % initialise variables
    FA = zeros(size(M,1), 1);
    ADC = zeros(size(M,1),1);
    VectorF = zeros(size(M,1), 3);
    VectorF2 = zeros(size(M,1), 3);
    VectorF3 = zeros(size(M,1), 3);
    EigVals = zeros(size(M,1), 3);
   
    parfor i = 1:size(M,1)
    
        if isscalar(thresh_val) % if thresh_val is a scalar, use it as a threshold
            if im(i) < thresh_val
                continue
            end
        else
            if thresh_val(i) == 0 % otherwise use it as a mask
                continue
            end
        end
        
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
    
    % quick hack because reshape sz needs to be at least length 2
    if length(sz) == 2
        sz = [sz(1), 1, sz(2)];
    end
    
    % reshape to match input dimensions
    FA = reshape(FA, sz(1:end-1));
    ADC = reshape(ADC, sz(1:end-1));
    VectorF = reshape(VectorF, [sz(1:end-1), 3]);
    VectorF2 = reshape(VectorF2, [sz(1:end-1), 3]);
    VectorF3 = reshape(VectorF3, [sz(1:end-1), 3]);
    
    VectorField = cat(ndims(FA)+1, VectorF, VectorF2, VectorF3);
    
    EigVals = reshape(EigVals, [sz(1:end-1), 3]);
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