function [ DKI ] = fit_DKI_model( im, b, thresh_val, method, weighting)
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
% An alternative method is implemented in fit_DT_YHd_method.
%
% See also dt2image.

    
% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright © 2014-2015 University of Oxford
% Version: 0.1.8
% $Rev: 1712 $
% $Date: 2015-04-02 11:33:57 +0100 (Thu, 02 Apr 2015) $
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


%%





% this script generates the bi-exponential parameters:
close all
clear all
clc

%cd /xraid2/irvin/vnmrsys/data/studies/s_20150312_01/data/fse3d_cs_06.fid/
cd /xraid2/irvin/vnmrsys/data/studies/s_20150306_01/data/fse3d_cs_05.fid/

I = load_untouch_nii('image_mag.nii');
I = I.img;

[I,B_mat,bvalue]=DWI_setup(I,'/xraid2/irvin/vnmrsys/data/studies/s_20150306_01/data/fse3d_cs_05.fid/');


B_mat = B_mat ./ repmat(permute(bvalue, [3,1,2]), [3 3 1]);

mask = I(:,:,:,1) > 1E4;

N = size(B_mat,3);
% define A matrix

A_D = zeros(N, 6);
A_K = zeros(N, 15);


A_D = squeeze([B_mat(1,1,:), B_mat(2,2,:), B_mat(3,3,:), 2*B_mat(1,2,:), 2*B_mat(1,3,:), 2*B_mat(2,3,:)])';
% A_D in the paper is not scaled by the b value

A_K = squeeze([B_mat(1,1,:).^2, B_mat(2,2,:).^2, B_mat(3,3,:).^2, ...
    4*B_mat(1,1,:).*B_mat(1,2,:), 4*B_mat(1,1,:).*B_mat(1,3,:), 4*B_mat(1,2,:).*B_mat(2,2,:), ...
    4*B_mat(1,1,:).*B_mat(1,3,:), 4*B_mat(2,2,:).*B_mat(2,3,:), 4*B_mat(2,3,:).*B_mat(3,3,:), ...
    2*B_mat(1,1,:).*B_mat(2,2,:) + 4*B_mat(1,2,:).*B_mat(1,2,:), ...
    2*B_mat(1,1,:).*B_mat(3,3,:) + 4*B_mat(1,3,:).*B_mat(1,3,:), ...
    2*B_mat(2,2,:).*B_mat(3,3,:) + 4*B_mat(2,3,:).*B_mat(2,3,:), ...
    8*B_mat(1,2,:).*B_mat(1,3,:) + 4*B_mat(1,1,:).*B_mat(2,3,:), ...
    8*B_mat(1,2,:).*B_mat(2,3,:) + 4*B_mat(2,2,:).*B_mat(1,3,:), ...
    8*B_mat(1,3,:).*B_mat(2,3,:) + 4*B_mat(3,3,:).*B_mat(1,2,:)])';
% this takes care of cross terms by saying that (for example) 6 x y^2
% = 2 (x^2) (y^2) + 4 (xy) (xy)
% so it uses these elements of the B matrix.




sz = size(I);
B = reshape(I, [prod(sz(1:3)), sz(4)]);
B = bsxfun(@rdivide, B, B(:,1));

B = log(B);

A = [-bsxfun(@times, bvalue', A_D), 1/6 * bsxfun(@times, bvalue'.^2, A_K)];

X = (pinv(A) * B');

C = [-A_D,                              zeros(size(A_D,1), size(A_K,2)); ...
    zeros(size(A_K,1), size(A_D,2)),    -A_K; ...
    3/max(bvalue(:))*A_D,               A_K];

d = zeros(size(C,1), size(X,1));


lb = [0 0 0 0 -1 -1 -1 ... % [S0 DTI
    0 0 0 -1 -1 -1 -1 -1 -1 0 0 0 -1 -1 -1] * 3E-3; % DKI ]
ub = [inf, ones(1,21) * 3E-3]; % [S0 DTI DKI]

x_lin = X(:,i);
Y = reshape(I, [prod(sz(1:3)), sz(4)]);
y = Y(i,:);
i = sub2ind(sz(1:3), 88, 57, 36);
x_lin = [y(1); x_lin];

x_nonlin = lsqcurvefit(@DKI_model, double(x_lin)', -A', double(y), lb, ub, options);


% when fitting statistical model at wolfram alpha:
% integrate (exp(-(((D-A)^2)/(2*S^2) + ((D-B)^2)/(2*T^2) + ((D-C)^2)/(2*U^2) ))) dD
