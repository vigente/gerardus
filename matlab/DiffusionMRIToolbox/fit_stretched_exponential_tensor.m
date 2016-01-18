function [ DT, I_fit] = fit_stretched_exponential_tensor( im, b, thresh_val)

    
% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright © 2014-2015 University of Oxford
% Version: 0.1.8
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

perform_weighting = false;
if nargin > 4
    perform_weighting = true;
    if length(weighting) ~= sz(end)
        disp('Weighting vector is wrong length')
        perform_weighting = false;
    end
end


% handle a vector with the wrong orientation
if (length(sz) == 2) && (sz(2) == 1)
    im = im';
    sz = size(im);
end

% take the log of the image to linearise the equation
imlog = log(abs(im));

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
%M = (Bv \ -imlog')';
M = weighted_linear_fit(-imlog, Bv', @(z)exp(-z));



% convert log(S0) to S0
M(:,7) = exp(M(:,7));

% add a column of ones
M = [M, ones(size(M,1),1)];

% get rid of any bad rows
M(isnan(M)) = 0;
M(isinf(M)) = 0;

% perform nonlinear fitting 
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

options = optimoptions('lsqcurvefit','Jacobian','off', 'DerivativeCheck', 'off', ...
    'display', 'off', 'TypicalX', mean(M_nl));
% M = [xx,xy,xz,yy,yz,zz,S0]
lb = [0 -3E-3 -3E-3 0 -3E-3 0 0 0]; % cross terms are allowed to be -ve
ub = [zeros(1,6)+3E-3, inf 1];


if (size(I_nlfit,1) > 1000) && (matlabpool('size') == 0)
    disp('Open matlabpool for parallel processing')
end

parfor i = 1:size(I_nlfit,1)

    M_nl(i,:) = lsqcurvefit(@stretch_exp_model, M_nl(i,:), Bv, I_nlfit(i,:), lb, ub, options);

end


M(thresh_val(:),:) = M_nl;

I_fit = bsxfun(@times, M(:,7), exp(-(bsxfun(@power, M(:,1:6) * Bv, M(:,8)))));

I_fit = reshape(I_fit, sz);

    
% return the diffusion tensor
DT = reshape(M, [sz(1:end-1), 8]);




% the function to fit inside the parfor loop
function [F, J] = stretch_exp_model(x, Bv)

F = x(7) * exp(-((x(1:6) * Bv)).^x(8));
    
% if nargout > 1 % Jacobian
%     J = zeros(length(F), length(x));
%     
%     J(:,1) =  x(7) * -Bv(1,:) .* exp(-(x(1:6) * Bv));
%     J(:,2) =  x(7) * -Bv(2,:) .* exp(-(x(1:6) * Bv));
%     J(:,3) =  x(7) * -Bv(3,:) .* exp(-(x(1:6) * Bv));
%     J(:,4) =  x(7) * -Bv(4,:) .* exp(-(x(1:6) * Bv));
%     J(:,5) =  x(7) * -Bv(5,:) .* exp(-(x(1:6) * Bv));
%     J(:,6) =  x(7) * -Bv(6,:) .* exp(-(x(1:6) * Bv));        
%     J(:,7) = exp(-(x(1:6) * Bv));
% end