function [ M_weighted] = weighted_linear_fit( Y, X, func, param )
%WEIGHTED_LINEAR_FIT fits the function Y = MX (matrices) such that the
% residuals are minimised after having @func applied to them. 
% I.e. func(Y) - func(MX) is the problem we are really trying to
% solve.
%
% Inputs:
%   Y: [m * n] array, the signal you are trying to fit (linearised)
%   X: [k * n] array, the x values for the signal. If you want a y
%       intercept, you can add a row of ones to X.
%   func: The function to transform Y back into the signal that you really 
%       want to fit. See example below.
%   param: A structure containing:
%       verbose (logical - default false)
%       unique_weights (logical - default true) computes a different 
%           weighting matrix for each signal (use when there is a diverse 
%           range of coefficients)
%       rician (logical - default false) defines the weighting matrix as
%           X.^2 (param.unique must be true). See "Weighted linear least
%           squares estimation of diffusion MRI parameters: strengths, 
%           limitations, and pitfalls"
%
% Output:
%   M_WEIGHTED: [m * k] array, containing the model coefficients.
%   
% Example usage: 
%     x = 0:10:1000;
%     y = abs(exp(-0.01 * x) + randn(size(x))/100);
%     % we have y = exp(m * x) where m = -0.01, plus some noise
%     % we can linearise this as log(y) = m*x
%     % using least squares, this is solved via m = x' \ Y'
%     % but this is susceptible to noise where signal is low
%     % instead, try a weighted fit
%     func = @(z) exp(z);
%     Y = log(y);
%     param.verbose = 1; param.unique_weights = 1; param.rician = 1;
%     M_weighted = weighted_linear_fit( Y, x, func, param);
%     M_weighted should equal -0.01.
% 
%
% The function assumes the same weighting matrix can be used for all of Y.
% If it can't, split up the data. This is similar to the inbuilt function
% robustfit, but handles arrays with multiple samples. Robustfun also has a
% limited number of weighting functions, rather than handling any possible
% function and trying to equalise the variance of the residuals.


% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright ? 2015 University of Oxford
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
narginchk(3, 4);
nargoutchk(0, 2);



if nargin < 4
    verbose = 0;
    unique_weights = 1;
    rician = 0;
else
    if isfield(param, 'verbose')
        verbose = param.verbose;
    else
        verbose = 0;
    end
    if isfield(param, 'unique_weights')
        unique_weights = param.unique_weights;
    else
        unique_weights = 1;
    end
    if isfield(param, 'rician')
        rician = param.rician;
    else
        rician = 0;
    end
end
    

% check if the variables are the right orientation
if size(Y,2) == size(X,2)
    % great, we have the right orientation
elseif size(Y,1) == size(X,2)
    if verbose, disp('Y has the wrong orientation'), end;
    Y = Y.';
elseif size(Y,2) == size(X,1)
    if verbose, disp('X has the wrong orientation'), end;
    X = X.';
elseif size(Y,1) == size(X,1)
    if verbose, disp('X and Y have the wrong orientation'), end;
    X = X.';
    Y = Y.';
else
    error('Incorrect variable sizes')
end
    
% turn off warnings (the / operator can generate them)
warn_state = warning;
warning off;
    
% fit the model without weighting
M = Y / X;

% get the residuals
MX = M*X;
R = Y - MX;
fY = func(Y);
Rf = fY - func(MX);

if unique_weights || (size(Y,1) == 1)
    
    v = ver('MATLAB');
    
    if str2double(v.Version) > 8.1 % R2015
        if ~isempty(gcp('nocreate'))
            pctRunOnAll warning('off')
        end
    else % R2013
        if matlabpool('size') ~= 0
            pctRunOnAll warning('off')
        end
    end
    
    M_weighted = zeros(size(M));
    parfor i = 1:size(Y, 1)
        
        if rician
            V = fY(i,:).^2; % or possibly not ^2
        else
            % weight by the ratio of the root mean square error in the linear space to the real space
            V = sqrt((Rf(i,:).^2) ./ (R(i,:).^2 + eps));
        end
        W = diag(V);
        M_weighted(i,:) = ((X * (W.^2) * X') \ (X * (W.^2) * Y(i,:)'))';
    end

else

    % just an average over the repetitions
    V = mean(sqrt((Rf.^2) ./ (R.^2 + eps)), 1);
    
    % define the weighting function
    W = diag(V);
    % fit the weighted model
    M_weighted = ((X * (W.^2) * X') \ (X * (W.^2) * Y'))';

end

% return warning state to original condition
warning(warn_state)


R_weighted = fY - func(M_weighted * X);

if verbose
    if sqrt(sum(Rf(:).^2)) < sqrt(sum(R_weighted(:).^2))
        disp(['Be careful, the weighting step increased your sum of squares from ' ...
            num2str(sqrt(sum(Rf(:).^2))) ' to ' num2str(sqrt(sum(R_weighted(:).^2)))])
        
    end
end

RMSE_unweighted = sum(Rf.^2, 2);
RMSE_weighted = sum(R_weighted.^2,2);
M_weighted(RMSE_unweighted < RMSE_weighted,:) = M(RMSE_unweighted < RMSE_weighted,:);
M_weighted(isnan(RMSE_weighted), :) = M(isnan(RMSE_weighted),:);

end