function [ M_weighted ] = weighted_linear_fit( Y, X, func, warnings_on )
%WEIGHTED_LINEAR_FIT fits the function Y = MX (matrices) such that the
% residuals follow a gaussian distribution after having @func applied to
% them. I.e. func(Y) - func(MX) follows a gaussian with the same variance
% for each column of Y
%
% Inputs:
%   Y: [m * n] array, the signal you are trying to fit
%   X: [k * n] array, the x values for the signal. If you want a y
%       intercept, you can add a row of ones to X.
%   func: The function to get from Y to some other signal that you want a
%       least square fit in. For example, func = @(z) exp(z);
%   WARNINGS_ON is optional (default true), it tells you if something bad
%       is happening. 
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
%     M_weighted = weighted_linear_fit( Y, x, func, 1 );
%     M_weighted should equal -0.01.
% 
%
% The function assumes the same weighting matrix can be used for all of Y.
% If it can't, split up the data. This is similar to the inbuilt function
% robustfit, but handles arrays with multiple samples. Robustfun also has a
% limited number of weighting functions, rather than handling any possible
% function and trying to equalise the variance of the residuals.


% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright © 2015 University of Oxford
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
narginchk(4, 6);
nargoutchk(0, 6);

if nargin < 4
    warnings_on = 1;
end

% check if the variables are the right orientation
if size(Y,2) == size(X,2)
    % great, we have the right orientation
elseif size(Y,1) == size(X,2)
    if warnings_on, disp('Y has the wrong orientation'), end;
    Y = Y';
elseif size(Y,2) == size(X,1)
    if warnings_on, disp('X has the wrong orientation'), end;
    X = X';
elseif size(Y,1) == size(X,1)
    if warnings_on, disp('X and Y have the wrong orientation'), end;
    X = X';
    Y = Y';
else
    error('Incorrect variable sizes')
end
    
    
% fit the model
M = (pinv(X') * Y')';

% get the residuals
R = func(Y) - func(M*X);

% if we have only 1 sample, weight by the ratio of the root mean square
% error in the linear space to the real space
if size(Y,1) == 1
    V = sqrt((R.^2) ./ ((Y-M*X).^2 + eps));
else
    % get the variance - we want the variance of the residuals to be the 
    % same for all points on the curve
    V = var(R);
end

% define the weighting function
W = diag(V);

% fit the weighted model
M_weighted = (pinv(X * (W.^2) * X') * X * (W.^2) * Y')';

if warnings_on
    R_weighted = func(Y) - func(M_weighted * X);

    if sqrt(sum(R(:).^2)) < sqrt(sum(R_weighted(:).^2))
        disp(['Be careful, the weighting step increased your sum of squares from ' ...
            num2str(sqrt(sum(R(:).^2))) ' to ' num2str(sqrt(sum(R_weighted(:).^2)))])
        
    end
end

end