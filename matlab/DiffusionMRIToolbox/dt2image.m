function [ img ] = dt2image( DT, b )
% DT2IMAGE Inverts the transform in fit_DT
% 
% The suggested usage of this file is to easily check how well a diffusion
% tensor fits to your data. Eg:
%
% DT = fit_DT( im, b );
% im_fitted = dt2image(DT, b);
% residuals = im_fitted - im;
%
% Ideally the residuals would be just noise - if there is any pattern to
% them, the model might not be a good fit.
%
%   DT is the diffusion tensor coefficients
%   b is the b values
%   img is the reconstructed image
%
% See also fit_DT

    
% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright � 2014 University of Oxford
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
narginchk(2, 2);
nargoutchk(0, 1);

sz = size(DT);

% reshape the tensor for efficient computing
DT = reshape(DT, [prod(sz(1:end-1)), sz(end)]);

% convert the S0 into log(S0) for linear fitting
DT(:,7) = log(DT(:,7));

% Sort all b matrices in to a vector Bv=[Bxx,2*Bxy,2*Bxz,Byy,2*Byz,Bzz];
Bv=squeeze([b(1,1,:),2*b(1,2,:),2*b(1,3,:),b(2,2,:),2*b(2,3,:),b(3,3,:)])';

% add the column of -1
Bv = [Bv, -ones(size(Bv,1), 1)];

% invert the model
logI = -Bv * DT';

% take the exponent (it was logged for the linear fitting)
img = exp(logI)';

sz_I = size(img);

% reshape back to the correct size
img = reshape(img, [sz(1:end-1), sz_I(end)]);

end

