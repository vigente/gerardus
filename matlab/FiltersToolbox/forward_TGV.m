function [ TGV, TGV_grad ] = forward_TGV( I )
%FORWARD_TGV Second order total generalised variation of an image (up to 3D)
%   
%   I is an image (2D or 3D)
%     
%   TGV is the second order total generalised variation of I
%   TGV_grad is the derivative (feed this into inverse_TGV)


% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright © 2015 University of Oxford
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

% Example usage:
% [X, Y] = ndgrid(1:100, 1:100);
% I = X.^2 + 3*Y.^2;
% [ TGV, TGV_grad ] = forward_TGV( I );
% res = inverse_TGV(TGV_grad);
% Note that the residuals are all zero (ignoring edge effects)
% this is because the second derivative of a quadratic is a constant
% so you don't need to subtract anything to make it sparser

% check arguments
narginchk(1,1);
nargoutchk(0, 2);

N = ndims(I);

Dx = I([2:end,1],:,:) - I;
Dy = I(:,[2:end,1],:) - I;

Dxx = Dx - Dx([2:end,1],:,:);
Dyy = Dy - Dy(:,[2:end,1],:);

% x derivative of Dy is the same as y derivative of Dx
Dxy = Dy - Dy([2:end,1],:,:); 

if N > 2
    
    Dz = I(:,:,[2:end,1]) - I;
    
    Dxz = Dz - Dz([2:end,1],:,:); % x derivative of Dz

    Dyz = Dz - Dz(:,[2:end,1],:); % y derivative of Dz

    Dzz = Dz - Dz(:,:,[2:end,1]); % z derivative of Dz
end

if N < 3
    % diagonals count double
    TGV = sum(abs(Dxx(:)) + abs(Dyy(:)) + 2 * abs(Dxy(:)));
    
    TGV_grad = cat(3,Dxx,Dxy,Dyy);
    
else
    
    TGV = sum(abs(Dxx(:)) + abs(Dyy(:)) + 2 * abs(Dxy(:)) + ... 
          2 * abs(Dxz(:)) + 2 * abs(Dyz(:)) + abs(Dzz(:)));

    TGV_grad = cat(N+1,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
end


end

