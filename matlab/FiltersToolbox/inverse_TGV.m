function [ res ] = inverse_TGV( TGV )
%INVERSE_TGV Returns the residuals of 2nd order total generalised variation
%
% TGV comes from forward_TGV (see example in that function)
% RES are the residuals, subtract them from the image to reduce the TGV
%
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

% check arguments
narginchk(1,1);
nargoutchk(0, 1);

N = ndims(TGV) - 1;

if N == 2
    
    Dxx = TGV(:,:,1);
    Dxy = TGV(:,:,2);
    Dyy = TGV(:,:,3);
    
    Dx = ((Dxx([end, 1:end-1], :, :) - Dxx) + ... % -ve x backwards difference of Dxx
          (Dxy(:, [end, 1:end-1], :) - Dxy)) / 2; % -ve y backwards difference of Dxy
    
    Dy = ((Dyy(:, [end, 1:end-1], :) - Dyy) + ...
          (Dxy([end, 1:end-1], :, :) - Dxy)) / 2;
   
    res = -((Dx([end, 1:end-1], :, :) - Dx) + ...
           (Dy(:, [end, 1:end-1], :) - Dy)) / 2;

elseif N == 3

    Dxx = TGV(:,:,1);
    Dxy = TGV(:,:,2);
    Dxz = TGV(:,:,3);
    Dyy = TGV(:,:,4);
    Dyz = TGV(:,:,5);
    Dzz = TGV(:,:,6);
    
    Dx = ((Dxx([end, 1:end-1], :, :) - Dxx) + ... % -ve x backwards difference of Dxx
          (Dxy(:, [end, 1:end-1], :) - Dxy) + ... % -ve y backwards difference of Dxy
          (Dxz(:, :, [end, 1:end-1]) - Dxz)) / 3; % -ve z backwards difference of Dxz
    
    Dy = ((Dxy([end, 1:end-1], :, :) - Dxy) + ... % -ve x backwards difference
          (Dyy(:, [end, 1:end-1], :) - Dyy) + ... % -ve y backwards difference
          (Dyz(:, :, [end, 1:end-1]) - Dyz)) / 3; % -ve z backwards difference
      
    Dz = ((Dxz([end, 1:end-1], :, :) - Dxz) + ... % -ve x backwards difference
          (Dyz(:, [end, 1:end-1], :) - Dyz) + ... % -ve y backwards difference
          (Dzz(:, :, [end, 1:end-1]) - Dzz)) / 3; % -ve z backwards difference
      
    res = -((Dx([end, 1:end-1], :, :) - Dx) + ...
           (Dy(:, [end, 1:end-1], :) - Dy) + ...
           (Dz(:, :, [end, 1:end-1]) - Dz)) / 3;
    
else
    error(['I can only handle 2D or 3D images, and the last dimension is ' ...
    'concatenated dxx, dxy, etc'])
end
