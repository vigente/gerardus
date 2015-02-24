function [I_interp, Tx, Ty, Tz, Txinv, Tyinv, Tzinv] = regdemons3d(I_rigid, I_mov, Tx, Ty, Tz, Txinv, Tyinv, Tzinv, iteration, range, sigma)

% REGDEMONS Image registration using Thirion's Demons algorithm.
%
% This function implements the Demons algorithm in
% Thirion, J.-P. (1998). Image matching as a diffusion process: an analogy
% with Maxwell's demons. Medical Image Analysis, 2(3), 243–260.
%
% [I_INTERP, TX, TY, TZ, TXINV, TYINV, TZINV] = REGDEMONS(I_RIGID, I_MOV, TX0, TY0, TZ0, TXINV0, TYINV0, TZINV0, ITERATION, HSIZE, SIGMA)
%
%   I_RIGID and I_MOV are matrices with the fixed and moving images,
%   respectively. They must be of the same size.
%
%   TX0, TY0 and TZ0 are three matrices with the same size as the images. They
%   contain the initial translation of each pixel, in the X, Y and Z
%   directions, respectively.

%   TXINV0, TYINV0 and TZINV0 are three matrices with the same size as the images. They
%   contain the initial inverse translation of each pixel, in the X, Y and Z
%   directions, respectively.
%
%   ITERATION is a scalar with the number of iterations for the algorithm.
%   It depends on the images size, but we usually work with values of up to
%   200.
%
%   HSIZE, SIGMA are the parameters of the Gaussian kernel. HSIZE is a
%   3-vector with the size of the kernel window. SIGMA is a scalar with the
%   standard deviation of the Gaussian.
%
%   I_INTERP is the registered image (i.e. I_MOV after applying the Demons
%   algorithm).
%
%   TX, TY, TY are the transform computed by the Demons algorithm.
%   TXINV, TYINV, TYINV are the inverse transform computed by the Demons algorithm.

% Author: Adam Szmul <aszmul@gmail.com>
% Copyright © 2015 University of Oxford
% Version: 0.1.1
% $Rev: 1644 $
% $Date: 2014-12-11 17:36:00 +0000 (Thu, 11 Dec 2014) $
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


narginchk(11, 11);
nargoutchk(0, 7);




[n m k ] = size(I_mov);
[x, y ,z]=meshgrid(1:m,1:n, 1:k);
I_interp = mirt3D_mexinterp(I_mov, x+Tx, y+Ty, z+Tz); % faster linear interpolation by mex file
        
  
[Gx Gy Gz] = gradient(I_rigid);

    for i=1:iteration
        
        
            
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (m-s)    
       Diff = (I_interp - I_rigid) ;  
     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %% Eq 4.
       
        Vy = -(Diff.* (Gy))./((Gx.*Gx + Gy.*Gy + Gz.*Gz) + Diff.*Diff + 0.0001);  
        Vx = -(Diff.* (Gx))./((Gx.*Gx + Gy.*Gy + Gz.*Gz) + Diff.*Diff + 0.0001);       
        Vz = -(Diff.* (Gz))./((Gx.*Gx + Gy.*Gy + Gz.*Gz) + Diff.*Diff + 0.0001);  
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Insterting zeros instead of Nan when divided by 0

   Vx(isnan(Vx))=0;            % to eliminate NaN instead of eps
   Vy(isnan(Vy))=0;            % to eliminate NaN instead of eps
   Vz(isnan(Vz))=0;            % to eliminate NaN instead of eps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        Tx = Tx + Vx;          
        Ty = Ty + Vy;          
        Tz = Tz + Vz;

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Inverse transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Txinv = Txinv - Vx;          
        Tyinv = Tyinv - Vy;          
        Tzinv = Tzinv - Vz;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%% Faster Gaussian smoothing by mex file

        Tx=imgaussian(Tx, sigma, range);   
        Ty=imgaussian(Ty, sigma, range);
        Tz=imgaussian(Tz, sigma, range);

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Inverse transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        Txinv=imgaussian(Txinv, sigma, range);   
        Tyinv=imgaussian(Tyinv, sigma, range);
        Tzinv=imgaussian(Tzinv, sigma, range);
 
%% Insterting zeros on borders to eliminate border problems       
        
Tx_zero = zeros(size(Tx));
Tx_zero(2:size(Tx,1)-1,2:size(Tx,2)-1,2:size(Tx,3)-1) = Tx(2:size(Tx,1)-1,2:size(Tx,2)-1,2:size(Tx,3)-1);
Tx=Tx_zero;

Ty_zero = zeros(size(Ty));
Ty_zero(2:size(Ty,1)-1,2:size(Ty,2)-1,2:size(Ty,3)-1) = Ty(2:size(Ty,1)-1,2:size(Ty,2)-1,2:size(Ty,3)-1);
Ty=Ty_zero;

Tz_zero = zeros(size(Tz));
Tz_zero(2:size(Tz,1)-1,2:size(Tz,2)-1,2:size(Tz,3)-1) = Tz(2:size(Tz,1)-1,2:size(Tz,2)-1,2:size(Tz,3)-1);
Tz=Tz_zero;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Inverse transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tx_zero = zeros(size(Txinv));
Tx_zero(2:size(Txinv,1)-1,2:size(Txinv,2)-1,2:size(Txinv,3)-1) = Txinv(2:size(Txinv,1)-1,2:size(Txinv,2)-1,2:size(Txinv,3)-1);
Txinv=Tx_zero;

Ty_zero = zeros(size(Tyinv));
Ty_zero(2:size(Tyinv,1)-1,2:size(Tyinv,2)-1,2:size(Tyinv,3)-1) = Tyinv(2:size(Tyinv,1)-1,2:size(Tyinv,2)-1,2:size(Tyinv,3)-1);
Tyinv=Ty_zero;

Tz_zero = zeros(size(Tzinv));
Tz_zero(2:size(Tzinv,1)-1,2:size(Tzinv,2)-1,2:size(Tzinv,3)-1) = Tzinv(2:size(Tzinv,1)-1,2:size(Tzinv,2)-1,2:size(Tzinv,3)-1);
Tzinv=Tz_zero;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Standard built in function       

[x, y ,z] = meshgrid(1:size(I_mov,2),1:size(I_mov,1), 1:size(I_mov,3));
I_interp = mirt3D_mexinterp(I_mov, x+Tx, y+Ty, z+Tz);  % faster linear interpolation by mex file



%       W = smooth3(I_mov);
%       I_interp(isnan(I_interp))=W(isnan(I_interp));  % replace NaN with mean near value
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

    end


end
