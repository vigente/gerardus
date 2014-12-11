function [I_interp, Tx, Ty] = regdemons(I_rigid, I_mov,Tx,Ty,iteration, range, sigma)
% REGDEMONS Image registration using Thirion's Demons algorithm.
%
% This function implements the Demons algorithm in
%
% Thirion, J.-P. (1998). Image matching as a diffusion process: an analogy
% with Maxwell's demons. Medical Image Analysis, 2(3), 243–260.
%
% [I_INTERP, TX, TY] = REGDEMONS(I_RIGID, I_MOV, TX0, TY0, ITERATION, HSIZE, SIGMA)
%
%   I_RIGID and I_MOV are matrices with the fixed and moving images,
%   respectively. They must be of the same size.
%
%   TX0 and TY0 are two matrices with the same size as the images. They
%   contain the initial translation of each pixel, in the X and Y
%   directions, respectively.
%
%   ITERATION is a scalar with the number of iterations for the algorithm.
%   It depends on the images size, but we usually work with values of up to
%   1,000.
%
%   HSIZE, SIGMA are the parameters of the Gaussian kernel. HSIZE is a
%   2-vector with the size of the kernel window. SIGMA is a scalar with the
%   standard deviation of the Gaussian.
%
%   I_INTERP is the registered image (i.e. I_MOV after applying the Demons
%   algorithm).
%
%   TX, TY are the transform computed by the Demons algorithm.

% Author: Adam Szmul <aszmul@gmail.com>
% Copyright © 2014 University of Oxford
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
narginchk(7, 7);
nargoutchk(0, 3);

  I_interp = I_mov;    
        
[Gx Gy] = gradient(I_rigid);

    for i=1:iteration
            
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (m-s)    
       Diff = (I_interp - I_rigid) ;  
     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %% Eq 4.
        Vy = -(Diff.* (Gx))./((Gx.^2+Gy.^2) + Diff.^2 + 0.0001);   % changed order during iteration of X and Y
 
        Vx = -(Diff.* (Gy))./((Gy.^2 +Gx.^2) + Diff.^2 + 0.0001);  % changed order during iteration of X and Y    
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Insterting zeros instead of Nan when divided by 0

   Vx(isnan(Vx))=0;            % to eliminate NaN instead of eps
   Vy(isnan(Vy))=0;            % to eliminate NaN instead of eps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Smoothing transformation field
        Hsmooth=fspecial('gaussian',range,sigma);
        Vx=imfilter(Vx,Hsmooth);
        Vy=imfilter(Vy,Hsmooth);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Tx = Tx + Vx;           % changed order during iteration of X and Y
        Ty = Ty + Vy;           % changed order during iteration of X and Y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Standard built in function       

[x,y]=ndgrid(1:size(I_rigid,1),1:size(I_rigid,2));
I_interp = interp2(I_mov, y+Ty,x+Tx, 'nearest');  % the other order of X and Y

     I_interp(isnan(I_interp))=mean(mean(I_interp(~isnan(I_interp))));  % eliminates NaN

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

    end


end
