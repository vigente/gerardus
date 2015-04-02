function [ TV, TV_grad ] = forward_TV( I )
%FORWARD_TV Total variation of a 3D image
%   Returns the sum of finite differences in 3D, and the residuals
%
%  [TV, TV_GRAD] = forward_TV(I)
% 
%   I is a 3D real image
%
%   TV is the total variation of I. That is, the sum of the magnitude of 
%   finite differences in x, y, and z
%   
%   TV_GRAD can be used to find the residuals between I and a smoothed
%   version of I. It contains Dx, Dy and Dz, concatenated in the 4th
%   dimension.


% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.1.2
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
narginchk(1,1);
nargoutchk(0, 2);


% Let mex do the work for you
if exist('forward_TV_aux', 'file')
    [Dx, Dy, Dz, TV] = forward_TV_aux(double(I), 1);  
else
    
    disp('Please compile the forward total variation transform for increased speed.')
    disp('Just find the file forward_TV_aux.cpp and type:')
    disp('mex forward_TV_aux.cpp')
    
    
    % This is the equivalent to the above mex function in Matlab, but the
    % mex file is faster

    Dx = I([2:end,end],:,:) - I;
    Dy = I(:,[2:end,end],:) - I;
    Dz = I(:,:,[2:end,end]) - I;
    

    TV = sum(abs(Dx(:)) + abs(Dy(:)) + abs(Dz(:)));

end

% concatenate the derivative
TV_grad = cat(4,Dx,Dy,Dz);


end

