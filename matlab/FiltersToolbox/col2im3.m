function [threeD] = col2im3(reshaped3D2D,blocksX,blocksY,blocksZ)
%COL2IM3 rearranges matrix columns into a 3D image

% 22 July 2015
%
% col2im3 the same as col2im, only for 3D arrays
% It reshapes each of the columns of the 2D array into cubes, then to stacks them back together layer by layer
% The inputs are a 2D array and numbers of blocks in each direction
% outputs from im2col3
% The output is a 3D array. The dimensions of which are
% blocksize*blocksX by blocksize*blocksY by blocksize*blocksZ
%
% Author: Nicolas Basty <nicolas.basty@eng.ox.ac.uk>
% Copyright ??? 2015 University of Oxford
% Version: 0.1.0
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
% _________________________________________________________________

blocksize = nthroot(size(reshaped3D2D,1),3); %cubic root gives size of a side of the cube

% blocksZ, number of slabs (blocks in Z) from output from previous function
% blocksX, number of blocks in x
% blocksZ, number of blocks in y

if  rem(blocksize,1) ~= 0 %checks if blocksize is a cubic number
    fprintf('Input array is supposed to have a cubic number as first dimension! \n')
    return
end

if blocksX*blocksY*blocksZ ~= size(reshaped3D2D,2)
    fprintf('The number of blocks does not match the input array. \n')
    return
    %checks if blocksize makes sense
end

blocknum = 1; %just initialising the block number

threeD = zeros(blocksize*blocksX,blocksize*blocksY,blocksize*blocksZ);
%initialising the final array

for slabn = 1:blocksZ %Z
    %loop goes in the Z direction slab by slab
    for xnum = 1:blocksX
        for ynum = 1:blocksY
            
            temp = reshaped3D2D(:,blocknum);%temporarily stores the n-th block in Ax1 form
            temp2 = reshape(temp,blocksize,blocksize,blocksize);%this reshapes the Ax1 into axaxa
            
            threeD(xnum+(blocksize*(xnum-1)-xnum+1):(xnum*blocksize),ynum+(blocksize*(ynum-1)-ynum+1):(ynum*blocksize),slabn+(blocksize*(slabn-1)-slabn+1):(slabn*blocksize)) = temp2;
            %stores the cube at the location within the big 3D array
            
            blocknum = blocknum +1; %move on to next block
        end
    end
end

end