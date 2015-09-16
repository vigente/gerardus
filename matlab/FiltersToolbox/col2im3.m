function [threeD] = col2im3_edit2(col,mat,block,kind)
%COL2IM3 rearranges matrix columns into a 3D image

% 16 September 2015
%
% col2im3 the same as col2im, only for 3D arrays
% It reshapes each of the columns of the 2D array into cubes, 
% and stacks them back together
% The inputs are a 2D array, the size of the output and the image blocks,
% and the kind of window (sliding/distinct).
% output from im2col3 is a 3D array.
%
% Author: Nicolas Basty <nicolas.basty@eng.ox.ac.uk>
% Copyright ??? 2015 University of Oxford
% Version: 0.2
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

if ~ischar(kind),
    error(message('images:col2im:wrongBlockType'));
end

% size of output 3D matrix
blocksize_x = block(1);
blocksize_y = block(2);
blocksize_z = block(3);

Number_of_blocks_X = mat(1)/block(1); %how many blocks there are along each dimension
Number_of_blocks_Y = mat(2)/block(2);
Number_of_blocks_Z = mat(3)/block(3);

blocknum = 1; % initialising the block number

threeD = zeros(blocksize_x*Number_of_blocks_X,blocksize_y*Number_of_blocks_Y,blocksize_z*Number_of_blocks_Z);
%initialising the final array

if strcmp(kind, 'distinct') % Distinct
    
for slabn = 1:Number_of_blocks_Z
    for ynum = 1:Number_of_blocks_X
        for xnum = 1:Number_of_blocks_Y
            
            temp = col(:,blocknum);%temporarily stores the n-th block in Ax1 form
            temp2 = reshape(temp,blocksize_x,blocksize_y,blocksize_z);%this reshapes the Ax1 into axaxa
            threeD(xnum+(blocksize_x*(xnum-1)-xnum+1):(xnum*blocksize_x),ynum+(blocksize_y*(ynum-1)-ynum+1):(ynum*blocksize_y),slabn+(blocksize_z*(slabn-1)-slabn+1):(slabn*blocksize_z)) = temp2;
            %stores the cube at the location within the big 3D array
            blocknum = blocknum +1; %move on to next block
        
        end
    end
end

elseif strcmp(kind, 'sliding') % sliding
    threeD = reshape(col(floor(end/2),:),mat(1)-block(1)+1,mat(2)-block(2)+1,(size(col,2)/(mat(1)-block(1)+1))/(mat(2)-block(2)+1));
    %might need to find another way of finding the representative intensity
else
    error(message('images:col2im:unknownBlockType', deblank(kind)))
end