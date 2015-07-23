function [reshaped3D2D,blocknumx,blocknumy,blocknumz] = im2col3(image3D,blocksize)
%IM2COL3 rearranges 3D image blocks into columns

% 21 July 2015
%
% im2col3 does the same as im2col, only for 3D arrays
% It decomposes the 3D array into cubes, then to 2D (layer by layer)
% The inputs are a 3D dataset and a number used for the wished blocksize
% Blocksize has to be a common divisor of the three dimensions
% The output is a 2D array of dimensions blocksize^3 by blocknumber
%
% Author: Nicolas Basty <nicolas.basty@eng.ox.ac.uk>
% Copyright 2015 University of Oxford
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

tic

if  numel(size(image3D)) ~= 3 %checks if input is 3D array
    fprintf('Input array is expected to be 3D, not %dD! \n',numel(image3D))
    return
end

if mod(size(image3D,1),blocksize)+ mod(size(image3D,2),blocksize)+mod(size(image3D,3),blocksize)~= 0
    fprintf('Blocksize is expected to be a common divisor of length of the sides of the array! \n')
    return
    %checks if blocksize is common divisor of the lengths of sides of the array
end

if  rem(blocksize,1) ~= 0 %checks if integer
    fprintf('Blocksize is expected to be an integer \n')
    return
end

blocknumx = length(1:blocksize:size(image3D,1)); %to find the number of blocks in each dimension
blocknumy = length(1:blocksize:size(image3D,2));
blocknumz = length(1:blocksize:size(image3D,3));

blocknum = blocknumx*blocknumy*blocknumz; %total number of cubes in the 3D array

fprintf('%d %dx%dx%d blocks inside the volume \n',blocknum, blocksize,blocksize,blocksize)

bl=1; %this is the n-th cube number, used to index the final array
reshaped3D2D = zeros(blocksize^3,blocknum); %initialising final array to store all the blocks in column form

for slicenum=1:blocksize:(size(image3D,3))
    
    slab = image3D(:,:,slicenum:slicenum+blocksize-1);
    % cuts out a slab of dimensions of a 2D slice and height of blocksize
    
    for vex=1:blocksize:(size(image3D,1)) %"ve" for vertex, x, y & z directions
        for vey=1:blocksize:(size(image3D,2))
            for vez=1:blocksize
                reshaped3D2D(blocksize^2*(vez-1)+1:blocksize^2*vez,bl)=...
                    reshape(slab(vex:(vex+blocksize-1),vey:(vey+blocksize-1),vez),[blocksize^2 1]);
                %reshape each slice of the block and store within the final array
            end
            bl=bl+1; %proceed to next block after the Z dimension loop is done
        end
    end
end

A = size(image3D,1); %for printout
B = size(image3D,2);
C = size(image3D,3);
D = blocksize^3;
E = blocknum;
T = toc;

fprintf('%dx%dx%d 3D-array is now %dx%d \nTime elapsed is %2.2f seconds\n',A,B,C,D,E,T)

end