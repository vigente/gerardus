function [reshaped3D2D] = im2col3(varargin)
%IM2COL3 rearranges 3D image blocks into columns. 

% 16 September 2015
%
% im2col3 does the same as im2col, only for 3D arrays.
% The inputs are a 3D dataset and a vector for the blocksize, [x y z], 
% and the kind of block ('sliding' or 'distinct')
% x,y&z have to be divisors of the respective dimension
% The output is a 2D array of dimensions xyz by the number of blocks 
% Each column of xyz elements represents a 3D element
%
% Author: Nicolas Basty <nicolas.basty@eng.ox.ac.uk>
% Copyright 2015 University of Oxford
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

[image3D, block, kind] = parse_inputs(varargin{:});

blocksize_x = block(1);
blocksize_y = block(2);
blocksize_z = block(3);

blocknumx = length(1:blocksize_x:size(image3D,1)); %to find the number of blocks in each dimension
blocknumy = length(1:blocksize_y:size(image3D,2));
blocknumz = length(1:blocksize_z:size(image3D,3));
blocknum = blocknumx*blocknumy*blocknumz; %total number of cubes in the 3D array


if  numel(size(image3D)) ~= 3 %checks if input is 3D array
    fprintf('Input array is expected to be 3D, not %dD! \n',numel(size(image3D)))
    return
end

if  rem(blocksize_x,1) ~= 0 || rem(blocksize_y,1) ~= 0 || rem(blocksize_z,1) ~= 0 %checks if integer
    fprintf('Blocksize is expected to be an integer \n')
    return
end

if strcmp(kind,'distinct')
    
    reshaped3D2D = zeros(blocksize_x*blocksize_y*blocksize_z,blocknum); %initialising final array to store all the blocks in column form
    
    tic
        
    bl=1; %this is the n-th cube number, used to index the final array
    
    for slicenum=1:blocksize_z:(size(image3D,3))
        
        slab = image3D(:,:,slicenum:slicenum+blocksize_z-1);
        % cuts out a slab of dimensions of a 2D slice and height of blocksize
        
        for edgx=1:blocksize_x:(size(image3D,1)) %"edg" for edge in x, y & z directions
            for edgy=1:blocksize_y:(size(image3D,2))
                for edgz=1:blocksize_z
                    
                    reshaped3D2D(blocksize_x*blocksize_y*(edgz-1)+1:blocksize_x*blocksize_y*edgz,bl)=...
                        reshape(slab(edgx:(edgx+blocksize_x-1),edgy:(edgy+blocksize_y-1),edgz),[blocksize_x*blocksize_y 1]);
                    %reshape each slice of the block and store within the final array
                end
                bl=bl+1; %proceed to next block after the Z dimension loop is done
            end
        end
    end
        
end

if strcmp(kind,'sliding')   
    
    bl=1; %this is the n-th cube number, used to index the final array

    for slicenum=1:((size(image3D,3))-block(3)+1)
        %loop ends when the last cube reaches the end
        
        slab = image3D(:,:,slicenum:slicenum+blocksize_z-1);
        % cuts out a slab of dimensions of a 2D slice and height of blocksize        
        
        for edgx=1:1:(size(image3D,1)- blocksize_x +1) 
            %"edg" for edge of the cube in x, y & z directions
            %loop ends when the last cube reaches the end

            for edgy=1:1:(size(image3D,2)- blocksize_y +1)
                for edgz=1:blocksize_z                    
        
                    reshaped3D2D(blocksize_x*blocksize_y*(edgz-1)+1:blocksize_x*blocksize_y*edgz,bl)=...
                        reshape(slab(edgx:(edgx+blocksize_x-1),edgy:(edgy+blocksize_y-1),edgz),[blocksize_x*blocksize_y 1]);
                    %reshape each slice of the block and store within the final array
                
                end
                
                bl=bl+1; %proceed to next block after the Z dimension loop is done
            
            end
        end
    end
end

%%% Function parse_inputs (taken from im2col)
function [a, block, kind, padval] = parse_inputs(varargin)

narginchk(2,4);

switch nargin
    case 2
        if (strcmp(varargin{2},'indexed'))
            error(message('images:im2col:tooFewInputs'))
        else
            % IM2COL(A, [M N])
            a = varargin{1};
            block = varargin{2};
            kind = 'sliding';
            padval = 0;
        end
        
    case 3
        if (strcmp(varargin{2},'indexed'))
            % IM2COL(A, 'indexed', [M N])
            a = varargin{1};
            block = varargin{3};
            kind = 'sliding';
            padval = 1;
        else
            % IM2COL(A, [M N], 'kind')
            a = varargin{1};
            block = varargin{2};
            kind = validatestring(varargin{3},{'sliding','distinct'},mfilename,'kind',3);
            padval = 0;
        end
        
    case 4
        % IM2COL(A, 'indexed', [M N], 'kind')
        a = varargin{1};
        block = varargin{3};
        kind = validatestring(varargin{4},{'sliding','distinct'},mfilename,'kind',4);
        padval = 1;
        
end

if (isa(a,'uint8') || isa(a, 'uint16'))
    padval = 0;
end