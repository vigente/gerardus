function im2 = blockproc3(im, blksz, fun, border, useparallel)
% BLOCKPROC3  Block processing for a 3D image
%
% Sometimes images are too large to be processed in one block. Matlab
% provides function blockproc(), that works on 2D images only. This
% blockproc3() function, on the other hand, works also with 3D images.
%
% IM2 = blockproc3(IM, BLKSZ, FUN)
%
%   IM is a 2D or 3D-array with the input grayscale image.
%
%   IM2 is the ouput image, processed by blocks, with the same size as IM.
%
%   BLKSZ is a 2- or 3-vector with the size of the blocks the image will be
%   split into. The last block in every row, column and slice will be
%   smaller, if the image cannot accommodate it.
%
%   FUN is a function handle. This is the processing applied to every
%   block.
%
% IM2 = blockproc3(IM, BLKSZ, FUN, BORDER)
%
%  BORDER is a 2- or 3-vector with the size of the border around each
%  block. Basically, each block size will be grown by this amount, the
%  block processed, and then the border cut off from the solution. This is
%  useful to avoid "seams" in IM2.
%
% Example:
%
%       fun = @(x) deconvblind(x, ones(20, 20, 10));
%       im2 = blockproc3d(im, [256 256 128], fun, [30 30 10]);
%
% IM2 = blockproc3(..., PARALLEL)
%
%   PARALLEL is a boolean flag to activate parallel processing if the
%   Parallel Computing Toolbox is available. When PARALLEL=true, all cores
%   available locally will start processing the blocks in the image (note
%   that the Parallel Computing Toolbox can only use up to 4 cores at the
%   same time, so if your machine has e.g. 6 cores, 4 cores will be used
%   and 2 cores will remain inactive). Default: PARALLEL=false.
%
%   Note that in parallel model, FUN gets wrapped in parallel jobs. Thus,
%   execution can be interrupted by pressing CTRL+C even when FUN is a C++
%   MEX file.
%
%   Parallel processing works by first splitting the image into blocks with
%   borders, assigning processing of each block to a separate job, then
%   processing each block, and finally trimming the borders and
%   reassembling the results into the output image.
%
%   This will increase the amount of necessary memory, but this function is
%   relatively efficient in the sense that only blocks being processed take
%   up extra memory (blocks waiting in the job queue don't). 
%
%   The speed improvement may not be as good as expected (e.g. 4 processor
%   will make the processing only 2x faster), as Matlab function can be
%   already more or less optimised to make use of several processors.
%
%   To use this option, it is important to NOT have created a pool of
%   workers in Matlab with "matlabpool open", as this would block all local
%   resources, and not leave any cores free for processing the image
%   blocks. A simple example:
%
%      % median filtering using a 19x19x19 voxel neighbourhood, processing
%      % the image by blocks, using parallel computations
%      matlabpool open
%      sz = [19 19 19];
%      fun = @(x) medfilt3(x, sz);
%      im2 = blockproc3(im, [128 128 64], fun, (sz+1)/2, true);
%      matlabpool close
%
%
% See also: blockproc, scimat_blockproc3.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.3.0
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
narginchk(3, 5);
nargoutchk(0, 1);

% defaults
if isempty(blksz)
    blksz = size(im);
end
if (nargin < 4) || isempty(border)
    border = [0 0 0];
end
if (nargin < 5) || isempty(useparallel)
    useparallel = false;
end

% for convenience, we need the size vector to have 3 components, even for a
% 2D image
if (length(blksz) < 3)
    blksz(3) = 1;
end

% image size
imsz = size(im);

% block limits without the extra borders...

% ... starting points
r0 = 1:blksz(1):imsz(1);
c0 = 1:blksz(2):imsz(2);
s0 = 1:blksz(3):imsz(3);

% ... end points
rx = r0 + blksz(1) - 1;
cx = c0 + blksz(2) - 1;
sx = s0 + blksz(3) - 1;
rx = min(rx, imsz(1));
cx = min(cx, imsz(2));
sx = min(sx, imsz(3));

% block limits with the extra borders...

% ... starting points
br0 = max(r0 - border(1), 1);
bc0 = max(c0 - border(2), 1);
bs0 = max(s0 - border(3), 1);

% ... end points
brx = min(rx + border(1), imsz(1));
bcx = min(cx + border(2), imsz(2));
bsx = min(sx + border(3), imsz(3));

% number of blocks in each dimension
NR = length(r0);
NC = length(c0);
NS = length(s0);

% init output
im2 = im;


if (useparallel) % parallel processing
    
    % build a cluster from the local profile
    c = parcluster;
    
    % number of blocks
    numblocks = NR * NC * NS;
    
    % process all input blocks (loops are in inverted order, so that
    % linear indices follow 1, 2, 3, 4...)
    job = cell(1, numblocks);
    for B = 1:numblocks
        
        % get r, c, s indices for current block
        [I, J, K] = ind2sub([NR, NC, NS], B);
        
        % create a new job for this block
        job{B} = createJob(c);
        
        % add processing of current block as a task
        createTask(job{B}, fun, 1, ...
            {im(br0(I):brx(I), ...
            bc0(J):bcx(J), ...
            bs0(K):bsx(K) ...
            )});
        
        % run job
        submit(job{B});
            
    end
    
    % wait on all the jobs to finish, and extract the outputs
    for B = 1:numblocks
        
        % wait until job finishes
        wait(job{B});
        
        % extract the output of this job
        aux = fetchOutputs(job{B});
        
        % assign result to output removing the borders
        im2(...
            r0(I):rx(I), ...
            c0(J):cx(J), ...
            s0(K):sx(K) ...
            ) ...
            = aux{1}(...
            r0(I)-br0(I)+1:rx(I)-br0(I)+1, ...
            c0(J)-bc0(J)+1:cx(J)-bc0(J)+1, ...
            s0(K)-bs0(K)+1:sx(K)-bs0(K)+1 ...
            );

    end
    
else % single processor (we save memory by not creating a cell vector with all the blocks)

    % iterate all image blocks
    for I = 1:NR
        for J = 1:NC
            for K = 1:NS

                % process current image block
                aux = fun(im(...
                    br0(I):brx(I), ...
                    bc0(J):bcx(J), ...
                    bs0(K):bsx(K) ...
                    ));
                
                
                % assign result to output removing the borders
                im2(...
                    r0(I):rx(I), ...
                    c0(J):cx(J), ...
                    s0(K):sx(K) ...
                    ) ...
                    = aux(...
                    (r0(I)-br0(I)+1):(end-(brx(I)-rx(I))), ...
                    (c0(J)-bc0(J)+1):(end-(bcx(J)-cx(J))), ...
                    (s0(K)-bs0(K)+1):(end-(bsx(K)-sx(K))) ...
                    );
                
            end
        end
    end

end

