function im2 = blockproc3(im, blksz, fun, border)
% BLOCKPROC3  Block processing for a 3D image
%
% Sometimes images are too large to be processed in one block. Matlab
% provides function blockproc(), that works on 2D images only. This
% blockproc3() function, on the other hand, works also with 3D images.
%
% The original Matlab function blockproc() can run processing in parallel
% (if the Parallel Computing Toolbox is available). This function doesn't
% support that functionality yet.
%
% IM2 = blockproc3(IM, BLKSZ, FUN)
%
%   IM is a 2D or 3D-array with a grayscale image.
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
%       fun = @(im) deconvblind(im, ones(20, 20, 10));
%       nrrd.data = blockproc3d(nrrd.data, [256 256 128], fun, [30 30 10]);
%
% See also: blockproc

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
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
error(nargchk(3, 4, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% defaults
if isempty(blksz)
    if isstruct(im)
        blksz = size(im.data);
    else
        blksz = size(im);
    end
end
if (nargin < 4) || isempty(border)
    border = [0 0 0];
end

% for convenience, we need the size vector to have 3 components, even for a
% 2D image
if (length(blksz) < 3)
    blksz(3) = 1;
end

% image size
if isstruct(im)
    imsz = size(im.data);
else
    imsz = size(im);
end

% starting points for blocks
r0 = 1:blksz(1):imsz(1);
c0 = 1:blksz(2):imsz(2);
s0 = 1:blksz(3):imsz(3);

% ending points for blocks
rx = r0 + blksz(1) - 1;
cx = c0 + blksz(2) - 1;
sx = s0 + blksz(3) - 1;
rx = min(rx, imsz(1));
cx = min(cx, imsz(2));
sx = min(sx, imsz(3));

% compute block limits taking into account the borders
br0 = max(r0 - border(1), 1);
bc0 = max(c0 - border(2), 1);
bs0 = max(s0 - border(3), 1);

brx = min(rx + border(1), imsz(1));
bcx = min(cx + border(2), imsz(2));
bsx = min(sx + border(3), imsz(3));

% init output
im2 = im;

% iterate all image blocks
for I = 1:length(r0)
    for J = 1:length(c0)
        for K = 1:length(s0)
            
            % process current image block
            if isstruct(im)
                aux = fun(im.data(...
                    br0(I):brx(I), ...
                    bc0(J):bcx(J), ...
                    bs0(K):bsx(K) ...
                    ));
            else
                aux = fun(im(...
                    br0(I):brx(I), ...
                    bc0(J):bcx(J), ...
                    bs0(K):bsx(K) ...
                    ));
            end
            
            
            % assign result to output removing the borders
            if isstruct(im)
                im2.data(...
                    r0(I):rx(I), ...
                    c0(J):cx(J), ...
                    s0(K):sx(K) ...
                    ) ...
                    = aux(...
                    (r0(I)-br0(I)+1):(end-(brx(I)-rx(I))), ...
                    (c0(J)-bc0(J)+1):(end-(bcx(J)-cx(J))), ...
                    (s0(K)-bs0(K)+1):(end-(bsx(K)-sx(K))) ...
                    );
            else
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
