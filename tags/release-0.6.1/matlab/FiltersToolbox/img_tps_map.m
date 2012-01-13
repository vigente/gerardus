function im = img_tps_map(s, t, im)
% IMG_TPS_MAP Warp an image using a thin-plate spline transformation
%
% IM2 = IMG_TPS_MAP(S, T, IM)
%
%   IM is an image.
%
%   S is a (P,Ds,N)-volume where each (:,:,i)-matrix has the coordinates
%   of the source points that define the warp.
%
%   T is a (P,Dt,N)-volume where each (:,:,i)-matrix has the coordinates
%   of the target points that define the warp.
%
%     (P is the number of points, Ds and Dt are the dimension and N is the
%     number of configurations).
%
%   Note that the thin-plate spline (TPS) warp is defined from the target
%   to the source image. So the warp applied to the source image is the
%   inverse of the TPS, that has no known explicit formulation (in
%   particular, the inverse of a TPS is not a TPS).

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
error(nargchk(3, 3, nargin));
error(nargoutchk(0, 1, nargout));

% get image size
[nr, nc] = size(im);

% create uniform grid for the target image
[gy, gx] = ndgrid(1:nr, 1:nc);

% warp grid points. The way to do this is warp the target uniform grid onto
% the source image. Ideally, we would like to warp with the inverse of the
% TPS, but there is no known expression for it. So the only thing we can do
% is to warp the grid with the TPS from target to source. Thus, the image
% will be warped from the source to the target using the implicit inverse
% of the TPS
gxyw = pts_tps_map(t, s, [gx(:) gy(:)]);

% resize grid points
gxw = reshape(gxyw(:, 1), [nr, nc]);
gyw = reshape(gxyw(:, 2), [nr, nc]);
clear gxyw;

% interpolate intensity values in the target image
im = interpn(gy, gx, im, gyw, gxw, 'linear');
