function h = nuimagesc(varargin)
% NUIMAGESC  imagesc() with non-uniform pixel spacing
%
% NUIMAGESC(X, Y, C) 
%
%   Diplays the C array as an image. C intensities will be scaled as in
%   IMAGESC().
%
%   X, Y are vectors with size(C, 2) and size(C, 1), respectively. X, Y
%   provide the axis values of the pixel centres. The vectors need to be
%   monotomic, but not uniform. The spacing between any consecutive pixels
%   needs to be a multiple of the smallest spacing.
%
%   When spacing is not uniform, X and Y will be extended at uniform
%   smallest spacing. Then, C will be linearly interpolated with INTERP2().
% 
% NUIMAGESC(X, Y, C, ...) 
%
%   The same extra arguments for IMAGE() can be passed to NUIMAGESC
%
% See also: interp2, imagesc, image.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2016 University of Oxford
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


% user has provided locations for the pixel centers
if (nargin >= 3 && isvector(varargin{1}) && isvector(varargin{2}))
    
    X = varargin{1}; % pixel centers X-axis
    Y = varargin{2}; % pixel centers Y-axis
    IM = varargin{3}; % image
    
    % user has provided one pixel center per pixel in the image, so
    % potentially the spacing is non-uniform
    if (length(X) > 2)
        
        % minimum spacing in the axis
        dX = min(diff(X));
        
        % recreate the axis with the minimum spacing
        X2 = min(X):dX:max(X);
        
    else
        
        X2 = X;
        
    end
    if (length(Y) > 2)
        
        % minimum spacing in the axis
        dY = min(diff(Y));
        
        % recreate the axis with the minimum spacing
        Y2 = min(Y):dY:max(Y);
        
    else
        
        Y2 = Y;
        
    end
    
    % if any axis has been changed, we need to interpolate the image
    if ((length(X) ~= length(X2)) || (length(Y) ~= length(Y2)))

        % grid from the axis pixel centers
        [Xg, Yg] = meshgrid(X, Y);
        [Xg2, Yg2] = meshgrid(X2, Y2);
        
        % loop each channel of the image
        IM2 = zeros(length(Y2), length(X2), size(IM, 3), 'like', IM);
        for ch = 1:size(IM, 3)
            IM2(:, :, ch) = interp2(Xg, Yg, single(IM(:, :, ch)), Xg2, Yg2);
        end
        
    else
        
        IM2 = IM;
        
    end
    
    % pass the new image to Matlab's imagesc
    if (nargout > 0)
        h = imagesc(X2, Y2, IM2, varargin{4:end});
    else
        imagesc(X2, Y2, IM2, varargin{4:end});
    end
    
else
    
    % pass the image to Matlab's imagesc
    if (nargout > 0)
        h = imagesc(varargin);
    else
        imagesc(varargin);
    end

end

end
