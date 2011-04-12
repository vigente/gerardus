function plot_dmatrix(x, d, LABELS)
% PLOT_DMATRIX  Plot the connections implied by a distance matrix
%
% PLOT_DMATRIX(X)
% PLOT_DMATRIX(X, [])
%
%   This function plots each column of matrix X as a point in 2D or 3D, and
%   lines connecting each point to all others.
%
% PLOT_DMATRIX(X, D)
%
%   If X has N points, D is an (N, N)-matrix with the distances between
%   pairs of points ("Inf" is taken as no connection).
%
% PLOT_DMATRIX(X, D, LABELS)
%
%   LABELS is a flag to indicate whether you want to display numbers close
%   to each point. By default LABELS=false.

% Author: Ramon Casero.
% Copyright © 2011 University of Oxford
% Version: 0.1.1
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
error(nargchk(1, 3, nargin, 'struct'));
error(nargoutchk(0, 0, nargout, 'struct'));

% defaults
if (nargin < 2 || isempty(d))
    d = dmatrix(x);
end
if (nargin < 3 || isempty(LABELS))
    LABELS = false;
end


if (size(x, 1) == 2) % 2D
    
    plot(x(1, :), x(2, :), '.')
    hold on
    plot(x(1, 1), x(2, 1), 'ko')
    for I = 1:size(d, 1)
        for J = 1:size(d, 2)
            if (~isinf(d(I, J)))
                plot([x(1, I) x(1, J)], [x(2, I) x(2, J)])
            end
        end
    end
    % add label to each contour point
    if LABELS
        for J = 1:size(x, 2)
            text(x(1, J), x(2, J), num2str(J))
        end
    end

elseif (size(x, 1) == 3) %3D
    
    plot3(x(1, :), x(2, :), x(3, :), '.')
    hold on
    plot3(x(1, 1), x(2, 1), x(3, 1), 'ko')
    for I = 1:size(d, 1)
        for J = 1:size(d, 2)
            if (~isinf(d(I, J)))
                plot3([x(1, I) x(1, J)], [x(2, I) x(2, J)], [x(3, I) x(3, J)])
            end
        end
    end
    % add label to each contour point
    if LABELS
        for J = 1:size(x, 2)
            text(x(1, J), x(2, J), x(3, J), num2str(J))
        end
    end
    
else
    error('Only 2D or 3D points can be plotted')
end


