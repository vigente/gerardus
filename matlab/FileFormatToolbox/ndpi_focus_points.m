function x = ndpi_focus_points(file)
% NDPI_FOCUS_POINTS  Read focus points coordinates from .ndpi image file
% obtained from a Hamamatsu microscope.
%
% X = ndpi_focus_points(FILE)
%
%   FILE is a string with the path and name of an .ndpi file. This is the
%   format used by Hamamatsu microscopes.
%
%   X is a 3-column matrix with the coordinates of the focus points used by
%   the microscope to acquire the image in FILE (in case of multi-image,
%   focus points for the first image, which is the one with the highest
%   resolution, are acquired). Units are meters.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

% check arguments
narginchk(1, 1);
nargoutchk(0, 1);

% read header from Hamamatsu .ndpi file
info = tiff_read_header(file);

% extract field with the coordinates of the focus points as a vector
x = info{1}.private_65439;

% number of focus points
N = length(x) / 3;

% reshape into xyz-coordinates and make them type double
x = reshape(double(x), 3, N)';

% remove unused focus points. Those are the ones with special values
% encoded in the z-coordinate:
%  0 Yellow (Unchecked) 
% -1 Red (Failed) 
% -2 Dark Green (Not used) 
% -3 Dark Yellow (Not checked, focus area not used) 
% -4 Dark Red (Failed, focus area not used) 
idx = ismember(x(:, 3), [-4:1:0]);
x(idx, :) = [];

% Coordinates are in nanometers. Convert the units to
% meters
x = x * 1e-9;

% % DEBUG: plot focus points
% hold off
% plot3(x(:, 1)*1e3, x(:, 2)*1e3, x(:, 3)*1e3, '.')
% xlabel('x (mm)')
% ylabel('y (mm)')
% zlabel('z (mm)')
