% test_scinrrd_surface_interpolation.m
%
% Script to test scinrrd_surface_interpolation(), that first creates an
% interpolated surface with surface_interpolation()

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.1.1
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

% create a NRRD struct skeleton to generate the segmentation mask
res = [.01 .01 .01];
offset = [0 0 0];
nrrd = scinrrd_im2nrrd(zeros(101, 101, 101, 'uint8'), res, offset);

% points that define the surface
x = [...
    .25 .25 .1; ...
    .75 .25 .5; ...
    .25 .75 .8; ...
    .75 .75 .2 ...
    ];

% interpolate surface
param.type = 'xy';
interp.type = 'mbae';
nrrd = scinrrd_surface_interpolation(nrrd, x, param, interp);

% plot isosurface of interpolated surface
[gx, gy, gz] = scinrrd_ndgrid(nrrd);
close
p = patch(isosurface(gx, gy, gz, nrrd.data));
set(gca, 'FontSize', 18)
set(p,'FaceColor','red','EdgeColor','none');
daspect([1 1 1])
view(3); axis tight
view(-84, 26)
camlight('headlight')
lighting gouraud
xlabel('x')
ylabel('y')
zlabel('z')
