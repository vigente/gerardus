function [xCoordsNew, yCoordsNew, zCoordsNew] = apply_rotation_matrix_to_vector(SA_Or_LA_World_Coordinates,R)
% apply_rotation_matrix_to_vector Applies a rotation matrix to short axis
% or long axis slice world coordinates.
%
% [XCOORDSNEW, YCOORDSNEW, ZCOORDSNEW] = apply_rotation_matrix_to_vector(SA_OR_LA_WORLD_COORDINATES,R)
% 
%   SA_OR_LA_WORLD_COORDINATES is a structure with x, y, z world
%   components. 
%   
%   R is the rotation matrix
%   
%   XCOORDSNEW is the new calculated  X coordinates
%
%   YCOORDSNEW is the new calculated  Y coordinates
% 
%   ZCOORDSNEW is the new calculated  Z coordinates


% Author: Benjamin Villard <b.016434@gmail.com>
% Copyright © 2014 University of Oxford
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


xCoords = SA_Or_LA_World_Coordinates.coords.x.W(:);
yCoords = SA_Or_LA_World_Coordinates.coords.y.W(:);
zCoords = SA_Or_LA_World_Coordinates.coords.z.W(:); 
coords = [xCoords';yCoords';zCoords'];
testNewCoords = R*coords;
xCoordsNew = reshape(testNewCoords(1,:),size( SA_Or_LA_World_Coordinates.coords.x.W));
yCoordsNew = reshape(testNewCoords(2,:),size( SA_Or_LA_World_Coordinates.coords.y.W));
zCoordsNew = reshape(testNewCoords(3,:),size(SA_Or_LA_World_Coordinates.coords.z.W));

end