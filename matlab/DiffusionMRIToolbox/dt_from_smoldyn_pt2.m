function [ MRI_signal, bvalue ] = dt_from_smoldyn_pt2( positions, voxel_size, step_time, delta, Delta, gradient_directions, gradient_strength, bvalue )
%DT_FROM_SMOLDYN_PT2 Calculates the MRI signal
%
% PART 2 OF THE CALCULATION OF THE DIFFUSION TENSOR FROM THE SMOLDYN MODEL DATA
% This part calculates the diffusion MRI signal from the positions data.
% Note that providing each set of positions data contains the same number
% of molecules, the MRI signal from different sets of positions data can be
% combined by just taking the mean.
%
% Inputs are:
%   positions: matrix of size no of molecules x 3 x no of timesteps generated
%       from the Smoldyn program and cut down in part 1
%   voxel_size: scalar. the total (isotropic) size of the volume simulated (um)
%   step_time: scalar. the time of each step as used to generate the Smoldyn diffusion (us)
%   delta: scalar. the gradient pulse time (us)
%   Delta: scalar. the time between gradient pulses (us)
%   gradient_directions: matrix of size no of directions x 3 which gives
%       the x,y,z components of the gradient directions
%   gradient_strength: IF UNKNOWN, SET TO 0 & IT WILL BE CALCULATED FROM
%       BVALUE. scalar. (T/um)
%   bvalue: IF UNKNOWN, SET TO 0 & IT WILL BE CALCULATED FROM
%       GRADIENT_STRENGTH. scalar. (um^2/us)
%
% Outputs are:
%   MRI_signal: relative signal for each of the gradient directions
%   bvalue: as input, included as an output as needed later

% Author: Jo Bates <jobates81@gmail.com>
% Copyright © 2014 University of Oxford
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

no_of_directions = size(gradient_directions,1);
no_of_molecules = size(positions,1);
% number of jumps
no_of_timesteps = size(positions, 3);

% define the gyromagnetic ratio of water (rad/us/T)
gyro_ratio = 267.5;

% calculate the bvalue or gradient_strength if required
if bvalue == 0;
    bvalue = gradient_strength^2 * gyro_ratio^2 * delta^2 *(Delta - (delta/3));
end
if gradient_strength == 0;
    gradient_strength = sqrt(bvalue/(gyro_ratio^2 * delta^2 *(Delta - (delta/3))));
end

% resolution & range of the co-ordinates
xres = 1; 
yres = 1;
zres = 1;
xmax = (voxel_size-1)/2;
ymax = (voxel_size-1)/2;
zmax = (voxel_size-1)/2;

% parameters independent of gradient direction or voxel
pulse_1_end_step = (delta/step_time) + 1;
pulse_2_start_step = (Delta/step_time)+1;
pulse_2_end_step = ((Delta+delta)/step_time) + 1;

% calc parameters which are independent of voxel
start_corner = zeros(no_of_directions, 3);
for g = 1:no_of_directions  
    % work out the starting corner of the image, where the MRI pulse with be 0
    sign_x = sign(gradient_directions(g,1));
    sign_y = sign(gradient_directions(g,2));
    sign_z = sign(gradient_directions(g,3));
    if sign_x < 0
        x_1 = xmax*2+1;
    else
        x_1 = 1;
    end
    if sign_y < 0
        y_1 = ymax*2+1;
    else
        y_1 = 1;
    end
    if sign_z < 0
        z_1 = zmax*2+1;
    else
        z_1 = 1;
    end
    start_corner (g,:) = [x_1, y_1, z_1];
end

MRI_signal = zeros(no_of_directions,1);
% for every gradient direction
for g = 1:no_of_directions
    x_1 = start_corner (g,1);
    y_1 = start_corner (g,2);
    z_1 = start_corner (g,3);
        % for each molecule, for each position, the position in gradient
        % direction is the dot product of the position times the gradient
        % direction. since A.B = a1*b1 + a2*b2 + a3*b3, calculate like this as
        % need to include the potential for different resolutions in each
        % direction and can't seem to do dot product of correct parts for these
        % 2 matrices.

        x_pos = (positions(:,1,:)-x_1).*xres.*gradient_directions(g,1);
        y_pos = (positions(:,2,:)-y_1).*yres.*gradient_directions(g,2);
        z_pos = (positions(:,3,:)-z_1).*zres.*gradient_directions(g,3);
        position_grad_dirn = x_pos + y_pos + z_pos;

        % calculate phase shift for each molecule
        mol_pulse = zeros(size(position_grad_dirn));
        mol_pulse(:,:,1:pulse_1_end_step) = position_grad_dirn(:,:,1:pulse_1_end_step).*gradient_strength;
        mol_pulse(:,:,pulse_2_start_step:pulse_2_end_step) = -position_grad_dirn(:,:,pulse_2_start_step:pulse_2_end_step).*gradient_strength;
        molpulsesum = sum(mol_pulse,3);
       
        % calc MRI signal
        phasech = gyro_ratio*molpulsesum*step_time;
        cosphase = cos(phasech);
        MRI_signal(g) = mean(cosphase);
end

end

