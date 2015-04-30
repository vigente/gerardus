function [ MRI_signal, gradient_directions ] = DT_from_smoldyn_pt1( gradient_strength, delta, Delta, positions )
%DT_FROM_SMOLDYN_PT1 calculates the simulated MRI signal from the molecule
%positions.
%
% inputs:
%   gradient_strength: in T/um
%   delta: in us
%   Delta: in us
%   positions: m by 3 by t matrix, where m is number of molecules, t is
%   number of timesteps, 3 gives x,y,z positions.
%
% outputs:
%   MRI_signal: vector of n values, n = no of directions
%   gradient_directions: defined within this function, needed for the
%   second part too
%
% notes:
%   

% Author: Joanne Bates <joanne.bates@eng.ox.ac.uk>
% Copyright (c) 2015 University of Oxford
% Version: 0.1.0
% Date: 28 April 2015
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


% DEFINE SOME PARAMETERS
% input these if haven't saved them
voxel_size = 600;           %um
actual_voxel_size = 200;    %um
% resolution & range of the co-ordinates
xres = 1; 
yres = 1;
zres = 1;
xmax = (voxel_size-1)/2;
ymax = (voxel_size-1)/2;
zmax = (voxel_size-1)/2;
% choose a step time (us)
step_time = 100;
% define the gyromagnetic ratio of water (rad/us/T)
gyro_ratio = 267.5;
% number of timesteps
no_of_timesteps = size(positions, 3);
% echo/measurement time
TE = ((Delta + delta)/step_time) + 1;

% CUT DOWN POSITIONS DATA TO ONLY REQUIRED MOLECULES
% exclude timesteps beyond the echo time
positions = positions(:, :, 1:TE);
% exclude those outside the voxels
% check for molecules finishing outside of voxel
keep1 = positions(:,1,TE)<= (voxel_size-actual_voxel_size);
keep2 = positions(:,1,TE)>= actual_voxel_size;
keep3 = positions(:,2,TE)<= (voxel_size-actual_voxel_size);
keep4 = positions(:,2,TE)>= actual_voxel_size;
keep5 = positions(:,3,TE)<= (voxel_size-actual_voxel_size);
keep6 = positions(:,3,TE)>= actual_voxel_size;
% create a vector of 0s (discard) & 1s (keep)
keep_molecule = keep1.*keep2.*keep3.*keep4.*keep5.*keep6;
% make 0s into nan
keep_molecule(keep_molecule==0) = nan;
positions(:,1,1) = positions(:,1,1).*keep_molecule;
% remove NaNs
positions = positions(isfinite(positions(:,1)),:,:);
% number of molecules
no_of_molecules = size(positions,1);

clear actual_voxel_size voxel_size wall_thickness
clear keep1 keep2 keep3 keep4 keep5 keep6 keep_molecule

% CALCULATION OF DIFFUSION SIGNAL

% as for diffusion signal, but for 6 directions
gradient_directions = [0.7071, 0.7071, 0; 0.7071, -0.7071, 0; 0.7071, 0, 0.7071; 0.7071, 0, -0.7071; 0, 0.7071, 0.7071; 0, 0.7071, -0.7071];
no_of_directions = size(gradient_directions,1);
% parameters independent of gradient direction or voxel
% b_value = gradient_strength^2 * gyro_ratio^2 * delta^2 *(Delta - (delta/3)); (us/um^2)

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
clear sign_x sign_y sign_z h g


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

        molecule_pulse_sum = zeros(no_of_molecules,1);

        % calculate phase shift for each molecule
        for k= 1:no_of_molecules
            molecule_pulse = zeros(no_of_timesteps,1);
            for m = 1:pulse_1_end_step
                molecule_pulse(m) = squeeze(position_grad_dirn(k,:,m))*gradient_strength;
            end
            for m = pulse_2_start_step:pulse_2_end_step
                molecule_pulse(m) = -squeeze(position_grad_dirn(k,:,m))*gradient_strength;
            end
            molecule_pulse_sum(k) = sum(molecule_pulse);
        end

    % calc MRI signal
    phase_change_molecules = gyro_ratio*molecule_pulse_sum*step_time;
    cos_phase_change = cos(phase_change_molecules);
    MRI_signal(g) = mean(cos_phase_change);
end
    

clear ADC_vox DT_vox FA_vox angle_vox eVal eVect step_time 
clear H Y d g v order k m
clear x_pos y_pos z_pos x_1 y_1 z_1 xmax ymax zmax xpos ypos zpos xres yres zres
clear no_of_molecules no_of_timesteps
clear b_value cos_phase_change molecule_pulse molecule_pulse_sum
clear phase_change_molecules position_grad_dirn positions
clear pulse_1_end_step pulse_2_end_step pulse_2_start_step start_corner
clear TE

end

