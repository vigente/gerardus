function [ positions, voxel_size, step_time, delta, Delta ] = dt_from_smoldyn_pt1( positions, voxel_size, actual_voxel_size, step_time, delta, Delta )
%DT_FROM_SMOLDYN_PT1 Cuts down the positions data to only that needed.
%
% PART 1 OF THE CALCULATION OF THE DIFFUSION TENSOR FROM THE SMOLDYN MODEL DATA
% This part cuts down the positions data to the actual voxel size, and to
% the echo time (TE), i.e. removes data not needed. Only those molecules
% within the voxel at the echo time contribute to the signal measured in
% that molecule. In addition, cutting down the data at this stage allows several
% positions files to be loaded and analysed in one go without reaching the
% memory capacity of Matlab and/or making the process too slow.
%
% Inputs are
%   positions: file of size no of molecules x 3 x no of timesteps generated
%       from the Smoldyn program and converted to matlab format
%   voxel_size: the total (isotropic) size of the volume simulated (um)
%   actual_voxel_size: the (isotropic) size of the voxel to be analysed, at
%       the centre of the total volume (um)
%   step_time: the time of each step as used to generate the Smoldyn diffusion (us)
%   delta: the gradient pulse time (us)
%   Delta: the time between gradient pulses (us)
%   If the data has already been cut down to the required voxel size, then
%   input actual voxel size and voxel size as the same number
%
%  Outputs are
%   positions: file of size no of molecules x 3 x no of timesteps with now
%       only those molecules that finish in the voxel included.
%   voxel_size: as input, included as an output as needed later
%   step_time: as input, included as an output as needed later
%   delta: as input, included as an output as needed later
%   Delta: as input, included as an output as needed later

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

% DEFINE SOME PARAMETERS
% echo/measurement time
TE = ((Delta + delta)/step_time) + 1;

% CUT DOWN POSITIONS DATA TO ONLY REQUIRED MOLECULES
positions = positions(:, :, (end-TE+1):end);

% IF THE VOXEL SIZE HAS ALSO CHANGED
if voxel_size ~= actual_voxel_size
% exclude those outside the voxels
% check for molecules finishing outside of voxel
    keep1 = positions(:,1,end)<= (voxel_size - ((voxel_size-actual_voxel_size)/2));
    keep2 = positions(:,1,end)>= (voxel_size - actual_voxel_size)/2;
    keep3 = positions(:,2,end)<= (voxel_size - ((voxel_size-actual_voxel_size)/2));
    keep4 = positions(:,2,end)>= (voxel_size - actual_voxel_size)/2;
    keep5 = positions(:,3,end)<= (voxel_size - ((voxel_size-actual_voxel_size)/2));
    keep6 = positions(:,3,end)>= (voxel_size - actual_voxel_size)/2;
    % create a vector of 0s (discard) & 1s (keep)
    keep_molecule = keep1.*keep2.*keep3.*keep4.*keep5.*keep6;
    % make 0s into nan
    keep_molecule(keep_molecule==0) = nan;
    positions(:,1,1) = positions(:,1,1).*keep_molecule;
    % remove NaNs
    positions = positions(isfinite(positions(:,1)),:,:);
end

end

