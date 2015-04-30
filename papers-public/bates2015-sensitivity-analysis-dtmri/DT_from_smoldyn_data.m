function [ DT, ADC, FA, angle ] = DT_from_smoldyn_data( gradient_strength, delta, Delta )
%DT_FROM_SMOLDYN_DATA calculates the diffusion tensor from the positions
%data.
%
% positions data needs to be done in 2 parts, as files are too large to
% import all at the same time.
%  
% inputs:
%   gradient_strength: in T/um
%   delta: in us
%   Delta: in us
%   
% outputs:
%   DT: diffusion tensor
%   ADC: apparant diffusion coefficient
%   FA: fractional anisotropy
%   angle: of 1st eigenvector wrt x axis, equivalent to fibre angle
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
     
% part 1
load('positions1.mat');
positions1 = positions;
load('positions2.mat');
positions = [positions1; positions];
clear positions1

[ MRI_signal, gradient_directions ] = DT_from_smoldyn_pt1( gradient_strength, delta, Delta, positions );
MRI_signal1 = MRI_signal;
clear MRI_signal

% part 2
load('positions3.mat');
positions1 = positions;
load('positions4.mat');
positions = [positions1; positions];
clear positions1

[ MRI_signal, gradient_directions ] = DT_from_smoldyn_pt1( gradient_strength, delta, Delta, positions );
MRI_signal2 = MRI_signal;
clear MRI_signal

MRI_signal = (MRI_signal1 + MRI_signal2)./2;
    
clear MRI_signal1 MRI_signal2
    
[ DT, ADC, FA, angle ] = DT_from_smoldyn_pt2( MRI_signal, gradient_directions, gradient_strength, delta, Delta);
        

end

