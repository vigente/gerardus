function [ MRI_signal_noisy ] = add_Rician_noise_to_simulated_dMRI_data( MRI_signal, SNR )
%ADD_RICIAN_NOISE_TO_SIMULATED_DMRI_DATA adds Rician noise to simulated
%MRI_signals.
%
% inputs: noisefree MRI_signal and SNR
% outputs: noisy MRI_signal
%
% noisy signal, S = |S* + N1 + iN2|
% where S* = noise-free signal
% N1 and N2 are random numbers, drawn from normal distribution, with mean=0, and std = 1/SNR.

% Author: Joanne Bates <joanne.bates@eng.ox.ac.uk>
% Copyright (c) 2015 University of Oxford
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

no_of_directions =  size(MRI_signal); 
N1 = (1/SNR).*randn(no_of_directions);
N2 = (1/SNR).*randn(no_of_directions);

MRI_signal_noisy = abs(MRI_signal + N1 + 1i.*N2);
end

