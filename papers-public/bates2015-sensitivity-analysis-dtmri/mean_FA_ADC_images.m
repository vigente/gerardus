function [ mean_FA, mean_ADC, mean_bvalue ] = mean_FA_ADC_images(S, S0, dirdata, bdata, mask )
%MEAN_FA_ADC_IMAGES calculation of FA and ADC from diffusion MRI dataset
% 
% Calculate DT, then use DT to get anisotropy measures
% Basic equation to be solved is Y = Hd (solving for d)
% 
% inputs:
%   S: n 3D images, where n is number of different acquisitions
%   S0: 3D image where bvalue = 0
%   dirdata: gradient directions (n by 3)
%   bdata: bvalues (vector of n)
%   mask: 3D mask of 0s and 1s, where 1 = within the left ventricle
%   
% outputs:
%   mean_ADC: mean apparant diffusion coefficient
%   mean_FA: mean fractional anisotropy
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

mean_bvalue = mean(bdata);

no_images = size(bdata,1);

% Create H matrix
H = zeros(no_images,6);
for i = 1:no_images
    h = [dirdata(i,1)^2, dirdata(i,2)^2, dirdata(i,3)^2, 2*dirdata(i,1)*dirdata(i,2), 2*dirdata(i,1)*dirdata(i,3), 2*dirdata(i,2)*dirdata(i,3)];
    H(i,:) = h;
end

% Initialise matrices
FA = zeros(size(S0,1), size(S0,2), size(S0,3));
ADC = zeros(size(S0,1), size(S0,2), size(S0,3));

% For each voxel
for i = 1:size(S0,1)
    for j = 1:size(S0,2)
        for k = 1:size(S0,3)
            % only for voxels within the mask
            if mask(i,j,k) == 1;
            Y = zeros(no_images,1);
            if S0(i,j,k) == 0
               S0(i,j,k) = 10; % to stop getting inf/NaN in later calculations
            end
            for m = 1:no_images
                if S(i,j,k,m) == 0
                   S(i,j,k,m) = 11; % to stop getting inf/NaN in later calculations
                end
                % calculate the Y matrix
                Y(m) = log(S0(i,j,k)./S(i,j,k,m))/bdata(m);
            end
            d = zeros(6,1);
            dv = zeros(3,3);
            % calculate the diffusion tensor
            d = H\Y;                                                
            % transform to a matrix
            dv = [d(1),d(4),d(5);d(4),d(2),d(6);d(5),d(6),d(3)];     
            % calculate the eigenvalues & vectors
            [eVect, eVal] = eig (dv);
            % eig returns a matrix with the eigenvalues on the diagonal, so just take these
            eVal = diag(eVal);
            % Sort eigenvalues from largest to smallest
            [eVal, order] = sort(eVal, 'descend');
            % Remove negative eigenvalues, make them equal to +eps
            % (smallest value in matlab)
            if eVal(3) <= 0
                eVal(3) = +eps;
            end
            if eVal(2) <= 0
                eVal(2) = +eps;
            end
            if eVal(1) <= 0
                eVal(1) = +eps;
            end
            % Order the eigenvectors by the eigenvalues too
            eVect = eVect(:,order);
            % Calculate the ADC
            ADCv = mean(eVal);
            % Calculate the FA
            FAv = sqrt(1.5)*sqrt(((eVal(1)-ADCv)^2+(eVal(2)-ADCv)^2+(eVal(3)-ADCv)^2)/(eVal(1)^2+eVal(2)^2+eVal(3)^2));
            % Store FA & ADC in full matrix
            FA(i,j,k) = FAv;
            ADC(i,j,k) = ADCv;
            else
                FA(i,j,k) = nan;
                ADC(i,j,k) = nan;
            end
        end
    end
end

% calculate means within the mask only    
mean_FA = nanmean(nanmean(nanmean(FA)));
mean_ADC = nanmean(nanmean(nanmean(ADC)));
end

