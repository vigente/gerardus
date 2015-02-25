function [DT, FA, ADC, EigVects, EigVals] = fit_DT_YHd_method(s, s0, dirdata, bdata )
%FIT_DT_YHD_METHOD Calculates the diffusion tensor to a set of diffusion
%MRI images
%
% inputs:
%   s = image data (diffusion value) 3D or 4D matrix, with final dimension
%   is the number of images/scans (N)
%   s0 = image data for b-value of 0 (or near 0), i.e. the reference image,
%   2D or 3D image
%   dirdata = unit vectors of the gradient direction [N 3]
%   bdata = b-value for each image [N 1]
%
% outputs:
%   DT = diffusion tensor components (dxx, dxy, dxz, dyy, dyz, dzz)
%   FA = fractional anisotropy
%   ADC = apparant diffusion coefficient
%   EigVects = eigenvectors of the DT
%   EigVals = eigenvalues of the DT
%
% Basic equation to be solved is Y = Hd + noise (solving for d), from
% Kingsley, 'Introduction to DTI mathematics: part III'. Concepts in
% Magnetic Resonance Park A, 2006
% Y is the matrix of the image magnitude, and b-value
% H is the matrix of the diffusion directions
% d is the vector of the components of the diffusion tensor
%
% for the alternative method of fitting DT see fit_DT.m

% Author: Joanne Bates <joanne.bates@eng.ox.ac.uk>
% Copyright � 2015 University of Oxford
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% calculate the number of images
no_images = size(bdata,1);

% Create H matrix from the direction data
H = zeros(no_images,6);
for i = 1:no_images
    h = [dirdata(i,1)^2, dirdata(i,2)^2, dirdata(i,3)^2, 2*dirdata(i,1)*dirdata(i,2), 2*dirdata(i,1)*dirdata(i,3), 2*dirdata(i,2)*dirdata(i,3)];
    H(i,:) = h;
end

% Initialise output matrices
DT = zeros(size(s0,1), size(s0,2), size(s0,3), 6);
FA = zeros(size(s0));
ADC = zeros(size(s0));
EigVects = zeros(size(s0,1), size(s0,2), size(s0,3), 3, 3);
EigVals = zeros(size(s0,1), size(s0,2), size(s0,3), 3);

% For each voxel
% (probably not the tidiest implementation!)
for i = 1:size(s0,1)
    for j = 1:size(s0,2)
        for k = 1:size(s0,3)
            Y = zeros(no_images,1);
            if s0(i,j,k) == 0
               s0(i,j,k) = 10; % to stop getting inf/NaN in later calculations
            end
            for m = 1:no_images
                if s(i,j,k,m) == 0
                   s(i,j,k,m) = 11; % to stop getting inf/NaN in later calculations
                end
                % calculate the Y matrix (log(So/S)/b)
                Y(m) = log(s0(i,j,k)./s(i,j,k,m))/bdata(m);
            end
            d = zeros(6,1);
            dv = zeros(3,3);
            % calculate the diffusion tensor
            d = H\Y;                                                
            % transform to a matrix
            dv = [d(1),d(4),d(5);d(4),d(2),d(6);d(5),d(6),d(3)];     
            % calculate the eigenvalues & vectors
            [eVect, eVal] = eig (dv);
            % eig returns a matrix with the eigenvalues on the diagonal, so
            % just take this diagonal
            eVal = diag(eVal);
            % Sort eigenvalues from largest to smallest
            [eVal, order] = sort(eVal, 'descend');
            % Remove negative eigenvalues, make them equal to +eps
            % they are only there because of noise.
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
            DT(i,j,k,:) = d;
            FA(i,j,k) = FAv;
            ADC(i,j,k) = ADCv;
            EigVects(i,j,k,:,:) = eVect;
            EigVals(i,j,k,:) = eVal;
        end
    end
end
    
end

