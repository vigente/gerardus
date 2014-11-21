function [ DT, FA, ADC, VectorF ] = fit_DT( im, b, thresh_val )
% FIT_DT    Fit the diffusion tensor model voxelwise to an image
%
%   IM is the input image, of any dimensionality, with the diffusion 
%   scans in the last dimension
%     
%   B is the b-matrix, of size [3 3 N]. If the measured b matrix isn't 
%   available, it can be approximated by:
%   b(:,:,n)=b_values(n)*unit_vectors(n,:)'*unit_vectors(n,:);
% 
%   THRESH_VAL is a threshold for skipping the FA, ADC and vector field 
%   computation in voxels x where im(x) < thresh_val.
%
%   DT is the diffusion tensor coefficients
% 
%   FA is the fractional anisotropy
% 
%   ADC is the apparent diffusion coefficient
% 
%   VECTORF is the vector field
%
% See also DT_to_image

    
% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.1.2
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

% check arguments
narginchk(2, 3);
nargoutchk(0, 4);

if nargin < 3
    thresh_val = -inf;
end
    
sz = size(im);

% take the log of the image to linearise the equation
imlog = log(im);

% Sort all b matrices in to a vector Bv=[Bxx,2*Bxy,2*Bxz,Byy,2*Byz,Bzz];
Bv=squeeze([b(1,1,:),2*b(1,2,:),2*b(1,3,:),b(2,2,:),2*b(2,3,:),b(3,3,:)])';

% If you have unit vectors and b values instead of b matrices, do this:
% b=zeros([3 3 size(unit_vectors,1)]);
% for i=1:size(unit_vectors,1),
%     b(:,:,i)=b_values(i)*unit_vectors(i,:)'*unit_vectors(i,:);
% end


% Add another column to Bv to handle the constant term:
% Slog = Bv * M + log(S0)
% becomes:
% Slog = [Bv, -1] * [M; -log(S0)]
Bv = [Bv, -ones(size(Bv,1), 1)];

% reshape Ilog so we can use the efficient \ operator
imlog = reshape(imlog, [prod(sz(1:end-1)), sz(end)]);

% Fit the tensor model
M = (Bv \ -imlog')';

% In certain cases, you might want to do a nonlinear fit instead:
% x are the coefficients of M
% modelfun = @(x,Bv) exp(-(x(1) * Bv(:,1) + x(2) * Bv(:,2) + x(3) * Bv(:,3) + x(4) * Bv(:,4) + x(5) * Bv(:,5) + x(6) * Bv(:,6) + x(7) * Bv(:,7)));
%x0 = M;
%mdl = nlinfit(Bv, squeeze(im),modelfun,x0);
% You might alternatively want to constrain the tensor to be positive
% (with either option, this will need to go into a for loop, and it will be really slow)
% M = lsqnonneg(Bv, -imlog')'

% return the diffusion tensor
DT = reshape(M, [sz(1:end-1), 7]);

M(isnan(M)) = 0;
M(isinf(M)) = 0;

if nargout > 1 % if you want the FA, ADC, etc.

    % initialise variables
    FA = zeros(size(M,1), 1);
    ADC = zeros(size(M,1),1);
    VectorF = zeros(size(M,1), 3);
   
    for i = 1:size(M,1)
    
        if im(i) < thresh_val
            continue
        end
        
        Mi = M(i,:);

        % The DiffusionTensor (Remember it is a symetric matrix,
        % thus for instance Dxy == Dyx)
        DiffusionTensor=[Mi(1) Mi(2) Mi(3); Mi(2) Mi(4) Mi(5); Mi(3) Mi(5) Mi(6)];

        % Calculate the eigenvalues and vectors, and sort the 
        % eigenvalues from small to large
        [EigenVectors,D]=eig(DiffusionTensor); 
        EigenValues=diag(D);
        [~,index]=sort(EigenValues); 
        EigenValues=EigenValues(index); 
        EigenVectors=EigenVectors(:,index);
        EigenValues_old=EigenValues;

        % Regulating of the eigen values (negative eigenvalues are
        % due to noise and other non-idealities of MRI)
        EigenValues=abs(EigenValues);
        
        % Apparent Diffuse Coefficient
        ADCv=(EigenValues(1)+EigenValues(2)+EigenValues(3))/3;

        % Fractional Anistropy (2 different definitions exist)
        % First FA definition:
        %FAv=(1/sqrt(2))*( sqrt((EigenValues(1)-EigenValues(2)).^2+(EigenValues(2)-EigenValues(3)).^2+(EigenValues(1)-EigenValues(3)).^2)./sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2) );
        % Second FA definition:
        FA(i)=sqrt(1.5)*( sqrt((EigenValues(1)-ADCv).^2+(EigenValues(2)-ADCv).^2+(EigenValues(3)-ADCv).^2)./sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2) );
        ADC(i)=ADCv;
        VectorF(i,:)=EigenVectors(:,end)*EigenValues_old(end);
        
    
    end
    
    % reshape to match input dimensions
    FA = reshape(FA, sz(1:end-1));
    ADC = reshape(ADC, sz(1:end-1));
    VectorF = reshape(VectorF, [sz(1:end-1), 3]);
    
end


