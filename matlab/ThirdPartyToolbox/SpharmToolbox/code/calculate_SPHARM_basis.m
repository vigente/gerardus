%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherical Harmonic Modeling and Analysis Toolkit (SPHARM-MAT) is a 3D 
% shape modeling and analysis toolkit. 
% It is a software package developed at Shenlab in Center for Neuroimaging, 
% Indiana University (SpharmMat@gmail.com, http://www.iupui.edu/~shenlab/)
% It is available to the scientific community as copyright freeware 
% under the terms of the GNU General Public Licence.
% 
% Copyright 2009, 2010, ShenLab, Center for Neuroimaging, Indiana University
% 
% This file is part of SPHARM-MAT.
% 
% SPHARM-MAT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SPHARM-MAT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with SPHARM-MAT. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ============================================
%
% Goal: Create canonical spherical harmonic bases
%
% Li Shen 
% 04/11/2002 - create
% 10/15/2002 - rename and modify
% 11/03/2008 - renamed by Sungeun Kim.

function Z = calculate_SPHARM_basis(vs, degree)

[PHI,THETA] = cart2sph(vs(:,1),vs(:,2),vs(:,3));
ind = find(PHI<0);
PHI(ind) = PHI(ind)+2*pi;
THETA = pi/2-THETA;
vertnum = size(THETA,1);

Z = spharm_basis(degree,THETA,PHI); 

return;


function Z = spharm_basis(max_degree,theta,phi)

Z = []; vnum = size(theta,1);

% save calculations for efficiency
for k = 0:(2*max_degree)
    fact(k+1) = factorial(k);
end
for m = 0:max_degree
    exp_i_m_phi(:,m+1) = exp(i*m*phi);
    sign_m(m+1) = (-1)^(m);
end

for n = 0:max_degree

	% P = legendre(n,X) computes the associated Legendre functions of degree n and 
	% order m = 0,1,...,n, evaluated at X. Argument n must be a scalar integer 
	% less than 256, and X must contain real values in the domain -1<=x<=1.
	% The returned array P has one more dimension than X, and each element
	% P(m+1,d1,d2...) contains the associated Legendre function of degree n and order
	% m evaluated at X(d1,d2...).

    Pn = legendre(n,cos(theta'))';
    
    posi_Y = [];
    nega_Y = [];
    
    m= 0:n;
    v = sqrt(((2*n+1)/(4*pi))*(fact(n-m+1)./fact(n+m+1)));
    v = v(ones(1,vnum),:).*Pn(:,m+1).*exp_i_m_phi(:,m+1);
    posi_Y(:,m+1) = v; % positive order;
    nega_Y(:,n-m+1) = sign_m(ones(1,vnum),m+1).*conj(v); % negative order
    
    Z(:,end+1:end+n) = nega_Y(:,1:n);
    Z(:,end+1:end+n+1) = posi_Y(:,1:(n+1));
end

return;

