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
% rotate_param_m03.m
%
% Goal: rotating parameter space
% Use : m03 < m02 < m01 < s03 < s01 < s00
% m03 : batch of s03
% m02 : modified m03
% m01 : batch of s01 
% s03 : use symmetry
% s01 : set up matrix D
% s00 : original version 
%
% SPHARM rotation using only coefficients
% -- factorial(x): fact(x+1)
% -- x0(m,n,l) = sqrt((l+n)!(l-n)!(l+m)!(l-m)!): 
%                <l(0:max_d), m(0:l), n(-m:m)> => x0{l+1}(m+1, n+l+1)
% -- d(l,m,n,beta) = d(l,-n,-m,beta): yes (but different from standard axes)
% -- d(l,m,n,beta) = (-1)^(m+n)d(l,n,m,beta): yes
% alpha, beta, gamma: vectors (t loop)
%
% Li Shen 
% 01/04/2007 - create

function [rvec] = rotate_param_m03(alpha, beta, gamma, fvec, max_degree)
global fact;
global sgm;

beta  = -beta; % caused by nTHETA = pi/2-THETA (original nTHETA = pi/2+THETA)?

fnum = (max_degree+1)^2;
for k=1:length(beta)
    D{k} = sparse(fnum,fnum);
end
for l = 0:max_degree
    M = [];
    for m = -l:l     % actual order (-l => l)
        ma = m+l; % l^2+ma+1
        for n = -m:l % actual order (-l => l) 
            d_beta = zeros(size(beta));
            na = n+l; % l^2+na+1
			for t = max(0,n-m):min(l+n,l-m)
              d_beta_temp = (-1)^(t)*(sqrt(fact(l+n+1)*fact(l-n+1)*fact(l+m+1)*fact(l-m+1))) ...
                            /(fact(l+n-t+1)*fact(l-m-t+1)*fact(t+m-n+1)*fact(t+1));
              d_beta_temp = d_beta_temp * (cos(beta/2)).^(2*l+n-m-2*t) .* (sin(beta/2)).^(2*t+m-n);
              d_beta = d_beta + d_beta_temp;
			end
            % use symmetry: d(l,m,n,beta) = d(l,-n,-m,beta); d(l,m,n,beta) = (-1)^(m+n)d(l,n,m,beta)
            M(ma+1,na+1,:) = exp(-i*m*alpha).*d_beta.*exp(-i*n*gamma); 
        end
    end 
    idx = (l^2+1):(l+1)^2; 
    for k=1:length(beta)
        if l>0
            Mx = M(:,end:-1:1,k).*(1-eye(size(M(:,:,k)))); Mx = Mx(end:-1:1,:);
            M(:,:,k) = M(:,:,k) - real(Mx).*sgm{l} + i*imag(Mx).*sgm{l};
        end       
        D{k}(idx,idx) = M(:,:,k); 
    end
end
for k=1:length(beta)
    rvec(:,:,k) = D{k}*fvec;
end

return;
