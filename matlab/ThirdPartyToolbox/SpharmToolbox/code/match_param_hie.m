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

%
% match fvec to atlas by rotating the parameterization of fvec
% - object does not transform in the object space
% - alpha, beta, gamma: samples in the rotation space
% - gran = # of alphas to process together
% - res: base res R + step of hierarchy Hs + depth of hierarchy Hd + top N 
%

function [best_fvec, best_angl] = match_param_hie(fvec,atlas,A,B,G,gran,res,max_d)

% base resolution R, step of hierarchy Hs, depth of hierarchy Hd, top N 
R = res(1); Hs = res(2); Hd = res(3); N = res(4);
disp(sprintf('----- Base Res %d, Step %d, Depth %d -----',R,Hs,Hd));
alpha = A{1}; beta = B{1}; gamma = G{1};
[fvec, Aprm] = match_param(fvec,atlas,alpha,beta,gamma,gran,N,max_d); 
fvec = fvec(:,:,1); Aprm = Aprm(1,:); % remove later on
best_angl = Aprm;
for i=1:Hd
    disp(sprintf('----- Res %d -----',R+Hs*i));
    alpha = A{i+1}; beta = B{i+1}; gamma = G{i+1};
    [fvec, Aprm] = match_param(fvec,atlas,alpha,beta,gamma,gran,1,max_d);
    best_angl(i+1,:) = Aprm;
end
best_fvec = fvec;
    
return;
