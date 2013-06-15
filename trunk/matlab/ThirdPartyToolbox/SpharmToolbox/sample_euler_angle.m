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
% hierarchy of euler angle samples in the rotation space
%

function [A, B, G, samples] = sample_euler_angle(res)

% base resolution R, step of hierarchy Hs, depth of hierarchy Hd 
R = res(1); Hs = res(2); Hd = res(3);
[alpha, beta, gamma] = euler_angle(R,0);
A{1} = alpha; B{1} = beta; G{1} = gamma;
samples = length(alpha)*length(gamma);
for i=1:Hd
    Rp = R;
    R = R + Hs;
    [alpha, beta, gamma] = euler_angle(R,Rp);
    A{i+1} = alpha; B{i+1} = beta; G{i+1} = gamma;
    samples = samples + length(alpha)*length(gamma);
end

return;

%
% euler angle samples in the rotation space
%

function [alpha, beta, gamma] = euler_angle(res,pres)

% % number of vertices at different level of subdivision
% vnum = [42 162 642  2562 10242 40962 163842]
% fnum = [80 320 1280 5120 20480 81920 327680];
% vnum = fnum/2+2

% gamma is set up in the way that high resolution samples contain low
% resolution samples
fnum = 80*4.^(0:10); vnum = fnum/2+2; knum = round(3.5*2.^(0:10)); k = knum(res+1);

gamma = 0:2*pi/k:2*pi; gamma = gamma(1:end-1); % alpha could be any value that requires gamma to rotate it back

[vs, fs] = sample_Param_Space(-res);

if pres==0 % include all thetas
    theta = Inf;
else % keep only the top thetas
	pvs = vs(1:vnum(pres),:);
	[PHI,THETA] = cart2sph(pvs(:,1),pvs(:,2),pvs(:,3));
	THETA = pi/2-THETA;
	theta = min(THETA(2:end));
end

[PHI,THETA] = cart2sph(vs(:,1),vs(:,2),vs(:,3));
ind = find(PHI<0); PHI(ind) = PHI(ind)+2*pi;
THETA = pi/2-THETA;

idx = find(THETA<=theta);
THETA = THETA(idx);
PHI = PHI(idx);

alpha = -PHI;
beta = -THETA;

%disp(sprintf('%d*%d = %d samples in rotation space',length(alpha),k,k*length(alpha)));

return;