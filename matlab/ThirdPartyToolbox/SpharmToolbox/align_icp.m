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
% align_icp.m 
%
% Goal: Align P to X using ICP
%
% Li Shen 
% 07/29/2003 - customized using Hany Farid's icp code.

function [P,M] = align_icp(P,X)

M = eye(4);

%%% 
%%% ICP - ALIGN P TO X
%%%
while(1)
   for k = 1 : length(P)
      dist = (X(:,1)-P(k,1)).^2 + (X(:,2)-P(k,2)).^2 + (X(:,3)-P(k,3)).^2;
      [val,ind] = min(dist);
      closest(k) = ind;
   end
   Xc = X(closest,:);
   
   muP   = mean(P);
   muXc  = mean(Xc);
   N     = length(P);
   sig   = zeros(3);
   for k = 1 : N
      sig = sig + (P(k,:)' * Xc(k,:));
   end
   sig   = sig/N;
   sig   = sig - (muP' * muXc); % cross-covariance matrix
   A     = sig - sig';
   delta = [A(2,3) ; A(3,1) ; A(1,2)];
   Q     = [trace(sig)  delta' ; delta  sig+sig'-trace(sig)*eye(3)];
   [V,D] = eig(Q);
   [val,ind] = max(diag(D));
   qR    = V(:,ind); % optimal rotation
   q0    = qR(1);
   q1    = qR(2);
   q2    = qR(3);
   q3    = qR(4);
   R     = [q0^2+q1^2-q2^2-q3^2  2*(q1*q2 - q0*q3)  2*(q1*q3 + q0*q2) ; ...
	    2*(q1*q2 + q0*q3)  q0^2+q2^2-q1^2-q3^2  2*(q2*q3 - q0*q1) ; ...
	    2*(q1*q3 - q0*q2) 2*(q2*q3 + q0*q1)  q0^2+q3^2-q1^2-q2^2];
   qT    = muXc' - R*muP'; % optimal translation
   
   Pnew  = (R*P')';
   diff  = mean( sum(((Pnew - P).^2)') );

   if( diff < 0.1 )
      fprintf( '%f\n', diff );
      break;
   end
   P = Pnew;
   
   P(:,1) = P(:,1) + qT(1);
   P(:,2) = P(:,2) + qT(2);
   P(:,3) = P(:,3) + qT(3);
   
   Mr = eye(4); Mr(1:3,1:3) = R;
   Mt = eye(4); Mt(1:3,4) = qT(1:3);
   M = Mt*Mr*M;
end

