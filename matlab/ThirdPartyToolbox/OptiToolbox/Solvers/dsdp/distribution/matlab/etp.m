%%***********************************************************
%% etp: Education testing problem.
%%
%% (D) maximize     e'*d
%%     subject to   B - diag(d) >= 0
%%                  d >= 0
%%
%% (P) minimize     Tr B*X
%%     subject to   X >= 0
%%                  diag(X) >= e
%%
%% Ref: M.T. Chu, J.W. Wright, IMA J. of Numerical Anal.,
%%      15 (1995), pp. 141--160.
%%-----------------------------------------------------------
%% [objval,d] = etp(B);
%%
%% B = nxn positive definite.
%%
%% DSDP5.6: version 5.6 
%% Copyright (c) 2005 by
%% S.J. Benson, Y. Ye
%% Last modified: 2 Feb 05
%%***********************************************************

   function [objval,d] = etp(B);

   if (~isreal(B))
      error('only real B allowed');
   elseif (norm(B-B','fro') > 1e-13);
      error(' B must be symmetric'); 
   end;
%%
%% validate B
%%
   n = length(B);  
   d = eig(B); d = real(d);
   if (min(d) < 0); 
      error('B must be positive def'); 
   end;
%%
%%
  
   AC = cell(2,3);
   AC{1,1}='SDP';
   AC{1,2}=n;
   AC{1,3}=sparse(n*(n+1)/2,n+1);
   AC{2,1}='LP';
   AC{2,2}=n;
   AC{2,3}=sparse(n,n+1);

   AC{1,3}(:,n+1)=dvec(B); 
   b = ones(n,1);  

   A = cell(2,n); 
   for k = 1:n 
       AC{1,3}(:,k) = dsparse(k,k,1,n,n); 
       AC{2,3}(:,k) = [zeros(k-1,1); -1; zeros(n-k,1)]; 
   end;  

   y0 = 0.9*min(d)*ones(n,1);
   OPTIONS=doptions;
%   OPTIONS.r0=0;
   OPTIONS.zbar=100*n;
   OPTIONS.print=5;
   [STAT,y,X] = dsdp(b,AC,OPTIONS,y0);
   [pobj,objval,err]=derror(STAT,y,X,b,AC);
   d = y; 
%%===========================================================

