%%*******************************************************
%% gpp: graph partitioning problem.
%%
%% (P) minimize    Tr C*X
%%     subject to  diag(X) = e; 
%%                 Tr (ee'*X) = gamma
%%                              
%% C = -(diag(B*e)-B)/4. 
%%
%% (D) maximize    e'*y + gamma*t
%%     subject to  diag(y) + t*ee' + S = C. 
%%-------------------------------------------------------
%%
%% [X,y,objval] = maxcut(B,gamma);
%%
%% B: weighted adjacency matrix of a graph.
%%
%% DSDP5.0
%% Copyright (c) 2004 by
%% S. Benson, Y. Ye
%% Last modified: Jan 2004
%%*******************************************************

 function [X,y,objval] = gpp(B,gamma);

   if ~isreal(B); error('only real B allowed'); end; 
 
   n = length(B); e = ones(n,1); 
   CC = -(spdiags(B*e,0,n,n)-B)/4.0; 
   b = [e; gamma];
  
   AC = cell(1,3);
   AC{1,1}='SDP';
   AC{1,2}=n;
   nn=n*(n+1)/2;

   AAC=sparse(nn,n);
   for k = 1:n; AAC(:,k) = dsparse(k,k,1,n,n); end; 
   AAC=[AAC dvec(ones(n)) dvec(CC)];
   AC{1,3}=AAC;
   y0 = -1.1*abs(CC)*e;
   y0=[y0; 0];
   OPTIONS=doptions;
   OPTIONS.gaptol=0.0001;
   OPTIONS.r0=0;
   OPTIONS.print=5;
   OPTIONS.rho=5;

   [STAT,y,Xv] = dsdp(b,AC,OPTIONS,y0);

   [objval,dobj,err] = derror(STAT,y,Xv,b,AC);
   X=dmat(Xv{1});
return;
%%=======================================================
