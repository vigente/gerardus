%%*******************************************************
%% maxcut: MAXCUT problem.
%%
%% (P)    min  Tr C*X
%%        s.t.  diag(X) = e; 
%%                              
%% C = -(diag(B*e)-B)/4. 
%%
%% (dual problem)   max  e'*y
%%                  s.t. diag(y) + S = C. 
%%-------------------------------------------------------
%%
%% [X,y,objval] = maxcut(B);
%%
%% B: weighted adjacency matrix of a graph.
%%
%% DSDP5.0
%% Copyright (c) 2004 by
%% S. Benson, Y. Ye
%% Last modified: Jan 2004
%%*******************************************************

 function [X,y,objval] = maxcut(B);

   if ~isreal(B); error('only real B allowed'); end; 
 
   n = length(B); e = ones(n,1); 
   CC = -(spdiags(B*e,0,n,n)-B)/4.0; 
   b = e;
  
   AC = cell(1,3);
   AC{1,1}='SDP';
   AC{1,2}=n;
   nn=n*(n+1)/2;

   AAC=sparse(nn,n);
   for k = 1:n; AAC(:,k) = dsparse(k,k,1,n,n); end; 
   AAC=[AAC sparse(dvec(CC))];
   AC{1,3}=AAC;
   y0 = -1.1*abs(CC)*e;
   OPTIONS=doptions;
   OPTIONS.gaptol=0.0001;
   OPTIONS.r0=0;
   OPTIONS.print=1;
   OPTIONS.rho=5;
   OPTIONS.zbar=norm(B,1);
   [STAT,y,Xv] = dsdp(b,AC,OPTIONS,y0);

   [objval,dobj,err] = derror(STAT,y,Xv,b,AC);
   X=dmat(Xv{1});
return;
%%=======================================================
