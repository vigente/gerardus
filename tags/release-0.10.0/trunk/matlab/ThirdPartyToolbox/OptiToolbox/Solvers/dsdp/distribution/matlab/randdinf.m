%%*******************************************************************
%% randdinfsdp.m : creates random dual infeasible SDP problems with 
%%                 various block diagonal structures. 
%%
%% [X,y] = randdinfsdp(blocksize,m);
%%
%% E.g.
%%      randdinfsdp([[10 5],10);
%%
%%  blocksize : for generating dense blocks, where
%%              blocksize(i) is the dimension of the ith block
%%
%% DSDP Version 5.6 
%% Copyright (c) 2005 by
%% S.J. Benson and Y. Ye
%% Last modified: 4 Feb 05
%%******************************************************************

   function  [X,y] = randinf(blocksize,m);

   if nargin < 2; error(' insufficient number of inputs '); end;   

   P=length(blocksize);
   AC=cell(P,3);

   for k = 1:P
       n=blocksize(k);
       AC{k,1}='SDP';
       AC{k,2}=n;
       AC{k,3}=sparse(n*(n+1)/2,m+1);
   end;


%% dual infeasible
%%
%% generate infeasibility certificate Xi, positive definite X0, and C
%%
   Xi = [];
   X0 = [];
   trXX=0;
   for k = 1:P;
              n = blocksize(k);
              tmp = randn(n); 
              tmp = tmp*tmp'; 
              tmp = 0.5*(tmp + tmp');
              %mineig = min(real(eig(tmp)));
              %if (mineig < 0); tmp = tmp - 1.1*mineig*eye(n); end;
              Xi{k} = tmp;
              trXX=trXX+sum(sum(tmp.*tmp));
              tmp = randn(n); 
              tmp = tmp*tmp'; 
              tmp = 0.5*(tmp + tmp');
              %mineig = min(real(eig(tmp)));
              %if (mineig < 0); tmp = tmp - 1.1*mineig*eye(n); end;
              X0{k} = tmp;
              tmp = randn(n); tmp = 0.5*(tmp + tmp');
              AC{k,3}(:,m+1)=dvec(tmp);
   end;

%%
%% set up the matrices A and b
%%
   A  = cell(1,m);
   b = ones(m,1);
   for k = 1:m;
       trAx=0;bb=0;
       for j = 1:P;
              n = blocksize(j);
              Aj  = randn(n);
              Ak = (Aj + Aj')/2;
              trAx = trAx + sum(sum(Ak .* Xi{j})); 
              Ak=Ak - trAx/trXX*Xi{j};
              Ak=(Ak+Ak')/2;
              AC{j,3}(:,k)=dvec(Ak);
              bb=bb+ sum(sum(Ak .* X0{j})); 
       end;
       b(k) = bb;
   end;

%   if (trCX>0) for j=1,P, A{j,3}(:,m+1)=-A{j,3}(:,m+1);  end;
%%
%% Xi is a dual infeasibility certificate
%%

%%
  OPTIONS=doptions;
  OPTIONS.print=1;
  [STAT,y,Xv] = dsdp(b,AC,OPTIONS);
  [pobj,dobj,err] = derror(STAT,y,Xv,b,AC);

%%=================================================
