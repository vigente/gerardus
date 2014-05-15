%%******************************************************
%% control: a SDP problem from control theory. 
%%
%%  (D)  max   t
%%       s.t   Bk'*P + P*Bk <= 0, k = 1:L
%%             P <= I,   
%%             P >= t*I,  P = P'.
%%
%%------------------------------------------------------
%%
%%  [objval,P] = control(B), 
%%
%%  where B{k} = Bk,  k = 1:L. 
%%  
%%  For example,  B1 = [-0.2    2.2    0.5;
%%                      -2.3   -1.0    2.5;
%%                      -2.1    1.2   -3.0];
%%
%%                B2 = [-1.5    0.6    2.2; 
%%                      -1.2   -2.1    3.3;
%%                       1.9   -3.2    2.9];
%%
%% DSDP: version 5.6 
%% Copyright (c) 2005 by
%% S.J. Benson, Y. Ye
%% Last modified: 2 Feb 05
%%******************************************************

  function [objval,P] = control(B); 

  if ~iscell(B); error('B must be a cell array such that B{k} = Bk'); end; 

  L = length(B); 
  N = length(B{1});

  nn=N*(N+1)/2;

  b = [zeros(nn,1); 1]; 

  AC = cell(L+2,3);
  for l = 1:L+2
     AC{l,1} = 'SDP'; AC{l,2} = N; 
     AC{l,3}=sparse(nn,nn+2);
  end;

  AC{L+1,3}(:,nn+2)=dvec(speye(N));
  AC{L+2,3}(:,nn+1)=dvec(speye(N));
          
  for l = 1:L 
      count = 1;
      for k = 1:N  
          for j = k:N
              G = B{l}; 
              ak = [G(k,:)'];
              aj = [G(j,:)'];
              col1 = [j*ones(N,1)];
              col2 = [k*ones(N,1)];
              row = [1:N]';
              if j ~= k;
                tmp = sparse([row; row],[col1; col2],[ak; aj],N,N);
              elseif (j == k);
                tmp = sparse(row,col1,ak,N,N);
              end;
              AC{l,3}(:,count)=dvec(tmp + tmp');
              count = count+1; 
          end;    
      end;
  end;

  count = 1; 
  for k = 1:N  
      for j = k:N
          ak = [zeros(k-1,1); 0.5 ;zeros(N-k,1)];
          aj = [zeros(j-1,1); 0.5 ;zeros(N-j,1)];

          col1 = [j*ones(N,1)];
          col2 = [k*ones(N,1)];
          row  = [1:N]';

          if j ~= k; 
             tmp1 = sparse([row; row],[col1; col2],[ak; aj],N,N);
             tmp2 = sparse([row; row],[col1; col2],-[ak; aj],N,N); 
          elseif (j == k);
             tmp1 = sparse(row,col1,ak,N,N); 
             tmp2 = sparse(row,col1,-ak,N,N); 
          end; 
          AC{L+1,3}(:,count) = dvec(tmp1 + tmp1');
          AC{L+2,3}(:,count) = dvec(tmp2 + tmp2');

          count = count+1; 
      end;
  end;

  [STAT,y,X]=dsdp(b,AC);
   P=dmat(y(1:nn));
   objval=dot(b,y);

%%======================================================
