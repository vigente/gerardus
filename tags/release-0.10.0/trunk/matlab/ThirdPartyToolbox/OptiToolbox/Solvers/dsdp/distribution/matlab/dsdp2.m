%*****************************************************************************
% DSDP5:  Dual Scaling Algorithm for Positive Semidefinite Programming
% Copyright (c) 2002 by
% S. J. Benson, Y. Ye
% Last modified: 20 January 2004
%*****************************************************************************
%
% Converts data from DSDP2 format to DSDP5 format and calls solver.
%
%
% > DSDP(C,A,b) attempts to solve the positive semidefinite program
%      MINIMIZE trace(C*X) SUCH THAT A(:,i)'*X*A(:,i) = b_i, i=1,...,m and X >= 0
%      using a dual scaling algorithm.  The first argument is an n x n symmetric
%      matrix.  The second argument is an n x m matrix, and the third 
%      argument b is a dense column vector of length m.  
%
% > DSDP(C,A,b,y0) specifies an initial dual vector y0.
%
% > DSDP(C,A,b,CC,AA) adds additional linear constraints to the dual problem.
%
% > DSDP(C,A,b,CC,AA,y0) specifies an initial dual vector y0.
%
% > [Y] = DSDP() returns the dual solution vector.
%
% > [Y,X] = DSDP() returns the primal and dual solution
%
% > [Y,X,V] = DSDP() returns the rank reduction vector.
%
%
%*****************************************************************************

function [Y,DY,K,V] = dsdp2(C,A,b,CC,AA,y0);

  m=length(b);

  AC=cell(2,3);
  n=size(C,1);
  nn=n*(n+1)/2;
  m=length(b);
  AAC=sparse(nn,m+1);
  for i=1:m, ai=A(:,i); AAC=[AAC dvec(ai*ai')]; end;
  AAC=[AAC dvec(C)];
  AC{1,1}='SDP';
  AC{1,2}=size(C,1);
  AC{2,1}='LP';
  AC{2,2}=0;

  if nrhs>4
     [n1,n2]=size(C{j});
     if (n1==1 | n2==1)
       AAC=sparse([]);
       for i=1:m, AAC=[AAC sparse(AA(:,i)]; end;
       AAC=[AAC sparse(CC)'];
       AC{2,1}='LP';
       AC{2,2}=length(CC);
       AC{2,3}=AAC;
       end;
  end;
       

  [STAT,y,X]=dsdp(b,AC,OPTIONS,y0);

  XX=cell(p,1);
  for j=1:p,
     [n1,n2]=size(C{j})
     if (n1==1 | n2==1)
       XX{j}=X{j};
     else
       XX{j}=dmat(X{j});
     end;
  end;

return;
