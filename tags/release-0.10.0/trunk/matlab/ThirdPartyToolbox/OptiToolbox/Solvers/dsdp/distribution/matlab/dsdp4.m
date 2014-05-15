%*****************************************************************************
% DSDP5:  Dual Scaling Algorithm for Positive Semidefinite Programming
% Copyright (c) 2002 by
% S. J. Benson, Y. Ye
% Last modified: 20 January 2004
%*****************************************************************************
%
% Converts data from DSDP4 format to DSDP5 format, calls solver, and 
% and converts solution to DSDP4 format.
%
%
% > DSDP(A,C,b) attempts to solve the positive semidefinite program
%      MINIMIZE trace(C*X) SUCH THAT trace(A_i*X) = b_i, i=1,...,m and X >= 0
%      using a dual scaling algorithm.  For a problem with p blocks and m
%      constraints, A is a p x m cell array and C is a p x 1 cell array.
%      One block may contain LP variables, and the cells corresponding to
%      this block should be a one dimensional array.  All other cells
%      must contain a square, symmetric, real valued matrix.  The third 
%      argument b is a dense column vector of length m.  
%
% > DSDP(A,C,b,OPTIONS,y0) specifies an initial dual vector y0.
%
% > [STAT,y,X] = DSDP() returns a structure containing relevant statistics, 
%                and approximate dual and primal solutions. 
%
%*****************************************************************************

function [STAT,y,XX] = dsdp4(A,C,b,OPTIONS,y0);

  p=length(C);
  m=length(b);

  AC=cell(p,3);
  for j=1:p,
     [n1,n2]=size(C{j});
     if (n1==1 | n2==1)
       AAC=sparse(n1*n2,m+1);
       for i=1:m, AAC=[AAC sparse(A{j,i})']; end;
       AAC=[AAC sparse(C{j})'];
       AC{j,1}='LP';
       AC{j,2}=length(C{j});
       AC{j,3}=AAC;
     else
       AAC=sparse(n1*(n1+1)/2,m+1);
       for i=1:m, AAC=[AAC dvec(A{j,i})]; end;
       AAC=[AAC dvec(C{j})];
       AC{j,1}='SDP';
       AC{j,2}=size(C{j},1);
       AC{j,3}=AAC;
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
