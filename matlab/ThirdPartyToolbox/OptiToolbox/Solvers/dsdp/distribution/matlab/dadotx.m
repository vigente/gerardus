%%*******************************************************************
%%  Apply the operators A to the matrix X.  Each column 
%%  of A and X is the vector representation of a square symmetric 
%%  matrix.
%%
%%  [v] = dadotx(A,X)
%%
%%  Input:  A = square symmetric matrices in DSDP dvec form
%%  Input:  X = square symmetric matrix in DSDP dvec form
%%
%%  Output: v = column vector representing the inner products
%%
%%  See also: DVEC, DMAT
%%
%% DSDP5.0 
%% Copyright (c) 2003 by
%% S. Benson and Y. Ye
%% Last modified: December 2003
%%******************************************************************

function [b]=dadotx(A,x); 

[m,m1]=size(x);
[p,p1]=size(A);
if (m1~=1) error('X Not a column vector.'); end;
if (m~=p) 
     error('Number of rows in first and second argument not equal.'); 
end;

n=floor(sqrt(2*m));
if (n*(n+1)/2 ~= m) 
     error('Impossible Dimension for square matrix.'); 
end;


xx=x';
for i=1:n, 
  k=i*(i+1)/2;
  xx(k) = xx(k)/2; 
end;

b=[xx]*A;
b=full(2*b');
