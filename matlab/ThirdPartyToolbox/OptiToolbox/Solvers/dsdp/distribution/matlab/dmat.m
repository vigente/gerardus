%%*******************************************************************
%%  Convert a vector matrix used by DSDP to a square symmetric matrix.
%%
%%  [A] = dmat(V)
%%
%%  Input: V = column vector representing upper triangle of a symmetric matrix.
%%
%%  Output:  A = square symmetric matrix
%%
%%  See also: DVEC (its inverse operator)
%%
%% DSDP5.0 
%% Copyright (c) 2003 by
%% S. Benson and Y. Ye
%% Last modified: December 2003
%%******************************************************************

function [A]=dmat(V); 

[m,m1]=size(V);
if (m1~=1) error('Not a column vector.'); end;

n=floor(sqrt(2*m));
if (n*(n+1)/2 ~= m) 
     error('Impossible Dimension. Cannot convert to square matrix.'); 
end;

nzV=nnz(V);
if (issparse(V) | nzV < n*(n+1)/2) A=sparse(n,n); else A=zeros(n,n);  end;

for i=1:n, 
    k1=i*(i-1)/2+1;  k2=i*(i+1)/2;
    A(1:i,i)= V(k1:k2); 
end;
A=A+triu(A,1)';
