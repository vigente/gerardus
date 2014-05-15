%%*******************************************************************
%%  Convert a square symmetric matrix to a vector 
%%  representation used by DSDP
%%
%%  [V] = dvec(A)
%%
%%  Input:  A = square symmetric matrix.  Lower triangular entries 
%%              are ignored.
%%
%%  Output: V = column vector representing half of a symmetric matrix.
%%              If A has n rows and columns, V is a column vector with
%%              n*(n+1)/2 rows.
%%
%%  See also: DSPARSE, DMAT (its inverse operator)
%%
%% DSDP5.0 
%% Copyright (c) 2003 by
%% S. Benson and Y. Ye
%% Last modified: December 2003
%%******************************************************************

function [V]=dvec(A); 

[n,n1]=size(A);

if (n~=n1) error('Matrix not square.'); end;

nn=n*(n+1)/2;

V=zeros(nn,1);

for i=1:n, 
    k1=i*(i-1)/2+1;  k2=i*(i+1)/2;
    V(k1:k2,1) = A(1:i,i); 
end;

V=sparse(V);
