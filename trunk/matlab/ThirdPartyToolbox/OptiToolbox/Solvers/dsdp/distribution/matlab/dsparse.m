%%*******************************************************************
%%  Create a column vector representing a symmetric matrix.
%% 
%%     V = DSPARSE(i,j,s,m,n) uses the same arguments as SPARSE(i,j,s,m,n)
%%
%%     However, this routine assumes the matrix will be symmetric, so
%%     that n must equal m.  The rows of [i,j,s] to generate an
%%     n-by-n sparse matrix with space allocated for length(s) nonzeros.  The
%%     two integer index vectors, i and j, and the real entries
%%     vector, s, all have the same length, which is the number of
%%     nonzeros in the resulting sparse matrix S .  When duplicate values
%%     for a matrix element are provided, that last one will be used.
%%
%%  Output: V = column vector representing half of a symmetric matrix.
%%              A matrix with n rows and columns, implies V is a column 
%%              vector with n*(n+1)/2 rows.  This column vector is in
%%              the format used by DSDP for symmetric data matrices.
%%
%%  Given the nonzero elements (i,j,s) that are on or above the diagonal, 
%%  this routine could be defined as:  V = dvec( sparse(i,j,s,n,n) );
%%   
%%  See also: DVEC
%%
%% DSDP5.6
%% Copyright (c) 2005 by
%% S. Benson and Y. Ye
%% Last modified: December 2005
%%******************************************************************
 

function [V]=dsparse(irow,jcol,nz,n1,n2); 

if (n1~=n2) error('Matrix not square.'); end;
if (size(jcol,1)~=size(irow,1)) error('Number of Row indices do not equal the number of column indices.'); end;
if (size(jcol,2)~=size(irow,2)) error('Number of Row indices do not equal the number of column indices.'); end;
if (size(jcol,1)~=size(nz,1)) error('Number of Row indices do not equal the number of nonzeros.'); end;
if (size(jcol,2)~=size(nz,2)) error('Number of Row indices do not equal the number of nonzeros'); end;

nn=n1*(n1+1)/2;
nnz=length(nz);

V=sparse(nn,1);
idx=zeros(nnz,1);
%% idx=(irow .* irow - irow)/2 + jcol
for k=1:nnz, 
   ii=irow(k);
   jj=jcol(k);
   if (jj>ii) tt=ii;ii=jj;jj=tt; end;
   kk=ii*(ii-1)/2+jj;
   idx(k)=kk;
end;
V(idx)=nz;
