function [sol] = ls_linsys(P,rhs,save_factor)
% LINSYS     - Solve a positive definite system of linear equations.
%        LINSYS(P,rhs) returns a "solution" to linear system Px = rhs,
%        where the matrix P is sparse symmetric and positive definite.
%        The function SYMBFCT must be called at least once before 
%        LINSYS is called.  LINSYS uses three Fortran MEX programs.
%
%        If P = [], LINSYS will try to use, if any, an existing Cholesky 
%        factor.

% Yin Zhang, July, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County
% Modified J.Currie AUT May 2013

% PATCH by Jos F. Sturm, May 1998.
% 1. If P is diagonal, then we don't have a symbolic factorization.
%    In this case, we solve the trivial system, without using mex-files.

global probData

% ----- check input rhs -----
if issparse(rhs)
   error('RHS vector must be in full matrix format.');
end;
m = length(probData.PERM);
if (max(size(rhs)) ~= m) || (min(size(rhs)) ~= 1)
   error('No symbolic factor or input sizes mismatch.');
end;

if ~isempty(P)
% ----- check input matrix P -----
   if (size(P,1) ~= m || size(P,2) ~= m || ~issparse(P))
      error('Input matrix must be square and sparse.');
   end;
% ----- Remove diagonal from P ------
   Pdiag = full(diag(P));
   P = P - sparse(1:m,1:m,Pdiag);
% ----- Cholesky factorization ------
    if nnz(P) > 0
       LNZ = ls_inpnv(Pdiag,P,probData.INVP,probData.PERM,probData.XLNZ,probData.XSUPER,probData.XLINDX,probData.LINDX,probData.NNZL);
       probData.LNZ = ls_blkfct(probData.XLNZ,probData.XSUPER,probData.SNODE,probData.SPLIT,probData.XLINDX,probData.LINDX,LNZ,probData.TMPSIZ,probData.LOOP_LEVEL);
    else
       probData.LNZ = Pdiag;
    end
end
if(nargin == 3 && save_factor), probData.LNZ0 = probData.LNZ; end

% ------ back solve ------
sol = zeros(m,1);
if length(probData.LNZ) == m        % diagonal solve
  sol(probData.LNZ>0) = rhs(probData.LNZ > 0) ./ probData.LNZ(probData.LNZ > 0);
else
  sol(probData.PERM) = ls_blkslv(probData.XLNZ,probData.XSUPER,probData.XLINDX,probData.LINDX,probData.LNZ,rhs(probData.PERM));
end
