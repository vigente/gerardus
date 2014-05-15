function ls_symbfct(P,opts)
% SYMBFCT     - Symbolic factorization for sparse Cholesky factorization.
%            symbfct(P) does a minimum degree ordering and a symbolic
%            factorization for a symmetric square sparse matrix P, and 
%            passes the results as global variables.
%            It calls two Fortran MEX files.

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County
% Modified J.Currie AUT May 2013

% PATCH by Jos F. Sturm, May 1998.
% 1. Use global variable "OUTFID" for sending messages. (NOT USED)
% 2. Don't call symfct.mex* if P is diagonal, since that would
%    cause a segmentation violation.

global probData

verb = opts.verb;

% ----- check input matrix P -----
m = size(P,1);
if(size(P,2) ~= m || ~issparse(P) || isempty(P))
   error('Matrix must be square, sparse and nonempty.');
end;
Pdiag = diag(diag(P)); P = P - Pdiag;

% ----- min-degree ordering -----
t0 = tic;
if(verb), fprintf('  Minimum-Degree Ordering ... '); end
[probData.PERM, probData.INVP] = ls_ordmmd(P);
if(verb), fprintf('Done. [%1.3fs]\n',toc(t0)); end

% ----- symbalic factorization -----
if(nnz(P) > 0)
    t0 = tic;
    if(verb), fprintf('  Symbolic Factorization ... '); end
    [probData.XLNZ,probData.NNZL,probData.XSUPER,probData.XLINDX,probData.LINDX,...
     probData.SNODE,probData.SPLIT,probData.TMPSIZ] = ls_symfct(P,probData.PERM,probData.INVP,opts.cache_size);
    if(verb), fprintf('Done. [%1.3fs]\n',toc(t0)); end
end
