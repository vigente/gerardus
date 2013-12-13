function [xsol,objp,info,msg,times,x,y,z,s,w] = lipsol(A,b,c,lb,ub,opts)
% LIPSOL        - Main program for LIPSOL

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County
% Modified J.Currie AUT May 2013

global probData

if(nargin < 6)
    opts = lipsolset; 
    opts.maxiter = 100;
    opts.maxtime = 1000;
    opts.tol = 1e-5;
end
if(nargin < 5), ub = []; end
if(nargin < 4), lb = []; end
if(nargin < 3), error('You must supply at least 3 arguments to lipsol(A,b,c)'); end

times = zeros(5,1);
t0 = tic;

%Problem Data Structure to avoid so many Global Variables
probData = struct('NNZA',0,'NNZL',0,'XLNZ',[],'PERM',[],'INVP',[],'LNZ',[],'LOOP_LEVEL',[],...
                  'XSUPER',[],'XLINDX',[],'LINDX',[],'SNODE',[],'SPLIT',[],'TMPSIZ',[],...
                  'LNZ0',[],'b_orig',b,'c_orig',c,'data_changed',0,'isFeasible',1,...
                  'Lbounds_non0',[],'Ubounds_exist',0,'nub',0,'Fixed_exist',0,...
                  'ifix',[],'infx',[],'xfix',[],'Zrcols_exist',0,'izrcol',[],...
                  'inzcol',[],'xzrcol',[],'Sgtons_exist',0,'isolved',[],...
                  'insolved',[],'xsolved',[],'nt',0,'message',[],...
                  'col_scaled',0,'colscl',[],'Dense_cols_exist',0,'Mdense',[],...
                  'idense',[],'ispars',[],'Hist',[],'Sherman_OK',[],'backed',[]);
%Initial Time              
opts.t0 = tic;

%Save loop_level into problem structure
probData.LOOP_LEVEL = opts.loop_level;
          
%Pretty Printing
if(opts.verb), fprintf('----------------------------------------------\n'); end
          
%Preprocess the Input Data
[A,b,c,lb,ub] = ls_preprocess(A,b,c,lb,ub,opts); 
times(1) = toc(t0);

%Check if preproccesor indicates infeasible
if(~probData.isFeasible)
    xsol = []; objp = []; info = [-1 0 0]; msg = 'Presolve Indicated Infeasible'; times = [];
    x = []; y = []; z = []; s = []; w = [];
    return;
end

%Scale the Input Data if Required
[A,b,c,ub] = scaling(A,b,c,ub,opts);
%Determine density of the problem
checkDense(A,opts);
times(2) = toc(t0) - times(1);

%Solve the Problem
[x,y,z,s,w,info] = ls_miip(A,b,c,ub,opts);  
times(3) = toc(t0) - times(2);

%Display Solve Status
if(opts.verb && ~isempty(probData.message)), fprintf('\n  Status: %s\n',probData.message); end

%Post Process the Solution
[xsol, objp, times ] = ls_postprocess(x,y,lb,opts.verb,times,t0);
msg = probData.message;

%Pretty Printing
if(opts.verb), fprintf('----------------------------------------------\n'); end


function [A,b,c,ub] = scaling(A,b,c,ub,opts)
% SCALING     - Scale matrix A. Usage: [A,b,c,ubounds] = scaling(A,b,c,ubounds)
% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global probData

[m, n] = size(A);
badscl = 1.e-4;

absnzs = abs(nonzeros(A));
thescl = min(absnzs)/max(absnzs);
clear absnzs

if (thescl < badscl)
    if(opts.verb), fprintf('Scaling ...\n'); end

    % ----- scaling vectors ------
    AA = abs(A);
    colscl = full(sqrt(max(AA)'));
    rowscl = full(sqrt(max(AA,[],2)));
    clear AA;

    % ----- column scaling -----
    if (probData.Ubounds_exist), ub = ub.*colscl; end
    colscl = reciprocal(colscl);
    A = A*sparse(1:n,1:n,colscl);
    c = c.*colscl;
    probData.col_scaled = 1;

    % ----- row scaling -----
    rowscl = reciprocal(rowscl);
    A = sparse(1:m,1:m,rowscl)*A;
    b = b.*rowscl;
    bnrm = norm(b);
    if (bnrm > 0) 
      q = median([1 norm(c)/bnrm 1.e+8]);
      if (q > 10), A = q*A; b = q*b; end
    end
    probData.data_changed = 1;
    probData.colscl = colscl;
end


function checkDense(A,opts) %#ok<INUSD>
% CHECKDENSE - Determine and locate dense columns of matrix A.
% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

%PROBLEM IN PCG - SKIPPING DENSE COLUMN MANAGEMENT [J.C. 19/5/13]

% global probData
% 
% [m, n] = size(A);
% probData.ispars = 1:n; 
% nzratio = 1;
% if (m >  500), nzratio = 0.20; end;
% if (m > 1000), nzratio = 0.10; end;
% if (m > 5000), nzratio = 0.05; end;
% 
% if (nzratio < 1)
%    checking = sum(spones(A))/m <= nzratio;
%    if any(checking == 0)
%       probData.Dense_cols_exist = 1;
%       probData.idense = find(sparse(1-checking));   % Dense  column indices
%       probData.ispars = find(checking);             % Sparse column indices
%       if(opts.verb), fprintf('Dense columns (nnz/m > %g): %i\n',nzratio,length(probData.idense)); end
%    else
%       if(opts.verb), fprintf('No Significant Dense Columns\n'); end
%    end
%    clear checking
% else
%     if(opts.verb), fprintf('Problem is Small Enough to Ignore Dense Columns\n'); end
% end


function Y = reciprocal(X)
% RECIPROCAL  - Invert the nonzero entries of a matrix elementwise.
%             Y = RECIPROCAL(X) has the same sparsity pattern as X
%	      (except possibly for underflow).

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

if issparse(X)
   [m, n]  = size(X);
   [i,j,Y] = find(X);
   Y = sparse(i,j,1./Y,m,n);
else
   Y = 1./X;
end
ceiling = 1.e+16; Y = min(ceiling,Y);