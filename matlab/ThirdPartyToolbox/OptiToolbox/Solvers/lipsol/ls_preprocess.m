function [A,b,c,lb,ub] = ls_preprocess(A,b,c,lb,ub,opts)
% PREPROCESS  - Preprocessing input data.
% Usage:  [A,b,c,lb,ub,probData] = preprocess(A,b,c,lb,ub,probData)

% Yin Zhang, last updated April, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County
% Modified J.Currie AUT May 2013

global probData

verb = opts.verb;

if(verb), fprintf('LIPSOL Preprocessing ...\n'); end

%Ensure columns
if(size(b,2) > 1), b = b'; end
if(size(c,2) > 1), c = c'; end
if(size(lb,2) > 1), lb = lb'; end
if(size(ub,2) > 1), ub = ub'; end

%Check problem sizes
[m,n] = size(A);
if(m ~= size(b,1)), error('A does not have the same number of rows as b'); end
if(n ~= size(c,1)), error('A does not have the same number of colums as c'); end
if(n ~= size(lb,1)), error('lb does not have the same number of elements as c'); end
if(n ~= size(ub,1)), error('ub does not have the same number of elements as c'); end

if any(lb > ub)
   if(verb), fprintf('\nPreprocessor: Lower bound exceeds upper bound\n'); end
   probData.isFeasible = 0; 
   return;
end

%Convert Infinite Bounds to "BIG"
ub(isinf(ub)) = opts.big;
%Check for infinite lower bounds (should be allowed, but no goood...)
if(any(isinf(lb))), throwAsCaller(MException('OPTI:LIPSOLINFLB','LIPSOL only solves problems with Finite Lower Bounds')); end

if(~issparse(A)), A = sparse(A); end;
b = sparse(b); c = sparse(c);
lb = sparse(lb);

%----- delete fixed variables -----
fixed = lb == ub;
probData.Fixed_exist = any(fixed);
if (probData.Fixed_exist)
   ifix = find(fixed); 
   infx = find(1 - fixed);
   xfix = lb(ifix);
   c    = c(infx);
   b    = b - A(:,ifix)*sparse(xfix);
   A    = A(:,infx);
   lb = lb(infx);
   ub = ub(infx);
   if(verb), fprintf(' - %i fixed vars\n', length(ifix)); end
   probData.data_changed = 1;
   probData.ifix = ifix;
   probData.infx = infx;
   probData.xfix = xfix;
end 

%----- delete zero rows -----
rnnzct = sum(spones(A'));
if any(rnnzct == 0)
   izrows = find(rnnzct == 0);
   if any(b(izrows) ~= 0) 
      if(verb), fprintf('\nPreprocessor: problem infeasible\n'); end
      probData.isFeasible = 0; 
      return;
   end;
   inzrows = find(rnnzct > 0);
   A = A(inzrows,:); b = b(inzrows); 
   rnnzct = rnnzct(inzrows);
   if(verb), fprintf(' - %i 0-rows', length(izrows)); end
   probData.data_changed = 1; m = size(A,1);
end

%----- make A structurally "full rank" -----
sprk = sprank(A');
if (sprk < m)
   [dmp, ~] = dmperm(A);
   irow = dmp(1:sprk);
   A = A(irow,:); b = b(irow); 
   rnnzct = rnnzct(irow);
   if(verb), fprintf(' - %i dep-rows\n', m-sprk); end
   probData.data_changed = 1;
end

%----- delete zero columns -----
zrcol = (max(abs(A)) == 0)';
if any(zrcol == 1)
   probData.Zrcols_exist = 1;
   izrcol = find(zrcol);
   if any(c(izrcol) < 0 & ub(izrcol) > opts.big-1)
      if(verb), fprintf('\nPreprocessor: problem unbounded below\n'); end
      probData.isFeasible = 0; 
      return;
   end
   xzrcol = zeros(size(izrcol))...
          + (c(izrcol) < 0).*ub(izrcol)...
          + (c(izrcol) > 0).*lb(izrcol);
   inzcol = find(1 - zrcol);
   A = A(:,inzcol);
   c = c(inzcol);
   lb = lb(inzcol);
   ub = ub(inzcol);
   if(verb), fprintf(' - %i 0-columns\n', nnz(zrcol)); end
   probData.data_changed = 1;
   probData.izrcol = izrcol;
   probData.xzrcol = xzrcol;
   probData.inzcol = inzcol;
end

%----- solve singleton rows -----
singleton = (rnnzct == 1);
nsgrows = nnz(singleton);
if nsgrows >= max(1, .01*size(A,1))
   probData.Sgtons_exist = 1;
   isgrows = find(singleton);
   iothers = find(1 - singleton);
   if(verb), fprintf(' - %i singletons\n',nsgrows); end

   Atmp = A(isgrows,:); Atmp1 = spones(Atmp); btmp = b(isgrows);
   if nsgrows == 1 
      isolved  = find(Atmp1); 
      insolved = find(Atmp1 == 0); 
      xsolved  = b(isgrows)/Atmp(isolved);
   else
      colnnzct = sum(Atmp1);
      isolved  = find(colnnzct);
      insolved = find(colnnzct == 0);
      [ii, jj] = find(Atmp); 
      Atmp = Atmp(ii,jj); btmp = btmp(ii);
      xsolved  = btmp./diag(Atmp);
      if any(colnnzct >  1)
         repeat = diff([0; jj]) == 0;
         for i = 1:length(xsolved) - 1
             if repeat(i+1) && xsolved(i+1) ~= xsolved(i)
                if(verb), fprintf('\nPreprocessor: problem infeasible\n'); end
                probData.isFeasible = 0; 
                return;
             end;
         end;
         ii = find(~repeat); jj = ii;
         Atmp = Atmp(ii,jj); btmp = btmp(ii);
         xsolved  = btmp./diag(Atmp);
      end;
   end;

   if (any(xsolved < lb(isolved)) || any(xsolved > ub(isolved)))
      if(verb), fprintf('\nPreprocessor: problem infeasible\n'); end
      probData.isFeasible = 0; 
      return;
   end;

   b = b(iothers) - A(iothers,isolved)*xsolved;
   A = A(iothers, insolved);
   c = c(insolved);
   lb = lb(insolved);
   ub = ub(insolved);
   probData.data_changed = 1;
   probData.isolved = isolved;
   probData.insolved = insolved;
   probData.xsolved = xsolved;
end

%----- shift nonzero lower bounds -----
probData.Lbounds_non0 = any(lb ~= 0);
if (probData.Lbounds_non0)
   b = b - A*lb;
   probData.data_changed = 1;
end

%----- find upper bounds -----
iubounds = ub < opts.big - 1;
probData.Ubounds_exist = full(any(iubounds));
if (probData.Ubounds_exist)
   ub = sparse(iubounds.*(ub-lb)); 
   probData.nub = nnz(ub);
end

[m, n] = size(A); probData.NNZA = nnz(A);
if(verb), fprintf('[m n] = [%d %d]\n',m,n); end
