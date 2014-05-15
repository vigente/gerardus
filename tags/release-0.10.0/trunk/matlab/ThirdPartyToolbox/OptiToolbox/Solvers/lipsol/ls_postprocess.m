function [x, objp, times] = ls_postprocess(x,y,lb,verb,times,t0)
% POSTPROCESS - Recover the original variables for the solution.
%	Usage: [x, objp] = postprocess(x,lbounds)

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County
% Modified J.Currie AUT May 2013

global NNZL probData

if(verb), fprintf('\nLIPSOL Postprocessing ...\n'); end

if (probData.col_scaled),   x = x.*probData.colscl;   end;
if (probData.Lbounds_non0), x = x + lb; end;

if (probData.Sgtons_exist)
   tmp(probData.insolved) = x;
   tmp(probData.isolved) = probData.xsolved;
   x = sparse(tmp);
end;

if (probData.Zrcols_exist)
   tmp(probData.inzcol) = x;
   tmp(probData.izrcol) = probData.xzrcol;
   x = sparse(tmp);
end;

if (probData.Fixed_exist)
   tmp(probData.infx) = x; 
   tmp(probData.ifix) = probData.xfix;
   x = sparse(tmp);
end;

if (size(x,1) < size(x,2))
   x = x';
end;

objp = full(probData.c_orig'*x);

%Finalize Times
times(4) = toc(t0) - times(3);

%Publish Problem Data
if(verb)    
    fprintf('[m n] = [%i %i], nnz(A) = %i, nnz(L) = %i\n', length(y),length(x),probData.NNZA,NNZL);

    if(probData.Lbounds_non0),      fprintf('Nonzero lower bounds exist.\n'); end
    if(probData.Ubounds_exist),     fprintf('Number of Upper bounds: %g\n',probData.nub); end
    if(probData.Fixed_exist),       fprintf('Number of fixed variables: %g\n',length(probData.ifix)); end
    if(probData.Zrcols_exist),      fprintf('Number of zero columns in A: %g\n',length(probData.izrcol)); end
    if(probData.Dense_cols_exist),  fprintf('Number of dense columns in A: %g\n',length(probData.idense)); end

    fprintf('CPU seconds: %6.2f ... preprocessing\n', times(1)+times(2));
    fprintf('             %6.2f ... solving\n', times(3));
    fprintf('             %6.2f ... postprocessing\n', times(4));
    fprintf('             %6.2f ... total\n', sum(times));
end