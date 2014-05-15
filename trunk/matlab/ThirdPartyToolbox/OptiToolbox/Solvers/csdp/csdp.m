% CSDP  Solve a SDP using CSDP [Standard Primal Form]
%
% SDP Form:
%  (P)  MINIMIZE   dot(f,x) SUCH THAT  SUM A_i*x_i - C >= 0. 
%  (D)  MAXIMIZE trace(C,X) SUCH THAT  trace(A_i*X) = f_i, i=1,..,m and X >= 0.
%   
% csdp uses the CSDP library, EPL.
%
%   [y,fvals,exitflag,stats,X] = csdp(f,A,b,lb,ub,sdcone,y0,opts)
%
%   Input arguments:
%       f - linear objective vector 
%       A - linear constraint matrix (sparse)
%       b - linear constraint rhs (A*y <= b)
%       lb - decision variable lower bounds
%       ub - decision variable upper bounds (lb <= y <= ub)
%       sdcone - Semidefinite Cone constraints (see below)
%       y0 - initial solution guess
%       opts - solver options (see below)
%
%   Return arguments:
%       y - solution vector
%       fval - structure with pval (primal) and dval (dual) objective values at the solution
%       exitflag - exit status (see below)
%       stats - solution status structure
%       X - cell array of dual solution vectors of the form [LB;UB;LP;SDP]
%
%   Semidefinte Cones:
%       sdcone is a cell array where each cell is a semidefinite cone. Each 
%       cone is a sparse matrix of the form:
%           sdcone{1} = [C A0 A1 A2]
%       where each matrix is assumed to be square and symmetric. For 
%       example, if C and A0 are 2x2 symmetric matrices in a problem with 
%       one variable, then supplied to csdp would be:
%           sdcone{1} = sparse([C(:) A0(:)]);
%
%   options Fields:
%       maxiter - maximum solver iterations
%       maxtime - maximum solver execution time
%       display - 0 to 2 display level
%       For the remainder see csdpset.
%
%   Return Status:
%       0 - csdp converged
%       1 - primal infeasible
%       2 - dual infeasible
%       3 - partial success (full accuracy not achieved)
%       4 - maximum iterations reached
%       5 - stuck at edge of primal feasbility (failure)
%       6 - stuck at edge of dual infeasiblity (failure)
%       7 - lack of progress (failure)
%       8 - X, Z, or O was singular (failure)
%       9 - Detected NaN or Inf (failure)
%      10 - General easy_sdp failure
%      11 - Failed check on C (check symmetry)
%      12 - Failed constraint check
%     -27 - maximum time reached
%     -50 - user termination
%
%   Copyright (C) 2013 Jonathan Currie (I2C2)