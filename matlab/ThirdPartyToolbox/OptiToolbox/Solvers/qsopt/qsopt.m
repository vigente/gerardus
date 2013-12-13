% QSOPT  Solve a LP using QSOPT
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_qsopt() INSTEAD!
%
% qsopt uses the QSOPT library.
%
%   [x,fval,exitflag,drow,dcol] = qsopt(f, A, b, lb, ub, opts)
%
%   Input arguments:
%       f - linear objective vector
%       A - linear inequality + equality matrix (sparse)
%       b - linear constraint rhs
%       lb - decision variable lower bounds
%       ub - decision variable upper bounds
%       opts - solver options (see below)
%
%   Return arguments:
%       x - solution vector
%       fval - objective value at the solution
%       exitflag - exit status (see below)
%       drow - dual row solution
%       dcol - dual col solution
%
%   Option Fields (all required):
%       maxiter - maximum solver iterations
%       maxtime - maximum execution time [s]
%       display - solver display level [0,1]
%       nin - number of inequality constraints (determines rows in A)
%
%   Return Status:
%       1 - optimal
%       0 - maximum iterations / maximum time exceeded
%      -1 - infeasible
%      -2 - inaccuracy / unbounded / other
%
%
%   Copyright (C) 2011 Jonathan Currie (I2C2)
%
%   See Also opti_qsopt