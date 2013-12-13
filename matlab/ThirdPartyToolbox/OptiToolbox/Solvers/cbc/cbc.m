% CBC  Solve a MILP using CBC
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_cbc() INSTEAD!
%
% cbc uses the Coin-Or Branch and Cut library.
%
%   [x,fval,exitflag,iter] = cbc(f, A, rl, ru, lb, ub, xtype, sos, opts)
%
%   Input arguments:
%       f - linear objective vector
%       A - linear constraint matrix (sparse)
%       rl - linear constraint lhs
%       ru - linear constraint rhs
%       lb - decision variable lower bounds
%       ub - decision variable upper bounds
%       xtype - decision variable integrality ('C', 'I' or 'B')
%       sos - SOS structure with fields type, index and weight
%       opts - solver options (see below)
%
%   Return arguments:
%       x - solution vector
%       fval - objective value at the solution
%       exitflag - exit status (see below)
%       iter - number of iterations taken by the solver
%
%   Option Fields (all optional):
%       tolint - integer tolerance
%       maxnodes - maximum nodes to explore
%       maxtime - maximum execution time [s]
%       display - solver display level [0,1]
%
%   Return Status:
%       1 - looks optimal
%       0 - maximum nodes / time exceeded, or gap cannot be reached
%      -1 - lp relaxation infeasible
%      -2 - inaccuracy / unbounded / other
%      -5 - user exit
%
%
%   Copyright (C) 2011 Jonathan Currie (I2C2)
%
%   See Also opti_cbc