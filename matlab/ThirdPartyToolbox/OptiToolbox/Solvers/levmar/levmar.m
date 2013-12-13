% LEVMAR  Solve a NLS using LEVMAR
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_levmar() INSTEAD!
%
% levmar uses the C/C++ Implementation of Levenberg-Marquardt library.
%
%   [x,fval,exitflag,iter] = levmar(fun,grad,x0,ydata,lb,ub,A,b,Aeq,beq,opts)
%
%   Input arguments:
%       fun - nonlinear fitting function handle
%       grad - gradient of nonlinear fitting function handle (optional)
%       x0 - initial solution guess
%       ydata - fitting data
%       lb - decision variable lower bounds
%       ub - decision variable upper bounds
%       A - linear inequality matrix (dense)
%       b - linear inequality rhs
%       Aeq - linear equality matrix (dense)
%       beq - linear equality rhs
%       opts - solver options (see below)
%
%   Return arguments:
%       x - solution vector
%       fval - objective value at the solution
%       exitflag - exit status (see below)
%       iter - number of iterations taken by the solver
%
%   Option Fields (all optional):
%       maxiter - maximum solver iterations
%
%   Return Status:
%       1 - optimal / stopped by small grad or Dp
%       0 - maximum iterations exceeded
%      -1 - could not converge / tolerance too small
%      -2 - nan or inf
%      -3 - unknown error
%
%   
%   Based in parts on the original MEX interface by Manolis Lourakis.
%
%   Copyright (C) 2011 Jonathan Currie (I2C2)