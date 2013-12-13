%*****************************************************************************
% DSDP5:  Dual-Scaling Algorithm for Positive Semidefinite Programming
% Copyright (c) 2004 by
% S. J. Benson, Y. Ye
% Last modified: 20 August 2004
%*****************************************************************************
%
% > [STAT,y] = DSDP() returns a structure STAT with relevant information 
%              concerning the performance of the solver and an approximate
%              dual solution y .  
%
%   The fields of STAT are:
%
%   Objective Value:
%       stype   = 'PDFeasible' if an feasible primal and dual solutions
%                were computed, 'Infeasible' if dual
%                infeasibility was detected, and 'Unbounded' if primal
%                infeasibility is detected.
%       obj    = objective value at solution 
%       pobj   = an approximately optimal objective value to (P)
%       dobj   = an approximately optimal objective value to (D)
%       stopcode = 0: convergence to prescribed accuracy, 
%                 ~0: termination for other reasons
%
%   Characteristics of Solution
%       tracex = if X was returned, this is the trace of it.  This number
%                also corresponds to the minimum penalty parameter that 
%                could solve this problem.  IMPORTANT: For improved 
%                performance, consider using penalty parameter (see DOPTIONS)
%                other than the default.
%       penalty = the penalty parameter used by the solver, which must be 
%                greater than the trace of the primal solution (see above).
%       errors = several error estimates to the solution. (See DERROR)
%       ynorm = the largest element of y (infinity norm).
%       boundy = the bounds placed on the magnitude of each variable y.
%       mu = final barrier parameter.
%       r      = the multiple of the identity matrix added to
%                C-A'(y) in the final solution to make S positive definite.
%                That is, S = C - A'y + r*I.
%       xy, xdy, xmu = values used to compute X.
%
%   Solver Statistics
%       iterates = number of iterations used by the algorithm.
%       pstep = final primal step size.
%       dstep = final dual step size.
%       pnorm = final norm of distance to central path.
%       gaphist = a history of the duality gap.
%       infhist = a history of the dual infeasibility.
%       datanorm = the Frobenius norm of C, A and b,
%
%  See also: DSDP
%
% DSDP5
% Copyright (c) 2004 by
% S. Benson and Y. Ye
% Last modified: October 2004
%*****************************************************************************

