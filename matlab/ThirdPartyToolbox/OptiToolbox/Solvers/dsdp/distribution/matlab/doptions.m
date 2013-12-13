%*****************************************************************************
% DSDP5:  Dual-Scaling Algorithm for Positive Semidefinite Programming
% Copyright (c) 2005 by
% S. J. Benson, Y. Ye
% Last modified: 05 Feb 2005
%*****************************************************************************
% [OPTIONS] = doptions;
%
% This script sets the default options for use in DSDP.
%
%   The OPTIONS structure may contain any of the following fields: 
%   Problem Formulation:  
%       r0 = multiple of the identity matrix added to the initial dual matrix.
%              S0 = C - A'(y0) + r0*I.  If r0 < 0, a dynamic selection will
%              be used.  Default value is usually very high(~1e12), but 
%              smaller numbers may significantly improve performance. IMPORTANT!
%              [default -1 (use dynamic strategy, usually very big)].
%       zbar = an upper bound for the dual solution.
%       penalty  = penalty parameter used to enforce feasibility (>=0)
%              This parameter can significantly infuence performance.
%              It must be larger than the trace of a solution X:
%              See STAT.tracex
%       boundy = a bound for the dual variables y. IMPORTANT!
%
%   Convergence:
%       gaptol = tolerance for duality gap as a fraction of the value of the 
%              objective functions (>=0).
%       maxit = maximum number of iterations allowed (>=0).
%       steptol = tolerance for stopping due to small steps (>=0).
%       pnormtol = terminate only if pnorm less than this number (>0).
%       dual_bound = stop solver when the objective (D) of a
%                    feasible iterate is greater than this value.
%
%   Detecting Unboundedness and Infeasibility:
%       inftol = maximum infeasibility in (D) allowed at solution.
%       infptol = maximum infeasibility in (P) allowed at solution.
%
%   Printing
%       print = k, if want to display result in each k iteration, else = 0.
%       logtime = 1, to profile DSDP subroutines, else = 0
%              Assumes proper compilation flags.
%       cc   = add this constant to the objective value (for display purposes only)
%
%   Other important options:
%       rho    = potential parameter as multiple of dimension.  Larger numbers
%                can make the algorithm go faster: [ >1]
%                IMPORTANT: Try 5 or 8.
%       reuse = reuse Newton matrix this many times (at most) per iteration.  
%               Applications requiring few iterations(<60) should consider 
%               setting this paramter to 0 or 1. 
%       mu0    = initial barrier parameter (>=0). Smaller numbers may improve 
%                performance. [-1 implies a dynamic strategy)]
%
%   Other less important options:
%       dyanamicrho = Use dynamic rho strategy to reduce barrier.
%       lp_barrier = scale the barrier on LP cones by this amount.
%       xmaker = 1 if return y and dy needed to construct X in STAT.xy, STAT.xdy,  STAT.xmu
%       maxtrustradius = maximum trust radius (>=0). Decrease to improve robustness.
%       bigM  = if > 0,make infeasibility a positive variable with a 
%               cost in the objective.
%*****************************************************************************
function [OPTIONS] = doptions();

OPTIONS.r0=-1;
OPTIONS.penalty=1e8;
OPTIONS.boundy=1e7;
OPTIONS.zbar=1e10;
OPTIONS.maxit=300;
OPTIONS.gaptol=0.000001;
OPTIONS.print=10;
OPTIONS.dual_bound=1e10;
OPTIONS.bigM=0;
OPTIONS.inftol=1.0e-10;
OPTIONS.cc=0;
OPTIONS.steptol=1.0e-2;
OPTIONS.mu0=-1;
OPTIONS.logtime=0;
OPTIONS.xmaker=0;
OPTIONS.dloginfo=0;
OPTIONS.dynamicrho=1;
