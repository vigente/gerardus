%LP_SOLVE  Solves mixed integer linear programming problems.
%
%  SYNOPSIS: [obj,x,duals,stat] = lp_solve(f,a,b,e,vlb,vub,xint,scalemode,keep)
%
%     solves the MILP problem
%
%             max v = f'*x
%               a*x <> b
%                 vlb <= x <= vub
%                 x(int) are integer
%
%  ARGUMENTS: The first four arguments are required:
%
%           f: n vector of coefficients for a linear objective function.
%           a: m by n matrix representing linear constraints.
%           b: m vector of right sides for the inequality constraints.
%           e: m vector that determines the sense of the inequalities:
%                     e(i) = -1  ==> Less Than
%                     e(i) =  0  ==> Equals
%                     e(i) =  1  ==> Greater Than
%         vlb: n vector of lower bounds. If empty or omitted,
%              then the lower bounds are set to zero.
%         vub: n vector of upper bounds. May be omitted or empty.
%        xint: vector of integer variables. May be omitted or empty.
%   scalemode: scale flag. Off when 0 or omitted.
%        keep: Flag for keeping the lp problem after it's been solved.
%              If omitted, the lp will be deleted when solved.
%
%  OUTPUT: A nonempty output is returned if a solution is found:
%
%         obj: Optimal value of the objective function.
%           x: Optimal value of the decision variables.
%       duals: solution of the dual problem.

% JONNY EDITs (below + mex): (THIS FILE IS NOW REDUNDANT)
% - Setup to include binary variables explicitly
% - Fixed sparse matrix handling on 64bit builds
% - Added returning of iter, exitflag + status message

function [obj, x, duals , iter, exitflag, statmsg] = lp_solve(f, a, b, e, vlb, vub, xint, scalemode, keep)

if nargin == 0
        help lp_solve;
        return;
end

[m,n] = size(a);
lp = mxlpsolve('make_lp', m, n);
mxlpsolve('set_verbose', lp, 10);
mxlpsolve('set_mat', lp, a);
mxlpsolve('set_rh_vec', lp, b);
mxlpsolve('set_obj_fn', lp, f);
mxlpsolve('set_minim', lp); % default is solving minimum lp.

for i = 1:length(e)
  if e(i) < 0
        con_type = 1;
  elseif e(i) == 0
        con_type = 3;
  else
        con_type = 2;
  end
  mxlpsolve('set_constr_type', lp, i, con_type);
end

if nargin > 4
  for i = 1:length(vlb)
    mxlpsolve('set_lowbo', lp, i, vlb(i));
  end
end

if nargin > 5
  for i = 1:length(vub)
    mxlpsolve('set_upbo', lp, i, vub(i));
  end
end

%JONNY EDIT
if nargin > 6
  for i = 1:length(xint)
      if(xint(i) == 'B')
          mxlpsolve('set_binary', lp, i, 1);
      elseif(xint(i) == 'I')
          mxlpsolve('set_int', lp, i, 1);
      end
  end
end

if nargin > 7
  if scalemode ~= 0
    mxlpsolve('set_scaling', lp, scalemode);
  end
end

result=mxlpsolve('solve', lp);
if result == 0 | result == 1 | result == 11 | result == 12
%  [obj, x, duals, stat] = mxlpsolve('get_solution', lp), result;
  [obj, x, duals] = mxlpsolve('get_solution', lp);
else
  obj = [];
  x = [];
  duals = [];  
end

%JONNY EDIT
iter = mxlpsolve('get_total_iter', lp);
exitflag = result;
statmsg = mxlpsolve('get_statustext', lp, result);

if nargin < 9
  mxlpsolve('delete_lp', lp);
end
