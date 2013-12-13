%% OPTI Toolbox Optimization Toolbox Overloads Demo
%
% This file illustrates the simplest use of OPTI toolbox where a user
% familiar with the Matlab Optimization Toolbox can use OPTI overloads.
%
% In this demo we will cover the following functions:
%
%   - linprog       (Linear Programming)
%   - bintprog      (Binary Integer Programming)
%   - mintprog      (Mixed Integer Linear Programming)
%   - quadprog      (Quadratic Programming)
%   - lsqcurvefit   (Nonlinear Least Squares)
%   - fsolve        (System of Nonlinear Equations)
%   - fminunc       (Unconstrained Nonlinear Optimization)
%   - fmincon       (Constrained Nonlinear Optimization)
%
% As well as setting up OPTI solver options which can be used with the
% overloaded functions.
%
%   Copyright (C) 2011 Jonathan Currie (I2C2)

%% Problem 1 - Linear Programming
% This is a simple two decision variable LP which we use the OPTI version
% of LINPROG to solve.

%Problem
f = -[6 5]';                %Objective Function (min f'x)
A = [1,4; 6,4; 2, -5];      %Linear Inequality Constraints (Ax <= b)
b = [16;28;6];    
lb = [0;0];                 %Bounds on x (lb <= x <= ub)
ub = [10;10];

%% Example 1 - opti_linprog()
% Solved using an OPTI LP solver. Note the function prototype is identical
% to that of the Optimization Toolbox function, but with 'opti_' in front
% of it:

[x,fval,exitflag,info] = opti_linprog(f,A,b,[],[],lb,ub)

%% Problem 2 - Binary Integer Programming
% Another simple two decision variable LP, this time with binary (0,1)
% decision variables. We will use the OPTI version of BINTPROG to solve
% this problem.

f = -[6 5]';
A = [-3,5; 6,4; 3, -5; -6, -4]; 
b = [6;9;1;3];

%% Example 2 - opti_bintprog()
% Solved using an OPTI BILP/MILP solver.

[x,fval,exitflag,info] = opti_bintprog(f,A,b)

%% Problem 3 - Mixed Integer Linear Programming
% A four decision variable problem with mixed continuous and integer
% decision variables. While an Optimization Toolbox overload does not
% exist, an OPTI version called MINTPROG is implemented.

f = [2, 3, 7, 7];
A = [-1, -1, 2, 5;1, -2, -1, -4];
b = [-2, -3]';
lb = zeros(4,1); 
ub = [30 100 20 1]';  
xint = 'CICI';

%% Example 3 - opti_minprog()
% Solved using an OPTI MILP solver.

[x,fval,exitflag,info] = opti_mintprog(f,A,b,[],[],lb,ub,xint)

%% Problem 4 - Quadratic Programming
% A three decision variable QP. We will use the OPTI version of QUADPROG to
% solve this problem:

H = eye(3);
f = -[2 3 1]';
A = [1 1 1;3 -2 -3; 1 -3 2]; 
b = [1;1;1];

%% Example 4 - opti_quadprog()
% Solved using an OPTI QP solver.

[x,fval,exitflag,info] = opti_quadprog(H,f,A,b)

%% Problem 5 - Nonlinear Least Squares
% Solving a nonlinear least squares (data fitting) problem with two
% decision variables. We will use the OPTI version of LSQCURVEFIT to solve
% this problem:

fun = @(x,xdata) x(1)*exp(x(2)*xdata);
x0 = [100; -1];
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];

%% Example 5 - opti_lsqcurvefit()
% Solved using an OPTI NLS solver.

[x,fval,exitflag,info] = opti_lsqcurvefit(fun,x0,xdata,ydata)

%% Problem 6 - Nonlinear Equation Solver
% Solving a system of nonlinear equations with two decision variables. We
% will use the OPTI version of FSOLVE to solve this problem:

fun = @(x) [2*x(1) - x(2) - exp(-x(1));
             -x(1) + 2*x(2) - exp(-x(2))];
x0 = [-5;5];

%% Example 6 - opti_fsolve()
% Solved using an OPTI SNLE solver.

[x,fval,exitflag,info] = opti_fsolve(fun,x0)

%% Problem 7 - Unconstrained Nonlinear Optimization
% Solving an unconstrained nonlinear optimization problem with two decision
% variables. We will use the OPTI version of FMINUNC to solve this problem:

fun = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
x0 = [0 0]';

%% Example 7 - opti_fminunc()
% Solved using an OPTI UNO solver.

[x,fval,exitflag,info] = opti_fminunc(fun,x0)

%% Problem 8 - Constrained Nonlinear Optimization
% Solving a box constrained nonlinear optimization problem with two
% decision variables. We will use the OPTI version of FMINCON to solve this
% problem:

fun = @(x) sin(x(1) + x(2)) + (x(1) - x(2))^2 - 1.5*x(1) + 2.5*x(2) + 1;    
lb = [-1.5;-3];
ub = [4;3];
x0 = [0;0];

%% Example 8 - opti_fmincon()
% Solved using an OPTI NLP solver.

[x,fval,exitflag,info] = opti_fmincon(fun,x0,[],[],[],[],lb,ub)

%% Adding Options
% You may also add OPTI options to the overloaded function. For example to
% choose an alternative NLP solver other than the default, examine the
% other solvers available:

checkSolver('NLP')

% Then select it via optiset:

opts = optiset('solver','lbfgsb')

%% Example 9 - NLP with L-BFGS-B
% Now the previous bounded problem can be solved with the L-BFGS-B solver
% by passing 'opts' to the routine:

[x,fval,exitflag,info] = opti_fmincon(fun,x0,[],[],[],[],lb,ub,[],opts)

%% Displaying solver information
% You can also display iteration by iteration information from some
% solvers by setting the 'display' field:

opts = optiset(opts,'display','iter')

%% Example 10 - Iteration Printing
% With the new options, the objective value can be printed at each
% iteration:

x = opti_fmincon(fun,x0,[],[],[],[],lb,ub,[],opts)

