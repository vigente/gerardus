%% OPTI Toolbox Differentiation Demo
%
% This file demonstrates each of the different differentiation algorithms
% supplied with the OPTI Toolbox. Note you should be familiar with the 
% operation of OPTI Toolbox by reading the accompanying examples, as well 
% as you should have completed the NLP demo.
%
%   Copyright (C) 2011 Jonathan Currie (I2C2)

%% Gradient vs Jacobian
% The terms Gradient and Jacobian are associated with differentiation
% functions, and for the purposes of this toolbox they are basically
% interchangeable:
%
% Gradient - The first derivative of the objective function (1 x n)
% Jacobian - The first derivative of the constraint function (m x n)
%
% Where you will note the main difference is the gradient is a vector, and
% the Jacobian is a matrix. However this is just a rule of thumb and will
% not hold in all instances (constraint functions with one row for
% example!). All functions are named as returning a Jacobian, which could
% also be a Gradient Vector.
 
%% Problem 1
% For the following examples we are going to use this nonlinear objective
% function:

fun = @(x) sin(pi*x(1)/12) * cos(pi*x(2)/16);

x = [0.75 0.25]';

%% Example 1 - Automatic Differentiation (AD)
% Automatic differentiation is one of the most powerful differentiation
% strategies which can provide error free gradients of a function by
% applying the chain rule to each operation. 
%
% The AD routines implemented in OPTI are supplied by adiff, a Matlab
% project by William McIlhagga. While most Matlab functions are overloaded,
% you cannot use external code (i.e. via MEX) or any toolbox or class
% functions.

dx = autoJac(fun,x)

%% Example 2 - Numerical Differentiation (ND)
% Numerical differentiation using finite differences is a computationally
% expensive procedure which can result in an inaccurate gradient if the
% internal perturbations are not chosen correctly. 
%
% The ND routine implemented in OPTI is the Intel MKL djacobi function
% which approximates the derivative using central differences. This is
% implemented via a MEX function which repeatedly calls the function in
% order to close in on the gradient.

dx = mklJac(fun,x)

%% Example 3 - Symbolic Differentiation (SD)
% Symbolic differentiation analytically differentiates the function as a
% symbolic expression, resulting in a single expression for the gradient.
% Complications occur if the function cannot be analytically differentiated
% or the symbolic routine cannot find a derivative.
%
% The SD routine implemented in OPTI uses the Matlab Symbolic Toolbox as
% well as two wrapper functions in order to generate the gradient. The
% wrapper functions convert the function handle to a symbolic expression
% and vice-versa.

grad = symJac(fun)

if(~isempty(grad)) %don't run if Symbolic Toolbox not installed
    dx = grad(x)
end

%% Problem 2
% For the next few examples we will be solving this NLP:

obj = @(x) log(1+x(1)^2) - x(2);

lb = [-2 -2]';
ub = [2 2]';

%% Example 4 - Applying AD to NLP Solving
% To use AD to generate the objective gradient for IPOPT create a function
% handle which calls autoJac, and then pass this to the optiprob function:

grad = @(x) autoJac(obj,x);

opts = optiset('solver','ipopt');
Opt = opti('obj',obj,'grad',grad,'bounds',lb,ub,'options',opts)

[x,fval,exitflag,info] = solve(Opt,[0;0]);
fval
info

%% Example 5 - Applying ND to NLP Solving
% To use ND to generate the objective gradient for IPOPT just replace the
% above with mklJac:

grad = @(x) mklJac(obj,x);

Opt = opti('obj',obj,'grad',grad,'bounds',lb,ub,'options',opts)

[x,fval,exitflag,info] = solve(Opt,[0;0]);
fval
info

%% Example 6 - Applying SD to NLP Solving
% To use SD we can generate the analytical gradient and use this as our
% gradient function:

grad = symJac(obj)

if(~isempty(grad))
    Opt = opti('obj',obj,'grad',grad,'bounds',lb,ub,'options',opts)

    [x,fval,exitflag,info] = solve(Opt,[0;0]);
    fval
    info
else
    display('Cannot run example without Symbolic Toolbox!');
end

%% Conclusion
% By default OPTI Toolbox will always use Numerical Differentiation
% (mklJac) for gradients / Jacobians if one has not been provided. This
% provides the best combination of flexibility with complex objective / 
% constraint functions and performance. However if you can use SD to
% generate a gradient this would be a preferred method. 
%
% If you find the optimizer is failing with ND and AD is a suitable 
% candidate you can try it to ensure your gradients are accurate.
