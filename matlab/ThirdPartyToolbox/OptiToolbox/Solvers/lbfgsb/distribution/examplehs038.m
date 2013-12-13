% Test the "lbfgs" Matlab interface on the Hock & Schittkowski test problem
% #38. See: Willi Hock and Klaus Schittkowski. (1981) Test Examples for
% Nonlinear Programming Codes. Lecture Notes in Economics and Mathematical
% Systems Vol. 187, Springer-Verlag.

% The starting point.
x0  = [-3  -1  -3  -1];   % The starting point.
lb  = [-10 -10 -10 -10];  % Lower bound on the variables.
ub  = [+10 +10 +10 +10];  % Upper bound on the variables.

x = lbfgsb(x0,lb,ub,'computeObjectiveHS038','computeGradientHS038',...
           [],'genericcallback','maxiter',80,'m',4,'factr',1e-12,...
           'pgtol',1e-5);
