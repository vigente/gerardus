%% OPTI Toolbox File IO Demo
%
% This file loads a number of supplied examples and shows how to solve
% them using the OPTI Toolbox. You should read and complete Basic_demo.m 
% BEFORE running the below examples.
%
% The underlying parser uses CoinUtils + GLPK for reading and writing all
% file types.
%
%   Copyright (C) 2011 Jonathan Currie (I2C2)

%% Loading a MPS Problem
% OPTI Toolbox is supplied with a number of example LP, MILP & QP problems
% in MPS, QPS, MOD and LP format, located in Test Problems/. To
% load a MPS problem, simply use the command below. Returned will be a 
% optiprob structure containing the data in the MPS file. Note all matrices 
% are returned as sparse matrices. If you do not specify a file extension 
% you must supply the  file type via a second argument.

prob = coinRead('testLP.mps')

%% Example 1 - Solving a Loaded MPS Problem
% Solving a loaded MPS problem is simple, just pass it to the opti
% constructor and call solve:

Opt = opti(prob) 

x = solve(Opt)

%% Loading a QPS Problem
% You can also load a problem by specifying the file type explicitly, as
% well as printing the headers + errors as they exist:

prob = coinRead('testQP2','qps',1)

%% Example 2 - Solving a Loaded QPS Problem
% Solving a loaded QPS problem is the same as an MPS problem:

Opt = opti(prob,optiset('solver','clp','maxiter',5e3)) 

x = solve(Opt)

%% Loading a GMPL Problem
% GMPL is a subset of AMPL used by GLPK, and has been implemented into the
% reader:

prob = coinRead('color','gmpl')

%% Example 3 - Solving a Loaded GMPL Problem
% As above, this is just as easy as any other file type:

Opt = opti(prob,optiset('solver','glpk')) 

[x,fval,ef,info] = solve(Opt);
info

%% Loading a LP Problem
% The LP problems compatible with this command are the CPLEX type problems,
% however using only a subset of the operators.

prob = coinRead('location','lp')

%% Example 4 - Solving a Loaded LP Problem
% As above, this is just as easy as any other file type:

Opt = opti(prob,optiset('solver','glpk')) 

[x,fval,ef,info] = solve(Opt);
info

%% Example 5 - Writing a MPS Problem
% You can also write an opti problem to an MPS file. The file will be 
% placed in the current Matlab directory.

write(Opt,'testMPS','mps')

%% Example 6 - Writing a LP Problem
% You can also write an opti problem to an LP file:

write(Opt,'testLP','lp')

%% Example 7 - Comparing Results!
% Check written problems vs original result:
probmps = coinRead('testMPS.mps');
problp = coinRead('testLP.lp');

[x,fval] = solve(opti(prob));
[x1,fval1] = solve(opti(probmps));
[x2,fval2] = solve(opti(problp));

%Below should all be the same!
fval
fval1
fval2

%% Remove written files (just a clean up - demo over!)
delete('testMPS.mps');
delete('testLP.lp');

