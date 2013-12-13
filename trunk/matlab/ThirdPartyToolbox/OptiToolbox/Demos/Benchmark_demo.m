%% OPTI Toolbox Benchmarking Demo
%
% This file demonstrates how to run optimization benchmarks using OPTI 
% Toolbox. Note you should be familiar with the operation of OPTI Toolbox
% by reading the accompanying examples (e.g. Basic_demo.m)
%
%   Copyright (C) 2011 Jonathan Currie (I2C2)

%% Determing Available Solvers
% OPTI Toolbox comes with a number of solvers, thus to determine which
% ones are available on your system for benchmarking, you can type:

checkSolver()

%% Example 1 - Automatic Benchmarking
% Running optiBench(problem type) will automatically run all available test
% problems for a given problem type, across all available solvers which can
% solve that type of problem. Once finished, it will print individual
% problem results, as well as generating a comparison plot.

optiBench('LP');

%% Example 2 - Customizing Automatic Benchmarking
% By suppling a second argument to optiBench, you can control the number of
% test problems to run.

optiBench('MILP',5);

%% Example 3 - Automatic NLP Benchmarking
% Supplied are the first 50 Hock-Schittkowski test NLPs for benchmarking
% NLP solvers. You can specify how many of these to run over each solver.

optiBench('NLP',10);

%% Example 4 - Manual Benchmarking
% You can also manually setup a benchmark test for an individual solver.
% Using test_solver(problem type, solver) you can specify a problem type
% and a solver to test, returning a vector of execution times and a result
% structure detailing the results of each solve.

[times,results] = testSolver('LP','clp') 

%% Example 5 - Customizing Manual Benchmarking
% As above, you can specify the number of problems to test the solver
% against. If you specify more than available, it will return a warning and
% execute the maximum problems available.

[times,results] = testSolver('QP','ooqp',5)

%% Example 6 - Determing the Best Solver for a OPTI Problem
% You can also supply an OPTI prolbem which will be run against each solver
% capable of solving it, to determine the best solver for your problem.

%Get a sample LP problem
Opt = opti(lp_prob(8));

%Test it against all LP solvers
[times,results] = testProblem(Opt)

%% Example 7 - Returning Pre-Generated Problems
% The lowest level available is the individual problem level, where you can
% obtain pre-written examples for each problem type.

%Get a LP problem
OptLP = opti(lp_prob(1))
%Get a MILP problem
OptMILP = opti(milp_prob(1))
%Get a QP problem
OptQLP = opti(qp_prob(1))
%Get a NLP problem
OptNLP = opti(nlp_prob(1))

