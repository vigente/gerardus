%% Few optimization test problems for SymBuilder
clear all

%Unfinished code, but sort of working! Requires the Symbolic Toolbox.

% Attempts to generate symbolic derivatives for optimization.

%% LP1 [fval = -31.4]
clc
B = SymBuilder();

B.AddObj('-6*x1 -5*x2');
B.AddCon('x1 + 4*x2 <= 16');
B.AddCon('6*x1 + 4*x2 <= 28');
B.AddCon('2*x1 - 5*x2 <= 6');
B.AddBound('0 <= x <= 10');

Build(B)

[x,fval,ef,info] = Solve(B)

%% MILP1 [fval = -15]
clc
B = SymBuilder();

B.AddObj('3*x1 -7*x2 -12*x3');
B.AddCon('-3*x1 + 6*x2 + 8*x3 <= 12');
B.AddCon('6*x1 - 3*x2 + 7*x3 <= 8');
B.AddCon('-6*x1 + 3*x2 + 3*x3 <= 5');
B.AddInteger('x = I');

Build(B)

[x,fval,ef,info] = Solve(B)

%% QP1 [fval = -6.4138]
clc
B = SymBuilder();

B.AddObj('x1^2 -2*x1*x2 + 2*x2^2 - 2*x1 - 6*x2');
B.AddCon('x1 + x2 <= 2');
B.AddCon('-x1 + 2*x2 <= 2');
B.AddCon('2*x1 + x2 <= 3');
B.AddCon('x1 + 1.5*x2 = 2');
B.AddBound('0 <= x <= 10');

Build(B)

[x,fval,ef,info] = Solve(B)

%% MIQP1 [fval = -2.75]
clc
B = SymBuilder();

B.AddObj('x1^2 + x2^2 + x3^2 - 2*x1 - 3*x2 - x3');
B.AddCon('x1 + x2 + x3 <= 1');
B.AddCon('3*x1 - 2*x2 - 3*x3 <= 1');
B.AddCon('x1 - 3*x2 + 2*x2 <= 1');
B.AddInteger('x2 = I');

Build(B)

[x,fval,ef,info] = Solve(B)

%% UNO [fval = 0]
clc
B = SymBuilder();

B.AddObj('(1-x1)^2 + 100 * (x2-x1^2)^2');

Build(B)

[x,fval,ef,info] = Solve(B,[0;0])

%% NLP1 [fval = 0]
clc
B = SymBuilder();

B.AddObj('(x1-x2)^2 + (x2+x3-2)^2 + (x4-1)^2 + (x5-1)^2');
B.AddCon('x1 + 3*x2 = 4');
B.AddCon('x3 + x4 - 2*x5 = 0');
B.AddCon('x2 - x5 = 0');

Build(B)

x0 = [ 2.5 0.5 2 -1 0.5 ];
[x,fval,ef,info] = Solve(B,x0)


%% MINLP1 [fval = -0.72007]
clc
clear
B = SymBuilder();

B.AddObj('sin(pi*x1/12)*cos(pi*x2/16)');
B.AddCon('-x1 + 2.5*x2 <= 1');
B.AddCon('x1 + 2.5*x2 <= -15');
B.AddInteger('x = I');

Build(B)

x0 = [0;0];
[x,fval,ef,info] = Solve(B,x0)

%% Changing Solver 
clc
B = SymBuilder();

B.AddObj('-6*x1 -5*x2');
B.AddCon('x1 + 4*x2 <= 16');
B.AddCon('6*x1 + 4*x2 <= 28');
B.AddCon('2*x1 - 5*x2 <= 6');
B.AddBound('0 <= x <= 10');

Build(B)
opts = symbset('solver','scip');

[x,fval,ef,info] = Solve(B,[],opts)

%% Approximate 2nd Derivatives
clc
B = SymBuilder();

B.AddObj('(x1-x2)^2 + (x2+x3-2)^2 + (x4-1)^2 + (x5-1)^2');
B.AddCon('x1 + 3*x2 = 4');
B.AddCon('x3 + x4 - 2*x5 = 0');
B.AddCon('x2 - x5 = 0');
B.AddCon('x2^2 + x5 = 4');
B.AddInteger('x3 = I');

Build(B)

x0 = [ 2.5 0.5 2 -1 0.5 ];
opts = symbset('use2ndDerivs','no');
[x,fval,ef,info] = Solve(B,x0,opts)