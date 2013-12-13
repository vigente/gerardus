f = [143 60];
A = [120 210; 110 30; 1 1];
b = [15000; 4000; 75];
lp = lp_maker(f, A, b, [-1; -1; -1], [], [], [], 1, 0);
solvestat = mxlpsolve('solve', lp)
format bank
obj = mxlpsolve('get_objective', lp)
format short
x = mxlpsolve('get_variables', lp)
mxlpsolve('delete_lp', lp);
