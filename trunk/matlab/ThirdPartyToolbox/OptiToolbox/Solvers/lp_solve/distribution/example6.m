f = [110*1.3 30*2.0 125*1.56 75*1.8 95*.95 100*2.25 50*1.35];
A = [120 210 150.75 115 186 140 85;
     110 30 125 75 95 100 50;
     1 1 1 1 1 1 1;
     1 -1 0 0 0 0 0;
     0 0 1 0 -2 0 0;
     0 0 0 -1 0 -1 1];
b = [55000;40000;400;0;0;0];
lp = lp_maker(f, A, b, [-1; -1; -1; -1; -1; -1], [10 10 10 10 20 20 20], [100 Inf 50 Inf Inf 250 Inf], [], 1, 0);
solvestat = mxlpsolve('solve', lp)
format bank
obj = mxlpsolve('get_objective', lp)
format short
x = mxlpsolve('get_variables', lp)
mxlpsolve('delete_lp', lp);