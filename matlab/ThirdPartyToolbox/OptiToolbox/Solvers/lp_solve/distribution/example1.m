% lpsolve example1 from the manual.

echo on
lp=mxlpsolve('make_lp', 0, 4);
mxlpsolve('set_verbose', lp, 3);
mxlpsolve('set_obj_fn', lp, [1, 3, 6.24, 0.1]);
mxlpsolve('add_constraint', lp, [0, 78.26, 0, 2.9], 2, 92.3);
mxlpsolve('add_constraint', lp, [0.24, 0, 11.31, 0], 1, 14.8);
mxlpsolve('add_constraint', lp, [12.68, 0, 0.08, 0.9], 2, 4);
mxlpsolve('set_lowbo', lp, 1, 28.6);
mxlpsolve('set_lowbo', lp, 4, 18);
mxlpsolve('set_upbo', lp, 4, 48.98);
mxlpsolve('set_col_name', lp, 1, 'COLONE');
mxlpsolve('set_col_name', lp, 2, 'COLTWO');
mxlpsolve('set_col_name', lp, 3, 'COLTHREE');
mxlpsolve('set_col_name', lp, 4, 'COLFOUR');
mxlpsolve('set_row_name', lp, 1, 'THISROW');
mxlpsolve('set_row_name', lp, 2, 'THATROW');
mxlpsolve('set_row_name', lp, 3, 'LASTROW');
mxlpsolve('write_lp', lp, 'a.lp');
mxlpsolve('get_mat', lp, 1, 2)
mxlpsolve('solve', lp)
mxlpsolve('get_objective', lp)
mxlpsolve('get_variables', lp)
mxlpsolve('get_constraints', lp)
mxlpsolve('delete_lp', lp);
echo off
