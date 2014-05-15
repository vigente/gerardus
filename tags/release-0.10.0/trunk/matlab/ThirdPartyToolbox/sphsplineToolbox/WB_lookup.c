/*--------------------------------------------------------------------
 *	$Id: WB_lookup.c,v 1.1.1.1 2008/05/09 21:34:52 myself Exp $
 *
 *	Copyright (c) 2008 by P. Wessel
 *	See COPYING file for copying and redistribution conditions.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation; version 2 of the License.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	Contact info: gmt.soest.hawaii.edu
 *--------------------------------------------------------------------*/
/*
 * C versions of the spherical spline functions described in the paper
 * Wessel, P., and J. M. Becker, 2008, Interpolation using a generalized
 * Green's function for a spherical surface spline in tension, G. J. Int.,
 * in press.
 *
 * The generation of suitable mex files for your system will require:
 * 1. You have matlab
 * 2. You have installed GMT 4.2 or later
 *
 * Author:      Paul Wessel
 * Date:        10-APR-2008
 *
 */
#include "SSST.h"

/* 2nd implementation of the surface spline in tension using a look-up table
 * of pre-computed values.
 */
 
double WB_lookup (double x, double par[], double *y, double xmin, double i_xinc)
{
	int k;
	double f, f0, df;

	f = (x - xmin) * i_xinc;	/* Floating point index */
	f0 = floor (f);
	df = f - f0;
	k = (int)f0;
	if (df == 0.0) return (y[k]);
	return (y[k]*(1.0 - df) + y[k+1] * df);
}

 /* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	struct GRD_HEADER grd;
	double p = 0.0, *z, *x, *y, *xm, *idx, *c, par[10], xmin, i_xinc;
	char *argv = "WB_lookup";
	int error, nx, ny, n, i;
 
	if (nrhs != 5 || nlhs < 1) {
		mexPrintf ("usage: z = WB_lookup(x,p,y,xmin,idx);\n");
		return;
	}

	nx  = mxGetN (prhs[1]);
	ny  = mxGetM (prhs[1]);
	if (nx > 1 || ny > 1) {
		mexPrintf ("WB_lookup: p must be a constant\n");
		return;
	}
	y = mxGetPr (prhs[1]);
	p = y[0];
	xm = mxGetPr (prhs[3]);
	xmin = xm[0];
	idx = mxGetPr (prhs[4]);
	i_xinc = idx[0];
	init_spline (p, par);
	
	c = mxGetPr (prhs[2]);

	x = mxGetPr (prhs[0]);
	nx  = mxGetN (prhs[0]);
	ny  = mxGetM (prhs[0]);
	n = nx * ny;

	/* Create a matrix for the return array */

	plhs[0] = mxCreateDoubleMatrix (ny, nx, mxREAL);
    
	z = mxGetPr (plhs[0]);

 	for (i = 0; i < n; i++) z[i] = WB_lookup (x[i], par, c, xmin, i_xinc);

	return;
}
