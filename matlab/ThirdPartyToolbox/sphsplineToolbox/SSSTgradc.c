/*--------------------------------------------------------------------
 *	$Id: SSSTgradc.c,v 1.1.1.1 2008/05/09 21:34:52 myself Exp $
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
 * C versions of the spherical spline gradient functions described in the paper
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

/* gradspline2d_Parker computes the gradient of the Green function for a 2-d surface spline
 * as per Parker [1994], G(x) = dilog(),
 * where x is cosine of distances. All x must be -1 <= x <= +1.
 * Parameters passed are:
 * par[0] = 6/M_PI^2 (to normalize results)
 */
 
double gradspline2d_Parker (double x, double par[])
{	/* Not normalized to 0-1 */
	if (x == +1.0 || x == -1.0) return (0.0);
	return (log(0.5 - 0.5 * x) * sqrt ((1.0 - x) / (1.0 + x)));
}

/* gradspline2d_Wessel_Becker computes the gradient of the Green function for a 2-d surface spline
 * in tension as per Wessel and Becker [2007], G(x) = M_PI * Pv(-x)/sin (v*x) - log (1-x),
 * where x is cosine of distances.  All x must be -1 <= x <= +1.
 * Parameters passed are:
 * par[0] = Real(nu)
 * par[1] = Imag(nu)
 * par[2] = G(-1)
 * par[3] = G(+1)
 * par[4] = Real {sin (nu * M_PI}
 * par[5] = Imag {sin (nu * M_PI)} == 0
 * par[6] = 1 / (par[3] - par[2])
 */
 
#define SLOPPY 1.0e-6

double gradspline2d_Wessel_Becker (double x, double par[])
{	/* g = -M_PI * (v+1)*[x*Pv(-x)+Pv+1(-x)]/(sin (v*x)*sin(theta)) - sqrt ((1.0 + x) / (1.0 - x)), normalized to 0-1 */
	int n;
	double z[2], v1[2], pq[4], s;
	
	if (fabs(x) > (1.0 - SLOPPY)) return (0.0);

	GMT_PvQv (-x, par, pq, &n);			/* Get P_nu(-x) */
	z[0] = pq[0] * x;	z[1] = pq[1] * x;	/* Get x*P_nu(-x) */
	v1[0] = par[0] + 1.0;	v1[1] = par[1];		/* Get nu+1 */
	GMT_PvQv (-x, v1, pq, &n);			/* Get P_(nu+1)(-x) */
	z[0] += pq[0];	z[1] += pq[1];			/* Get x*P_nu(-x) + P_(nu+1)(-x) */
	Cdiv (z, &par[4], pq);				/* Get ---"--- / sin (nu*M_PI) */
	Cmul (pq, v1, z);				/* Mul by nu + 1 */
	s = M_PI / sqrt (1.0 - x*x);			/* Mul by pi/sin(theta) */
	z[0] *= s;
	z[0] += sqrt ((1.0 + x)/(1.0 - x));		/* Add in last term */
	
	return (-z[0]);
}

/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	struct GRD_HEADER grd;
	double p = 0.0, *z, *x, *y, par[7];
	char *argv = "SSSTgradc";
	int error, nx, ny, n, i;
 
	if ((nrhs < 1 || nrhs > 2) || nlhs < 1) {
		mexPrintf ("usage: z = SSSTgradc(x[,p]);\n");
		return;
	}

	if (nrhs == 2) {
		nx  = mxGetN (prhs[1]);
		ny  = mxGetM (prhs[1]);
		if (nx > 1 || ny > 1) {
			mexPrintf ("SSSTgradc: p must be a constant\n");
			return;
		}
		y = mxGetPr (prhs[1]);
		p = y[0];
	}
	init_spline (p, par);
	x = mxGetPr (prhs[0]);
	nx  = mxGetN (prhs[0]);
	ny  = mxGetM (prhs[0]);
	n = nx * ny;

	/* Create a matrix for the return array */

	plhs[0] = mxCreateDoubleMatrix (ny, nx, mxREAL);
    
	z = mxGetPr (plhs[0]);

	if (p == 0.0)
 		for (i = 0; i < n; i++) z[i] = gradspline2d_Parker (x[i], par);
	else
 		for (i = 0; i < n; i++) z[i] = gradspline2d_Wessel_Becker (x[i], par);

	return;
}
