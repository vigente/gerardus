/*--------------------------------------------------------------------
 *	$Id: SSST.h,v 1.1.1.1 2008/05/09 21:34:52 myself Exp $
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

#include "gmt.h"
#include "mex.h"

#ifndef M_LOG_2
#define M_LOG_2 0.69314718055994530942
#endif
#ifndef M_GAMMA
#define M_GAMMA 0.577215664901532860606512
#endif

/* Functions for complex math */

static void Cdiv (double A[], double B[], double C[])
{	/* Complex division */
	double i_denom;
	i_denom = 1.0 / (B[0]*B[0] + B[1]*B[1]);
	C[0] = (A[0]*B[0] + A[1]*B[1]) * i_denom;
	C[1] = (A[1]*B[0] - A[0]*B[1]) * i_denom;
}

static void Cmul (double A[], double B[], double C[])
{	/* Complex multiplication */
	C[0] = A[0]*B[0] - A[1]*B[1];
	C[1] = A[0]*B[1] + A[1]*B[0];
}

static void Ccot (double Z[], double cotZ[])
{	/* Complex cot(z) */
	double sx, cx, e, A[2], B[2];
	
	sincos (2.0*Z[0], &sx, &cx);
	e = exp (-2.0*Z[1]);
	A[0] = -e * sx;		A[1] = B[0] = e * cx;
	A[1] += 1.0;	B[0] -= 1.0;	B[1] = -A[0];
	Cdiv (A, B, cotZ);
}

void init_spline (double p, double *par) {
	if (p == 0.0)	/* No tension: Set up required parameters */
		par[0] = 6.0 / (M_PI*M_PI);
	else {	/* General tension case */
		par[0] = -0.5;
		if (p <= 0.5) {	/* nu is real */
			double z[2];
			par[0] += sqrt (0.25 - p * p);
			par[1] = 0.0;
			par[4] = sin (M_PI * par[0]);
			z[0] = par[0] + 1.0;	z[1] = 0.0;
			par[3] = (M_PI / tan (M_PI * par[0])) - M_LOG_2 + 2.0 * (M_GAMMA + GMT_psi (z, NULL));
		}
		else {	/* nu is complex */
			double z[2], cot_piv[2], psi[2];
			par[1] = sqrt (p * p - 0.25);
			par[4] = -cosh (M_PI * par[1]);
			z[0] = par[0] * M_PI;	z[1] = par[1] * M_PI;
			Ccot (z, cot_piv);
			cot_piv[0] *= M_PI;	cot_piv[1] *= M_PI;
			z[0] = par[0] + 1.0;	z[1] = par[1];
			(void) GMT_psi (z, psi);
			psi[0] += M_GAMMA;
			psi[0] *= 2.0;	psi[1] *= 2.0;
			z[0] = cot_piv[0] + psi[0] - M_LOG_2;	z[1] = cot_piv[1] + psi[1];
			par[3] = z[0];	/* Ignore complex parts which cancels out */
		}
		par[2] = (M_PI / par[4]) - M_LOG_2;
		par[6] = 1.0 / (par[3] - par[2]);
		par[7] = p;
	}
}
