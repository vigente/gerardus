/* Copyright (C) 2011 Stefan Vigerske
 * All Rights Reserved.
 * This file is distributed under the Eclipse Public License.
 */

/* $Id$ */

/* this file provides dummy implementations of the method metis_nodend as expected by the HSL codes if Metis is not available
 * as in Metis, we implement the method in several naming variants to copy with C and fortran naming style conventions
 */

typedef int idxtype;

void METIS_NODEND(int * a, idxtype * b, idxtype * c, int * d, int * e, idxtype * f, idxtype * perm)
{
  perm[0] = -1;
}

void metis_nodend(int * a, idxtype * b, idxtype * c, int * d, int * e, idxtype * f, idxtype * perm)
{
  perm[0] = -1;
}

void metis_nodend_(int * a, idxtype * b, idxtype * c, int * d, int * e, idxtype * f, idxtype * perm)
{
  perm[0] = -1;
}

void metis_nodend__(int * a, idxtype * b, idxtype * c, int * d, int * e, idxtype * f, idxtype * perm)
{
  perm[0] = -1;
}

void METIS_NodeND(int * a, idxtype * b, idxtype * c, int * d, int * e, idxtype * f, idxtype * perm)
{
  perm[0] = -1;
}
