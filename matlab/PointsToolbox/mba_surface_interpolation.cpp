/**
 * mba_surface_interpolation.cpp
 *
 * MBA_SURFACE_INTERPOLATION  Scattered data Multilevel B-spline interpolation
 *
 * This MEX-function uses the MBA library [1] to compute a Multilevel
 * B-spline interpolated surface from a scattered set of points.
 *
 * ZI = MBA_SURFACE_INTERPOLATION(X, Y, Z, XI, YI)
 *
 *   X, Y, Z are column vectors of the same size with the 3D
 *   coordinates of a scattered set of points.
 *
 *   XI, YI are column vectors of the same size (but can have a
 *   different size from X, Y and Z) with the 2D coordinates of the
 *   locations that we want to interpolate.
 *
 *   ZI is a vector of the same length as XI and YI, with the
 *   interpolated values.
 *
 * To compile this MEX-file in a 64-bit linux architecture, run from
 * your gerardus/matlab directory
 *
 * >> mex -v -largeArrayDims -outdir PointsToolbox/ -f ./engopts.sh PointsToolbox/mba_surface_interpolation.cpp
 *
 * To compile in a 32-bit linux architecture, run (untested)
 *
 * >> mex -v -outdir PointsToolbox/ -f ./engopts.sh PointsToolbox/mba_surface_interpolation.cpp
 *
 * For Windows, Mac or other architectures, the corresponding section
 * in the ../engopts.sh file will need to be edited before it
 * compiles. I cannot test them, so I haven't touched those sections.
 *
 * [1] MBA - Multilevel B-Spline Approximation
 * Library. http://www.sintef.no/Projectweb/Geometry-Toolkits/MBA/
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  *
  * University of Oxford means the Chancellor, Masters and Scholars of
  * the University of Oxford, having an administrative office at
  * Wellington Square, Oxford OX1 2JD, UK. 
  *
  * This file is part of Gerardus.
  *
  * This program is free software: you can redistribute it and/or modify
  * it under the terms of the GNU General Public License as published by
  * the Free Software Foundation, either version 3 of the License, or
  * (at your option) any later version.
  *
  * This program is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  * GNU General Public License for more details. The offer of this
  * program under the terms of the License is subject to the License
  * being interpreted in accordance with English Law and subject to any
  * action against the University of Oxford being under the jurisdiction
  * of the English Courts.
  *
  * You should have received a copy of the GNU General Public License
  * along with this program.  If not, see
  * <http://www.gnu.org/licenses/>.
  */

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

/* mex headers */
#include <math.h>
#include <matrix.h>
#include <mex.h>

/* C++ headers */
#include <iostream>
#include <limits>

// MBA libary
#include <MBA.h>
#include <UCButils.h>
#include <PointAccessUtils.h>

#ifndef MBA_SURFACE_INTERPOLATION_CPP
#define MBA_SURFACE_INTERPOLATION_CPP

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // check number of input and output arguments
  if (nrhs != 5) {
    mexErrMsgTxt("Five input arguments required.");
  }
  else if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments.");
  }

  // check size of input arguments
  mwSize Mx = mxGetM(prhs[0]); // number of scattered points
  mwSize Mxi = mxGetM(prhs[3]); // number of interpolated (grid) points
  if (mxGetN(prhs[0]) != 1) mexErrMsgTxt( "X must have one column." );
  if (mxGetN(prhs[1]) != 1) mexErrMsgTxt( "Y must have one column." );
  if (mxGetN(prhs[2]) != 1) mexErrMsgTxt( "Z must have one column." );
  if (mxGetN(prhs[3]) != 1) mexErrMsgTxt( "XI must have one column." );
  if (mxGetN(prhs[4]) != 1) mexErrMsgTxt( "YI must have one column." );
  if (Mx != mxGetM(prhs[1]) || Mx != mxGetM(prhs[2])) 
    mexErrMsgTxt( "X, Y and Z must have the same number of points (rows)." );
  if (mxGetM(prhs[3]) != mxGetM(prhs[4])) 
    mexErrMsgTxt( "XI and YI must have the same number of points (rows)." );

  // create output vector and pointer to populate it
  plhs[0] = mxCreateDoubleMatrix(Mxi, 1, mxREAL);
  double *zi = mxGetPr(plhs[0]);

  // create pointers to input vectors
  double *x = mxGetPr(prhs[0]);
  double *y = mxGetPr(prhs[1]);
  double *z = mxGetPr(prhs[2]);
  double *xi = mxGetPr(prhs[3]);
  double *yi = mxGetPr(prhs[4]);
  
  // duplicate the input data in vector format so that we can pass it
  // to the MBA library

  typedef std::vector<double> VectorType;

  // scattered points
  boost::shared_ptr<VectorType> xv(new VectorType);
  boost::shared_ptr<VectorType> yv(new VectorType);
  boost::shared_ptr<VectorType> zv(new VectorType);

  double xmin = std::numeric_limits<double>::max();
  double xmax = std::numeric_limits<double>::min();
  double ymin = std::numeric_limits<double>::max();
  double ymax = std::numeric_limits<double>::min();

  for (mwSize i = 0; i < Mx; i++) {
    xv->push_back(x[i]);
    yv->push_back(y[i]);
    zv->push_back(z[i]);
    
    // keep track of the interpolation domain boundaries. We are going
    // to need them to decide on the relative scale when computing
    // mba.MBAalg(). This will happen before we can use
    // e.g. surf.umin() to obtain that information
    xmin = std::min(xmin, x[i]);
    xmax = std::max(xmax, x[i]);
    ymin = std::min(ymin, y[i]);
    ymax = std::max(ymax, y[i]);
  }
  // interpolation points
  // boost::shared_ptr<VectorType> xiv(new VectorType);
  // boost::shared_ptr<VectorType> yiv(new VectorType);

  // for (mwSize i = 0; i < Mxi; i++) {
  //   xiv->push_back(xi[i]);
  //   yiv->push_back(yi[i]);
  // }

  // create the Multilevel B-spline object
  MBA mba(xv, yv, zv);

  // compute the interpolat
  if ((xmax-xmin)/(ymax-ymin) > 1.0) {
    mba.MBAalg((xmax-xmin)/(ymax-ymin), 1.0, 7);
  } else {
    mba.MBAalg(1.0, (ymax-ymin)/(xmax-xmin), 8);
  }

  // get the surface object
  UCBspl::SplineSurface surf = mba.getSplineSurface();

  // compute the interpolated surface value for each grid point, being
  // careful to return a NaN if the grid point is outside the
  // interpolation domain, because otherwise the MBA library seg faults
  for (mwSize i = 0; i < Mxi; i++) {
    if (xi[i] < xmin || xi[i] > xmax 
	|| yi[i] < ymin || yi[i] > ymax) {
      zi[i] = mxGetNaN();
    } else {
      zi[i] = surf.f(xi[i], yi[i]);
    }
  }

}

#endif /* MBA_SURFACE_INTERPOLATION_CPP */
