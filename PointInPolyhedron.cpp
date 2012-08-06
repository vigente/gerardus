#ifndef POINTINPOLYHEDRON_CPP
#define POINTINPOLYHEDRON_CPP

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

/* mex headers */
#include <mex.h>

/* C++ headers */
#include <iostream>
// #include <cmath>
// #include <matrix.h>
// #include <vector>

/* ITK headers */
#include "MatlabImportFilter.h"

/* CGAL headers */
#include <CGAL/Simple_cartesian.h>

/*
 * mexFunction(): entry point for the mex function
 */
void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {

  typedef CGAL::Simple_cartesian<double> K;

  // interface to deal with input arguments from Matlab
  MatlabImportFilter::Pointer matlabImport = MatlabImportFilter::New();
  matlabImport->SetMatlabArgumentsPointer(nrhs, prhs);

  // check that we have a list of vertex coordinates, a list of
  // triangles that conform the surface, and a list of query points
  matlabImport->CheckNumberOfArguments(3, 3);

  typedef K::Point_3 Point;
  std::vector<Point> node;
  // matlabImport->
  // 		      GetVectorArgument<typename BoxFilterType::RadiusType, 
  // 					typename BoxFilterType::RadiusValueType>
  // 		      (2, "RADIUS", radius)

}

#endif /* POINTINPOLYHEDRON_CPP */
