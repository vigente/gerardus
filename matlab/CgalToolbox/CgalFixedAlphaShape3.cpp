/* 
 *
 * CGAL_FIXED_ALPHA_SHAPE3  Individual alpha-shapes of a 3D set of points
 *
 * TRI = cgal_fixed_alpha_shape3(X, ALPHA)
 *
 *   X is a 3-column matrix. X(i, :) contains the xyz-coordinates of a 3D
 *   point.
 *
 *   ALPHA is a vector of scalar alpha values, alpha=R^2, where R is the
 *   probe radius.
 *
 *   TRI is a cell array of the same length as ALPHA. Cell TRI{i} contains
 *   the alpha shape triangulation for ALPHA{i}. Each row contains the 3
 *   nodes that form one triangular facet in the mesh. The i-th mesh can be
 *   visualised running:
 *
 *     >> trisurf(tri{i}, x)
 *
 *
 * This function uses CGAL's implementation of fixed alpha shapes [1]. Fixed
 * alpha shapes are more efficient when only the shape for one or a few
 * alpha values is required. When many alpha values are required, it may be
 * faster to use cgal_alpha_shape3().
 *
 * However, note that Matlab function alphavol() implemented by Jonas
 * Lundgren and provided as a third-party function in Gerardus seems to be
 * faster than either of the CGAL MEX functions, at least for a single alpha
 * value.
 *
 * [1] http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Alpha_shapes_3/Chapter_main.html
 *
 * See also: alphavol, cgal_alpha_shape3, scimat_lconvhull_smoothing
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2013 University of Oxford
  * Version: 0.3.1
  * $Rev$
  * $Date$
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

#ifndef CGALFIXEDALPHASHAPE3
#define CGALFIXEDALPHASHAPE3

/* mex headers */
#include <mex.h>

/* C++ headers */
#include <iostream>
#include <ctime> // DEBUG

/* Gerardus headers */
#include "MatlabImportFilter.h"
#include "MatlabExportFilter.h"

/* CGAL headers */
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
// vertex
typedef CGAL::Triangulation_vertex_base_with_info_3<mwSize, K>  Vb;
typedef CGAL::Fixed_alpha_shape_vertex_base_3<K, Vb> AsVb;
// cell
typedef CGAL::Fixed_alpha_shape_cell_base_3<K>       Fb;
// triangulation structure: vertex and cell
typedef CGAL::Triangulation_data_structure_3<AsVb, Fb> Tds;

// with an example of 412,068 points, using CGAL::Fast_location makes
// it slightly slower
//typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Delaunay;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef CGAL::Fixed_alpha_shape_3<Delaunay>          Alpha_shape_3;

typedef K::Point_3                                   Point;
typedef std::pair<Point, mwIndex>                    PointWithIndex;

typedef Alpha_shape_3::Facet                         Facet;

/*
 * mexFunction(): entry point for the mex function
 */
void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {

  // // DEBUG:
  // std::cout << "LLInitializing and reading data from Matlab" << std::endl;
  // clock_t time0 = clock();

  // interface to deal with input arguments from Matlab
  enum InputIndexType {IN_X, IN_ALPHA, InputIndexType_MAX};
  MatlabImportFilter::Pointer matlabImport = MatlabImportFilter::New();
  matlabImport->ConnectToMatlabFunctionInput(nrhs, prhs);

  // check the number of input arguments
  matlabImport->CheckNumberOfArguments(1, InputIndexType_MAX);

  // register the inputs for this function at the import filter
  typedef MatlabImportFilter::MatlabInputPointer MatlabInputPointer;
  MatlabInputPointer inX = matlabImport->RegisterInput(IN_X, "X");
  MatlabInputPointer inALPHA = matlabImport->RegisterInput(IN_ALPHA, "ALPHA");

  // interface to deal with outputs to Matlab
  enum OutputIndexType {OUT_TRI, OutputIndexType_MAX};
  MatlabExportFilter::Pointer matlabExport = MatlabExportFilter::New();
  matlabExport->ConnectToMatlabFunctionOutput(nlhs, plhs);

  // check number of outputs the user is asking for
  matlabExport->CheckNumberOfArguments(0, OutputIndexType_MAX);

  // register the outputs for this function at the export filter
  typedef MatlabExportFilter::MatlabOutputPointer MatlabOutputPointer;
  MatlabOutputPointer outTRI = matlabExport->RegisterOutput(OUT_TRI, "TRI");

  // if the set of points is empty, the outputs are empty too
  if (mxIsEmpty(prhs[IN_X])) {
    matlabExport->CopyEmptyArrayToMatlab(outTRI);
    return;
  }

  // default coordinates are NaN values, so that the user can spot
  // whether there was any problem reading them
  Point xDef(mxGetNaN(), mxGetNaN(), mxGetNaN());

  // get size of input matrix with the points
  mwSize nrowsX = mxGetM(prhs[IN_X]);
  mwSize ncolsX = mxGetN(prhs[IN_X]);
  if (ncolsX != 3) {
    mexErrMsgTxt("X must have 3 columns");
  }

  // read points from function
  std::vector<PointWithIndex> x(nrowsX);
  for (mwIndex i = 0; i < nrowsX; ++i) {

    // exit if user pressed Ctrl+C
    ctrlcCheckPoint(__FILE__, __LINE__);

    // read i-th row of input matrix as a point
    x[i] = std::make_pair(
			  matlabImport->ReadRowVectorFromMatlab<void, Point>(inX, i, xDef),
			  i+1 // because this will be a row index in Matlab, 1, ..., Nrows
			  );

    // don't accept NaNs or Infs, otherwise they will give a segfault
    if(mxIsNaN(x[i].first[0])
       || mxIsNaN(x[i].first[1])
       || mxIsNaN(x[i].first[2])
       || mxIsInf(x[i].first[0])
       || mxIsInf(x[i].first[1])
       || mxIsInf(x[i].first[2])
       ) {
      mexErrMsgTxt("X contains NaN or Inf values" );
    }
    
  }

  // // DEBUG:
  // for (std::vector<PointWithIndex>::iterator it = x.begin(); it != x.end(); ++it) {
  //   std::cout << "point " << it->second << " = " << it->first << std::endl;
  // }

  // // DEBUG
  // std::cout << "time = " << (double(clock() - time0) / CLOCKS_PER_SEC) << " sec" << std::endl;
  // time0 = clock();
  // std::cout << "Computing Delaunay triangulation" << std::endl;

  // Delaunay triangulation
  // http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Triangulation_3/Chapter_main.html#Subsection_39.5.3
  Delaunay delaunay(x.begin(), x.end());
  CGAL_assertion(delaunay.number_of_vertices() == nrowsX);

  // // DEBUG
  // std::cout << "Delaunay triangulation computed" << std::endl;
  // std::cout << "time = " << (double(clock() - time0) / CLOCKS_PER_SEC) << " sec" << std::endl;
  // time0 = clock();

  // read vector of alpha values provided by the user
  std::vector<double> alphaDef(1, 0.0);
  std::vector<double> alpha = matlabImport
    ->ReadArrayAsVectorFromMatlab<double, std::vector<double> >(inALPHA, alphaDef);

  // create output for surface triangulation that we are going to
  // extract
  const mwSize triDims[2] = {1, alpha.size()};
  plhs[OUT_TRI] = mxCreateCellArray(2, triDims);
  if (plhs[OUT_TRI] == NULL) {
    mexErrMsgTxt("Cannot allocate memory for output TRI" );
  }

  // for each alpha value provided by the user, compute the
  // corresponding alpha shape and extract its surface triangulation
  for (mwIndex i = 0; i < alpha.size(); ++i) {

    // exit if user pressed Ctrl+C
    ctrlcCheckPoint(__FILE__, __LINE__);
    
    // // DEBUG
    // std::cout << "Computing alpha shape in REGULARIZED mode by default"
    //           << std::endl;

    // the alpha shape destroys the Delaunay triangulation, so we need
    // to make a copy each time
    Delaunay delaunayCopy = delaunay;

    // compute alpha shape
    Alpha_shape_3 as(delaunayCopy, alpha[i]);
    
    // // DEBUG
    // std::cout << "Alpha shape computed in REGULARIZED mode by default"
    //           << std::endl;
    // std::cout << "time = " << (double(clock() - time0) / CLOCKS_PER_SEC) << " sec" << std::endl;
    // time0 = clock();
    
    // // DEBUG:
    // std::cout << "Extracting surface triangulation for alpha = " << alpha[i]
    // 	      << std::endl;

    // get alpha-shape surface
    std::list<Facet>       facets;
    as.get_alpha_shape_facets(std::back_inserter(facets),
			      Alpha_shape_3::REGULAR);

    // // DEBUG:
    // std::cout << "Number of facets = " << facets.size() << std::endl;

    // allocate memory in the current cell for the surface 
    double *triOut = matlabExport->AllocateMatrixInCellInMatlab<double>(outTRI, i, facets.size(), 3);

    // write facets to Matlab output
    mwSize row = 0; // row index of Matlab output
    for (std::list<Facet>::iterator it = facets.begin(); it != facets.end(); ++it) {

      // exit if user pressed Ctrl+C
      ctrlcCheckPoint(__FILE__, __LINE__);

      // vertex 1 of the triangle
      triOut[row] 
	= it->first->vertex(Delaunay::vertex_triple_index(it->second,0))->info();
      
      // vertex 2 of the triangle
      triOut[row+facets.size()] 
	= it->first->vertex(Delaunay::vertex_triple_index(it->second,1))->info();
      
      // vertex 3 of the triangle
      triOut[row+2*facets.size()] 
	= it->first->vertex(Delaunay::vertex_triple_index(it->second,2))->info();
      
      // // DEBUG:
      // std::cout << "facet = "
      // 		<< it->first->vertex(Delaunay::vertex_triple_index(it->second,0))->info()
      // 		<< ", "
      // 		<< it->first->vertex(Delaunay::vertex_triple_index(it->second,1))->info()
      // 		<< ", "
      // 		<< it->first->vertex(Delaunay::vertex_triple_index(it->second,2))->info()
      // 		<< std::endl; 
      
      // increase row counter
      ++row;

    } // end loop for each facet
    
    // // DEBUG
    // std::cout << "Surface triangulation extracted"
    //           << std::endl;
    // std::cout << "time = " << (double(clock() - time0) / CLOCKS_PER_SEC) << " sec" << std::endl;
    // time0 = clock();
    

  } // end loop for each alpha
  
}

#endif /* CGALFIXEDALPHASHAPE3 */
