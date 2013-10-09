/*
 * CgalSurfaceSubdivision.cpp
 *
 * CGAL_SURFSUBDIVISION  3D Surface Approximating Subdivision Methods.
 *
 * This function is a Matlab wrapper of the CGAL 3D Surface Subdivision
 * Methods by Le-Jeng Andy Shiue:
 *
 * http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Subdivision_method_3/Chapter_main.html
 *
 * [TRI2, X2] = cgal_surfsubdivision(TRI, X, METHOD, ITER)
 *
 *   TRI is a 3-column matrix. Each row represents the indices of the three
 *   vertices that form a triangle. TRI as a whole represents the closed
 *   surface.
 *
 *   X is a 3-column matrix. Each row represents the Cartesian coordinates
 *   of a vertex on the surface, indexed by TRI values.
 *
 *   METHOD is a string to select one of the following subdivision schemes
 *   (all of them are approximating, not interpolating):
 *
 *     'CatmullClark': Catmull-Clark.
 *     'Loop':         Loop
 *     'DooSabin':     Doo-Sabin
 *     'Sqrt3':        √3
 *
 *   ITER is an integer with the number of iterations to run. With ITER=0,
 *   the output mesh is the same as the input mesh. With larger values of
 *   ITER, the output mesh will be smoother, at the cost of having more
 *   vertices. Usually, ITER=8 will produce already a very smooth mesh.
 *
 *   TRI2, X2 contain the output triangular mesh. Note that most subdivision
 *   methods will produce some non-triangular facets, even when the starting
 *   mesh is a triangulation. This function runs internally a triangulation
 *   process to convert the generalized output polyhedron into a triangular
 *   mesh.
 */

/*
 * Author: Ramon Casero <rcasero@gmail.com>
 * Copyright © 2013 University of Oxford
 * Version: 0.1.1
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

#ifndef CGALSURFACESUBDIVISION
#define CGALSURFACESUBDIVISION

/* mex headers */
#include <mex.h>

/* C++ headers */
#include <iostream>

/* Gerardus headers */
#include "MatlabImportFilter.h"
#include "MatlabExportFilter.h"

/* CGAL headers */
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Subdivision_method_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/triangulate_polyhedron.h>
#include "PolyhedronBuilder.h"

typedef MatlabImportFilter::MatlabInputPointer               MatlabInputPointer;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef CGAL::Polyhedron_3<Kernel>                           Polyhedron;
typedef CGAL::Point_3<Kernel>                                Point;
typedef Polyhedron::Facet                                    Facet;
typedef Polyhedron::Facet_iterator                           Facet_iterator;
typedef Polyhedron::Vertex_handle                            Vertex_handle;
typedef Polyhedron::Vertex_iterator                          Vertex_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator         Halfedge_around_facet_circulator;

/*
 * mexFunction(): entry point for the mex function
 */
void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {

  // interface to deal with input arguments from Matlab
  enum InputIndexType {IN_TRI, IN_X, IN_METHOD, IN_ITER, InputIndexType_MAX};
  MatlabImportFilter::Pointer matlabImport = MatlabImportFilter::New();
  matlabImport->ConnectToMatlabFunctionInput(nrhs, prhs);

  // check that we have all input arguments
  matlabImport->CheckNumberOfArguments(4, InputIndexType_MAX);

  // register the inputs for this function at the import filter
  MatlabInputPointer inTRI =        matlabImport->RegisterInput(IN_TRI, "TRI");
  MatlabInputPointer inX =          matlabImport->RegisterInput(IN_X, "X");
  MatlabInputPointer inMETHOD =     matlabImport->RegisterInput(IN_METHOD, "METHOD");
  MatlabInputPointer inITER =       matlabImport->RegisterInput(IN_ITER, "ITER");

  // interface to deal with outputs to Matlab
  enum OutputIndexType {OUT_TRI, OUT_X, OutputIndexType_MAX};
  MatlabExportFilter::Pointer matlabExport = MatlabExportFilter::New();
  matlabExport->ConnectToMatlabFunctionOutput(nlhs, plhs);

  // check number of outputs the user is asking for
  matlabExport->CheckNumberOfArguments(0, OutputIndexType_MAX);

  // register the outputs for this function at the export filter
  typedef MatlabExportFilter::MatlabOutputPointer MatlabOutputPointer;
  MatlabOutputPointer outTRI = matlabExport->RegisterOutput(OUT_TRI, "TRI");
  MatlabOutputPointer outX = matlabExport->RegisterOutput(OUT_X, "X");

  // if any of the inputs is empty, the output is empty too
  if (mxIsEmpty(prhs[IN_TRI]) || mxIsEmpty(prhs[IN_X])) {
    matlabExport->CopyEmptyArrayToMatlab(outTRI);
    matlabExport->CopyEmptyArrayToMatlab(outX);
    return;
  }

  // polyhedron to contain the input mesh
  Polyhedron mesh;
  PolyhedronBuilder<Polyhedron> builder(matlabImport, inTRI, inX);
  mesh.delegate(builder);

  // get size of input matrix with the points
  mwSize nrowsTri = mxGetM(inTRI->pm);
  mwSize nrowsX = mxGetM(inX->pm);
  
  // // DEBUG:
  // std::cout << "Number of facets read = " << mesh.size_of_facets() << std::endl;
  // std::cout << "Number of vertices read = " << mesh.size_of_vertices() << std::endl;

  if (nrowsTri != mesh.size_of_facets()) {
    mexErrMsgTxt(("Input " + inTRI->name + ": Number of triangles read into mesh different from triangles provided at the input").c_str());
  }
  if (nrowsX != mesh.size_of_vertices()) {
    mexErrMsgTxt(("Input " + inX->name + ": Number of vertices read into mesh different from vertices provided at the input").c_str());
  }

  // read input parameters
  unsigned int iter = matlabImport->ReadScalarFromMatlab<unsigned int>(inITER, 1);
  std::string method = matlabImport->ReadStringFromMatlab(inMETHOD, "CatmullClark");

  // compute subdivision surface
  if (method == "CatmullClark") {
    CGAL::Subdivision_method_3::CatmullClark_subdivision(mesh, iter);
  } else if (method == "Loop") {
    CGAL::Subdivision_method_3::Loop_subdivision(mesh, iter);
  } else if (method == "DooSabin") {
    CGAL::Subdivision_method_3::DooSabin_subdivision(mesh, iter);
  } else if (method == "Sqrt3") {
    CGAL::Subdivision_method_3::Sqrt3_subdivision(mesh, iter);
  } else {
    mexErrMsgTxt(("Input " + inMETHOD->name + ": Invalid method").c_str());
  } 

  // the subdivision mesh may have non-triangular facets. Split all facets to triangles
  CGAL::triangulate_polyhedron<Polyhedron>(mesh);

  // allocate memory for Matlab outputs
  double *tri = matlabExport->AllocateMatrixInMatlab<double>(outTRI, mesh.size_of_facets(), 3);
  double *x = matlabExport->AllocateMatrixInMatlab<double>(outX, mesh.size_of_vertices(), 3);

   // extract the vertices and triangles of the solution
  // snippet adapted from CgalMeshSegmentation.cpp

  // vertices coordinates. Assign indices to the vertices by defining
  // a map between their handles and the index
  std::map<Vertex_handle, int> V;
  int inum = 0;
  for(Vertex_iterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit) {

    // save to Matlab output
    x[inum] = vit->point().x();
    x[inum + mesh.size_of_vertices()] = vit->point().y();
    x[inum + 2*mesh.size_of_vertices()] = vit->point().z();

    // save to internal list of vertices
    V[vit] = inum++;

  }  

  // triangles given as (i,j,k), where each index corresponds to a vertex in x
  mwIndex row = 0;
  for (Facet_iterator fit = mesh.facets_begin(); fit != mesh.facets_end(); ++fit, ++row) {

    if (fit->facet_degree() != 3) {
      std::cerr << "Facet has " << fit->facet_degree() << " edges" << std::endl;
      mexErrMsgTxt("Facet does not have 3 edges");
    }

    // go around the half-edges of the facet, to extract the vertices
    Halfedge_around_facet_circulator heit = fit->facet_begin();
    int idx = 0;
    do {
      
      // extract triangle indices and save to Matlab output
      // note that Matlab indices go like 1, 2, 3..., while C++ indices go like 0, 1, 2...
      tri[row + idx * mesh.size_of_facets()] = 1 + V[heit->vertex()];
      idx++;

    } while (++heit != fit->facet_begin());
    
  }
}

#endif /* CGALSURFACESUBDIVISION */
