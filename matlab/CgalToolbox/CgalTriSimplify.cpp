/*
 * CGAL_TRI_SIMPLIFY  Reduce number of faces in triangular mesh using edge
 * collapse.
 *
 * This function provides a Matlab interface to Fernando Cacciola's edge
 * collapse method to simplify a triangulated surface mesh available in
 * CGAL, using the default Lindstrom-Turk cost-strategy:
 *
 * http://doc.cgal.org/latest/Surface_mesh_simplification/index.html
 *
 * From the documentation: "Roughly speaking, the method consists of
 * iteratively replacing an edge with a single vertex, removing 2 triangles
 * per collapse."
 *
 * "Edges are collapsed according to a priority given by a user-supplied
 * cost function, and the coordinates of the replacing vertex are determined
 * by another user-supplied placement function. The algorithm terminates
 * when a user-supplied stop predicate is met, such as reaching the desired
 * number of edges."
 *
 * [TRI2, X2] = cgal_tri_simplify(TRI, X, RATIO)
 *
 *   TRI is a 3-column matrix. Each row contains the 3 nodes that form one
 *   triangular facet in the mesh.
 *
 *   X is a 3-column matrix. X(i, :) contains the xyz-coordinates of the
 *   i-th node in the mesh.
 *
 *   RATIO is a scalar with the stop criterion. The algorithm will stop when
 *   the number of current undirected edges is RATIO * number of original
 *   undirected edges.
 *
 *   TRI2, X2 is the description of the simplified output mesh.
 *
 * See also: cgal_surfsubdivision.
 */

/*
 * Author: Ramon Casero <rcasero@gmail.com>
 * Copyright © 2013 University of Oxford
 * Version: 0.1.0
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

//#define DEBUG

/* mex headers */
#include <mex.h>

/* C++ headers */
#include <iostream>

/* Gerardus headers */
#include "MatlabImportFilter.h"
#include "MatlabExportFilter.h"

/* CGAL headers */
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include "PolyhedronBuilder.h"

typedef MatlabImportFilter::MatlabInputPointer               MatlabInputPointer;

// Adaptor for Polyhedron_3
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Stop-condition policy
//#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef CGAL::Polyhedron_3<Kernel>                           Polyhedron; 
namespace SMS = CGAL::Surface_mesh_simplification;

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
  enum InputIndexType {IN_TRI, IN_X, IN_R, InputIndexType_MAX};
  MatlabImportFilter::Pointer matlabImport = MatlabImportFilter::New();
  matlabImport->ConnectToMatlabFunctionInput(nrhs, prhs);

  // check that we have all input arguments
  matlabImport->CheckNumberOfArguments(2, InputIndexType_MAX);

  // register the inputs for this function at the import filter
  MatlabInputPointer inTRI =        matlabImport->RegisterInput(IN_TRI, "TRI");
  MatlabInputPointer inX =          matlabImport->RegisterInput(IN_X, "X");
  MatlabInputPointer inR =          matlabImport->RegisterInput(IN_R, "R");

  // interface to deal with outputs to Matlab
  enum OutputIndexType {OUT_TRI, OUT_X, OutputIndexType_MAX};
  MatlabExportFilter::Pointer matlabExport = MatlabExportFilter::New();
  matlabExport->ConnectToMatlabFunctionOutput(nlhs, plhs);

  // check number of outputs the user is asking for
  matlabExport->CheckNumberOfArguments(0, OutputIndexType_MAX);

  // register the outputs for this function at the export filter
  typedef MatlabExportFilter::MatlabOutputPointer MatlabOutputPointer;
  MatlabOutputPointer outTRI = matlabExport->RegisterOutput(OUT_TRI, "TRI2");
  MatlabOutputPointer outX = matlabExport->RegisterOutput(OUT_X, "X2");

  // if any of the inputs is empty, the output is empty too
  if (mxIsEmpty(prhs[IN_TRI]) || mxIsEmpty(prhs[IN_X])) {
    matlabExport->CopyEmptyArrayToMatlab(outTRI);
    matlabExport->CopyEmptyArrayToMatlab(outX);
    return;
  }

  // read input parameters
  double ratio = matlabImport->ReadScalarFromMatlab<double>(inR, 0.1);

  // polyhedron to contain the input mesh
  Polyhedron mesh;
  PolyhedronBuilder<Polyhedron> builder(matlabImport, inTRI, inX);
  mesh.delegate(builder);

  // get size of input matrix with the points
  mwSize nrowsTri = mxGetM(inTRI->pm);
  mwSize nrowsX = mxGetM(inX->pm);

#ifdef DEBUG  
  std::cout << "Number of facets read = " << mesh.size_of_facets() << std::endl;
  std::cout << "Number of vertices read = " << mesh.size_of_vertices() << std::endl;
#endif

  if (nrowsTri != mesh.size_of_facets()) {
    mexErrMsgTxt(("Input " + inTRI->name + ": Number of triangles read into mesh different from triangles provided at the input").c_str());
  }
  if (nrowsX != mesh.size_of_vertices()) {
    mexErrMsgTxt(("Input " + inX->name + ": Number of vertices read into mesh different from vertices provided at the input").c_str());
  }

  // the simplification stops when the number of undirected edges
  // drops below some portion, e.g. 10%, of the initial count
  SMS::Count_ratio_stop_predicate<Polyhedron> stop(ratio);
  
  // This the actual call to the simplification algorithm.
  // The surface and stop conditions are mandatory arguments.
  // The index maps are needed because the vertices and edges
  // of this surface lack an "id()" field.
  SMS::edge_collapse(mesh, stop,
		     CGAL::vertex_index_map(boost::get(CGAL::vertex_external_index, mesh)) 
		     .edge_index_map(boost::get(CGAL::edge_external_index, mesh)));

  // allocate memory for Matlab outputs
  double *tri = matlabExport->AllocateMatrixInMatlab<double>(outTRI, mesh.size_of_facets(), 3);
  double *x = matlabExport->AllocateMatrixInMatlab<double>(outX, mesh.size_of_vertices(), 3);

   // extract the triangles of the solution
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
