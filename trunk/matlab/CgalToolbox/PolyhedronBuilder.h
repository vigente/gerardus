/*
 * PolyhedronBuilder: class to read from Matlab two inputs TRI, X into
 * a surface CGAL::Polyhedron_3.
 *
 * This class has code specific to Gerardus and to CGAL.
 *
 * An example of how to use this class in a MEX Matlab function:
 *
 * #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
 * #include <CGAL/Polyhedron_3.h>
 * #include "PolyhedronBuilder.h"
 *
 * typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
 * typedef CGAL::Polyhedron_3<Kernel>                           Polyhedron;
 *
 * void mexFunction(int nlhs, mxArray *plhs[], 
 *		    int nrhs, const mxArray *prhs[]) {
 *
 *   // interface to deal with input arguments from Matlab
 *   enum InputIndexType {IN_TRI, IN_X, InputIndexType_MAX};
 *   MatlabImportFilter::Pointer matlabImport = MatlabImportFilter::New();
 *   matlabImport->ConnectToMatlabFunctionInput(nrhs, prhs);
 *
 *   // check that we have all input arguments
 *   matlabImport->CheckNumberOfArguments(2, InputIndexType_MAX);
 *
 *   // register the inputs for this function at the import filter
 *   MatlabInputPointer inTRI =    matlabImport->RegisterInput(IN_TRI, "TRI");
 *   MatlabInputPointer inX =      matlabImport->RegisterInput(IN_X, "X");
 * 
 *   // read input mesh into a polyhedron surface
 *   Polyhedron mesh;
 *   PolyhedronBuilder<Polyhedron> builder(matlabImport, inTRI, inX);
 *   mesh.delegate(builder);
 *
 * }
 *
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

#ifndef POLYHEDRONBUILDER_H
#define POLYHEDRONBUILDER_H

/* Gerardus headers */
#include "MatlabImportFilter.h"

/* C++ headers */
#include <exception>

/* CGAL headers */
#include <CGAL/Polyhedron_incremental_builder_3.h>

template <class Polyhedron>
class PolyhedronBuilder : public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS> {
public:
  typedef MatlabImportFilter::MatlabInputPointer               MatlabInputPointer;
  typedef typename Polyhedron::HalfedgeDS                      HalfedgeDS;

  MatlabImportFilter::Pointer matlabImport;
  MatlabInputPointer inTri;
  MatlabInputPointer inX;

  PolyhedronBuilder(MatlabImportFilter::Pointer _matlabImport, 
		    MatlabInputPointer _inTri, MatlabInputPointer _inX) 
    : matlabImport(_matlabImport), inTri(_inTri), inX(_inX) { }
  void operator()(HalfedgeDS& hds) {

    typedef typename HalfedgeDS::Vertex   Vertex;
    typedef typename Vertex::Point Point;

    // instantiate incremental builder
    CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> builder(hds, true);

    // get size of input matrix with the points
    mwSize nrowsTri = mxGetM(inTri->pm);
    mwSize ncolsTri = mxGetN(inTri->pm);
    mwSize nrowsX = mxGetM(inX->pm);
    mwSize ncolsX = mxGetN(inX->pm);
    if ((ncolsTri != 3) || (ncolsX != 3)) {
      mexErrMsgTxt("TRI and X inputs must have 3 columns");
    }

    // allocate space in the polyhedron
    builder.begin_surface(nrowsX, nrowsTri);

    // default coordinates are NaN values, so that the user can spot
    // whether there was any problem reading them
    Point xDef(mxGetNaN(), mxGetNaN(), mxGetNaN());

    // add mesh vertices
    for (mwIndex i = 0; i < nrowsX; ++i) {

      // exit if user pressed Ctrl+C
      ctrlcCheckPoint(__FILE__, __LINE__);
      
      // get coordinates of the vertex
      Point x = matlabImport->ReadRowVectorFromMatlab<void, Point>(inX, i, xDef);

      if (mxIsNaN(x.x()) || mxIsNaN(x.y()) || mxIsNaN(x.z())) {
	mexErrMsgTxt(("Input " + inX->name + ": Vertex coordinates are NaN").c_str());
      }
      
      // add vertex to the mesh
      builder.add_vertex(x);
    }
    
    // add mesh triangles
    for (mwIndex i = 0; i < nrowsTri; ++i) {
      
      // exit if user pressed Ctrl+C
      ctrlcCheckPoint(__FILE__, __LINE__);
      
      // get indices of the 3 vertices of each triangle. These indices
      // follow Matlab's convention v0 = 1, 2, ..., n
      mwIndex v0 = matlabImport->ReadScalarFromMatlab<mwIndex>(inTri, i, 0, mxGetNaN());
      mwIndex v1 = matlabImport->ReadScalarFromMatlab<mwIndex>(inTri, i, 1, mxGetNaN());
      mwIndex v2 = matlabImport->ReadScalarFromMatlab<mwIndex>(inTri, i, 2, mxGetNaN());
      if (mxIsNaN(v0) || mxIsNaN(v1) || mxIsNaN(v2)) {
	mexErrMsgTxt(("Input " + inTri->name + ": Triangle indices are NaN").c_str());
      }

      // add triangle. Note that we have to substract 1 from the index
      // to change from Matlab index convention (1, 2, 3, ...) to C++
      // index convention (0, 1, 2, ...)
      builder.begin_facet();
      builder.add_vertex_to_facet(v0-1);
      builder.add_vertex_to_facet(v1-1);
      builder.add_vertex_to_facet(v2-1);
      builder.end_facet();
      
      // if the facet couldn't be added successfully, we exit with an
      // error. Otherwise, when we try builder.end_surface(), Matlab
      // will throw a segmentation error
      if (builder.error()) {
	mexErrMsgTxt(("Inputs " + inTri->name + " and " + inX->name 
		      + " do not form a valid polyhedron").c_str());
      }
 
    }
     
    // finish up surface
    builder.end_surface();

  }
};

#endif /* POLYHEDRONBUILDER_H */
