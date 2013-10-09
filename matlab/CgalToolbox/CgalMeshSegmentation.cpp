/* CgalMeshSegmentation.cpp
 *
 * CGAL_MESHSEG  Surface meshing of an isosurface from a segmentation or
 * grayscale image
 *
 * This function is a Matlab wrapper of the CGAL 3D Surface Mesh Generation
 * for a grayscale input image.
 *
 * http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Surface_mesher/Chapter_main.html
 *
 * Note that even though CGAL internally associates x <-> rows, y <-> cols,
 * this functions presents to the user the usual Matlab convention of 
 * x <-> cols, y <-> rows.
 *
 * [TRI, X] = cgal_meshseg(IM, ISOVAL)
 *
 *   IM is a 3D array or a SCIMAT struct with a segmentation or a grayscale
 *   image.
 *
 *   ISOVAL is a scalar value that defines the isosurface to be meshed.
 *
 *   Note that if you have a binary segmentation (background=0, segmented
 *   voxels=1), ISOVAL=1 will give a very tight surface connecting the
 *   centres of the segmented voxels on the boundary. With small values,
 *   ISOVAL=0.01, the surface will be close to the centres of the adjacent
 *   background voxels. With ISOVAL=0.5, the surface will be halfway between
 *   the centres of the background and segmented voxels. The latter is
 *   usually the desired result.
 *
 *   TRI is a 3-column matrix. Each row represents the indices of the three
 *   vertices that form a triangle. TRI as a whole represents the closed
 *   surface.
 *
 *   X is a 3-column matrix. Each row represents the Cartesian coordinates
 *   of a vertex on the surface, indexed by TRI values.
 *
 * ... = cgal_meshseg(IM, ISOVAL, MINALPHA, MAXRAD, MAXD, C, MANIFOLD)
 *
 *   MINALPHA, MAXRAD, MAXD are scalars that implement the three meshing
 *   criteria in CGAL::Surface_mesh_default_criteria_3
 *
 *   http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Surface_mesher_ref/Class_Surface_mesh_default_criteria_3.html#Cross_link_anchor_1519
 *
 *   MINALPHA is "a lower bound on the minimum angle in degrees of the
 *   surface mesh facets". By default, MINALPHA = 30.0.
 *
 *   MAXRAD is "an upper bound on the radius of surface Delaunay balls. A
 *   surface Delaunay ball is a ball circumscribing a facet, centered on the
 *   surface and empty of vertices. Such a ball exists for each facet of the
 *   current surface mesh. Indeed the current surface mesh is the Delaunay
 *   triangulation of the current sampling restricted to the surface which
 *   is just the set of facets in the three dimensional Delaunay
 *   triangulation of the sampling that have a Delaunay surface ball". By
 *   default, MAXRAD is 1/2 of the minimum voxel size dimension. For
 *   example, if voxels have size [0.1 0.2 0.5], then by default
 *   MAXRAD=0.05.
 *
 *   MAXD is "an upper bound on the center-center distances of the surface
 *   mesh facets. The center-center distance of a surface mesh facet is the
 *   distance between the facet circumcenter and the center of its surface
 *   Delaunay ball". By default, MAXD is computed the same as MAXRAD.
 *
 *   C is a 3-vector with the coordinates of the centre of the bounding
 *   sphere used by the meshing algorithm. This is an important parameter.
 *   If C is close to the surface, it can produce lots of little triangles
 *   in that area. If C is outside the segmentation, the computed mesh may
 *   be incomplete.
 *
 *   MANIFOLD is a boolean flag.
 *
 *     MANIFOLD=false, the mesher uses CGAL::Non_manifold_tag: "When
 *     instantiated with the tag Non_manifold_tag the function template
 *     make_surface_mesh does not ensure that the output mesh is a manifold
 *     surface. The manifold property of output mesh may nevertheless result
 *     from the choice of appropriate meshing criteria".
 *
 *     MANIFOLD=true, the mesher uses CGAL::Manifold_tag: "When instantiated
 *     with the tag Manifold_tag the function template make_surface_mesh
 *     ensures that the output mesh is a manifold surface without boundary".
 *
 * Important!
 *
 * Note that this function can produce meshes with (1) stray vertices that
 * belong to no triangle, (2) triangles not oriented with respect to the
 * manifold, and (3) holes in the surface. These problems can be solved with
 * the following functions available in Gerardus:
 *
 *   [tri, x] = tri_squeeze(tri, x); % remove stray vertices
 *   [x, tri] = meshcheckrepair(x, tri, 'deep'); % correct orientation
 *   tri = cgal_tri_fillholes(tri, x); % fill holes in the surface
 *
 * See also: bwmesh, tri_squeeze, meshcheckrepair, cgal_tri_fillholes.
 */

/*
 * Author: Ramon Casero <rcasero@gmail.com>
 * Copyright © 2013 University of Oxford
 * Version: 0.1.4
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

// Warning! Within CGAL, x <-> rows, y <-> cols, whereas in Matlab, x
// <-> cols, y <-> rows. We maintain the interface with the user
// consistent with Matlab. That is, if the user provides a point's
// coordinates, he is providing (x,y,z) <-> (col,row,slice). It is up
// to us within this program to swap things around
// appropriately. Wherever we swap things so that the user can work
// consistently with Matlab, we add the notice:
// *swap to Matlab convention*

#ifndef CGALMESHSEGMENTATION
#define CGALMESHSEGMENTATION

/* mex headers */
#include <mex.h>

/* C++ headers */
#include <iostream>
#include <algorithm>
#include <limits>
#include <math.h> // DEBUG

/* Gerardus headers */
#include "MatlabImageHeader.h"
#include "MatlabImportFilter.h"
#include "MatlabExportFilter.h"

/* CGAL headers */
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Gray_level_image_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/ImageIO.h>

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

// labeled image
typedef Tr::Geom_traits GT;
typedef CGAL::Gray_level_image_3<GT::FT, GT::Point_3> Gray_level_image;
typedef CGAL::Implicit_surface_3<GT, Gray_level_image> Surface_3;

/*
 * mexFunction(): entry point for the mex function
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // interface to deal with input arguments from Matlab
  enum InputIndexType {IN_IM, IN_ISO, 
		       IN_MINALPHA, IN_MAXRAD, IN_MAXD, IN_C, IN_MANIFOLD, InputIndexType_MAX};
  MatlabImportFilter::Pointer matlabImport = MatlabImportFilter::New();
  matlabImport->ConnectToMatlabFunctionInput(nrhs, prhs);

  // check the number of input arguments
  matlabImport->CheckNumberOfArguments(2, InputIndexType_MAX);

  // register the inputs for this function at the import filter
  typedef MatlabImportFilter::MatlabInputPointer MatlabInputPointer;
  MatlabInputPointer inIM       = matlabImport->RegisterInput(IN_IM, "IM");
  MatlabInputPointer inISO      = matlabImport->RegisterInput(IN_ISO, "ISO");
  MatlabInputPointer inMINALPHA = matlabImport->RegisterInput(IN_MINALPHA, "MINALPHA");
  MatlabInputPointer inMAXRAD   = matlabImport->RegisterInput(IN_MAXRAD, "MAXRAD");
  MatlabInputPointer inMAXD     = matlabImport->RegisterInput(IN_MAXD, "MAXD");
  MatlabInputPointer inC        = matlabImport->RegisterInput(IN_C, "C");
  MatlabInputPointer inMANIFOLD = matlabImport->RegisterInput(IN_MANIFOLD, "C");

  // get input parameters
  double isoval   = matlabImport->ReadScalarFromMatlab<double>(inISO, 0.5);
  double minalpha = matlabImport->ReadScalarFromMatlab<double>(inMINALPHA, 30.0);
  bool asManifold = matlabImport->ReadScalarFromMatlab<bool>(inMANIFOLD, false);

  // interface to deal with outputs to Matlab
  enum OutputIndexType {OUT_TRI, OUT_X, OutputIndexType_MAX};
  MatlabExportFilter::Pointer matlabExport = MatlabExportFilter::New();
  matlabExport->ConnectToMatlabFunctionOutput(nlhs, plhs);

  // check number of outputs the user is asking for
  matlabExport->CheckNumberOfArguments(0, OutputIndexType_MAX);

  // register the outputs for this function at the export filter
  typedef MatlabExportFilter::MatlabOutputPointer MatlabOutputPointer;
  MatlabOutputPointer outTRI = matlabExport->RegisterOutput(OUT_TRI, "TRI");  
  MatlabOutputPointer outX   = matlabExport->RegisterOutput(OUT_X, "X");

  // if the image is empty, the output is empty
  if (mxIsEmpty(inIM->pm)) {
    matlabExport->CopyEmptyArrayToMatlab(outTRI);
    matlabExport->CopyEmptyArrayToMatlab(outX);
    return;
  }  

  // if the image is of type bool, the mesher enters an infinite loop
  // for some reason
  if ((mxGetClassID(inIM->pm) == mxLOGICAL_CLASS) || (mxGetClassID(inIM->pm) == mxINT64_CLASS)) {
      mexErrMsgTxt(("Input " + inIM->name + " has invalid type.").c_str());
  }

  // get Matlab image pointer and metadata in CGAL format
  _image *im = matlabImport->ReadCgalImageFromMatlab(inIM);

  // // DEBUG
  // //
  // // format     |   extension(s)   |  lecture  | ecriture 
  // // INRIMAGE   |  .inr[.gz]       |     X     |     X     -> + .gradient[.gz] + .gradient_direction[.gz]
  // // GIS        |  .dim, .ima[.gz] |     X     |     X
  // // ANALYZE    |  .hdr, .img[.gz] |     X     |     X
  // // PNM        |  .ppm, .pgm      |     X     |     X
  // // GIF        |  .gif            |     X     |     
  // // BMP        |  .gif            |     X     |     
  // _writeImage(im, "/tmp/foo.inr");

  // get two of the input arguments now that we know what the image
  // voxel size is
  double defRadAndD = std::min(std::min(im->vx, im->vy), im->vz);
  double maxrad   = matlabImport->ReadScalarFromMatlab<double>(inMAXRAD, defRadAndD * 0.5);
  double maxd     = matlabImport->ReadScalarFromMatlab<double>(inMAXD, defRadAndD * 0.5);

  // wrap as a Gray image so that we can pass it to the mesher
  Gray_level_image image(im, isoval);

  // compute the centre of mass of the segmentation. This is performed
  // computing the mean of the coordinates of all voxels in the
  // segmentation
  //
  // compute also the minimum and maximum voxels: these two voxels
  // form a tight rectangular box around the segmentation, and will be
  // used to define the radius of the boundary sphere for the meshing
  double xc = 0.0; // centroid coordinates
  double yc = 0.0;
  double zc = 0.0;
  mwSize nnz = 0; // number of segmented voxels
  double xmin = std::numeric_limits<double>::max(); // coordinates of minimum box boundary
  double ymin = std::numeric_limits<double>::max();
  double zmin = std::numeric_limits<double>::max();
  double xmax = std::numeric_limits<double>::min(); // coordinates of maximum box boundary
  double ymax = std::numeric_limits<double>::min();
  double zmax = std::numeric_limits<double>::min();

  for (mwIndex s = 0; s < image.zdim(); ++s) {
    for (mwIndex c = 0; c < image.ydim(); ++c) {
      for (mwIndex r = 0; r < image.xdim(); ++r) {
	
	if (image.value(r, c, s) != 0.0) {
	  
	  // contribution to the segmentation centroid (add the
	  // coordinates of the current voxel, as if the image had
	  // offset 0 and voxel size 1)
	  xc += r;
	  yc += c;
	  zc += s;

	  // count another segmented voxel
	  nnz++;

	  // contribution to the enclosing box boundaries
	  xmin = std::min(xmin, (double)r);
	  ymin = std::min(ymin, (double)c);
	  zmin = std::min(zmin, (double)s);
	  xmax = std::max(xmax, (double)r);
	  ymax = std::max(ymax, (double)c);
	  zmax = std::max(zmax, (double)s);
	}
      } // end for r
    } // end for c
  } // end for s

  // scale coordinates by voxel size
  xc *= image.vx();
  yc *= image.vy();
  zc *= image.vz();
  xmin *= image.vx();
  ymin *= image.vy();
  zmin *= image.vz();
  xmax *= image.vx();
  ymax *= image.vy();
  zmax *= image.vz();

 // divide by number of voxels in the segmentation to find centroid
 // coordinates if the image had offset (0,0,0)
  if (nnz != 0) {
    xc /= (double)nnz;
    yc /= (double)nnz;
    zc /= (double)nnz;
  }
  
  // DO NOT add image offset to the coordinates of the centre. CGAL
  // seems to do the computations ignoring the image offset and
  // referring everything to an offset of (0,0,0). The code that adds
  // the image offset is commented out in the next lines instead of
  // deleted, so that we have it for future reference.
  //
  xc += im->tx;
  yc += im->ty;
  zc += im->tz;

  // put them together
  GT::Point_3 defCentroid(yc, xc, zc); // *swap to Matlab convention*
  GT::Point_3 centroid = matlabImport->ReadRowVectorFromMatlab<void, GT::Point_3>
    (inC, defCentroid); // *centroid read with Matlab x/col, y/row convention*

  // // DEBUG
  // // note that centroid follows Matlab convention
  // std::cout << "Segmentation centre of mass (x/col, y/row, z/slice) without image offset = (" 
  // 	    << centroid.x() << ", " << centroid.y() << ", " << centroid.z() << ")" << std::endl;
  // // note that ymin, xmin follow CGAL convention
  // std::cout << "Segmentation box boundary (x/col, y/row, z/slice) without image offset = (" 
  // 	    << ymin << ", " << xmin << ", " << zmin << ")"
  // 	    << " to (" << ymax << ", " << xmax << ", " << zmax << ")"
  // 	    << std::endl;

  // the centre of the bounding sphere is the centroid (centre of
  // mass) of the segmentation
  // *swap centroid to CGAL convention*
  GT::Point_3 bounding_sphere_center(centroid.y(), centroid.x(), centroid.z());

  // the radius of the bounding sphere is the distance from the
  // centroid to the furthest boundary box corner. Then, we add 5%
  // more radius to avoid potential finite precision problems if the
  // radius is too tight.
  //
  // The padding is because the isosurface can be anywhere between the
  // segmented voxel's centre and the adjacent background voxel's
  // centre. With the padding, the bounding sphere will always contain
  // the isosurface, because it makes it go throught the adjacent
  // background voxel.
  double padding = 2 * std::max(std::max(im->vx, im->vy), im->vz);
  GT::FT bounding_sphere_squared_radius = 1.05 * std::max(
   (xc - xmin + padding)*(xc - xmin + padding) 
   + (yc - ymin + padding)*(yc - ymin + padding) 
   + (zc - zmin + padding)*(zc - zmin + padding),
   (xc - xmax + padding)*(xc - xmax + padding) 
   + (yc - ymax + padding)*(yc - ymax + padding) 
   + (zc - zmax + padding)*(zc - zmax + padding));

  // // DEBUG
  // std::cout << "Squared bounding sphere radius = " << bounding_sphere_squared_radius << std::endl;
  // std::cout << "Bounding sphere radius = " << sqrt(bounding_sphere_squared_radius) << std::endl;

  // instantiate the bounding sphere
  GT::Sphere_3 bounding_sphere(bounding_sphere_center,
			       bounding_sphere_squared_radius);

  // definition of the surface, with 10^-5 as relative precision
  Surface_3 surface(image, bounding_sphere, 1e-5);

  // variables to store the mesh as a triangulation
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3(tr);    // 2D-complex in 3D-Delaunay triangulation

  // Mesh criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(minalpha, maxrad, maxd);

  // Meshing
  if (asManifold) {
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
  } else {
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
  }

  // allocate memory for Matlab outputs
  double *tri = matlabExport->AllocateMatrixInMatlab<double>(outTRI, c2t3.number_of_facets(), 3);
  mwSize numOfVertices = std::distance(tr.finite_vertices_begin(), tr.finite_vertices_end());
  double *x = matlabExport->AllocateMatrixInMatlab<double>(outX, numOfVertices, 3);

  // extract the vertices and triangles of the solution
  // snippet copied from include/CGAL/IO/Complex_2_in_triangulation_3_file_writer.h

  typedef Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef Tr::Facet Facet;
  typedef Tr::Edge Edge;
  typedef Tr::Vertex_handle Vertex_handle;

  // vertices coordinates. Assign indices to the vertices by defining
  // a map between their handles and the index
  std::map<Vertex_handle, int> V;
  int inum = 0;
  for(Finite_vertices_iterator vit = tr.finite_vertices_begin();
      vit != tr.finite_vertices_end();
      ++vit) {

    // save to Matlab output
    // *swap to Matlab convention*
    //
    // Note that we need to swap because 
    //
    // vit->point().x() <--> rows
    // vit->point().y() <--> columns
    //
    // but in Matlab
    //
    // x(:, 1)          <--> columns
    // x(:, 2)          <--> rows
    //
    // Note also that the output of vit->point().[xyz]() ignores the
    // image offset, and is referred to an offset of (0,0,0)
    x[inum] = vit->point().y() + im->ty;
    x[inum + numOfVertices] = vit->point().x() + im->tx;
    x[inum + 2*numOfVertices] = vit->point().z() + im->tz;

    // save to internal list of vertices
    V[vit] = inum++;
  }

  // triangles given as (i,j,k), where each index corresponds to a vertex in x
  mwIndex row = 0;
  for (C2t3::Facet_iterator fit = c2t3.facets_begin(); fit != c2t3.facets_end(); ++fit, ++row) {

    // extract triangle indices and save to Matlab output
    // note that Matlab indices go like 1, 2, 3..., while C++ indices go like 0, 1, 2...
    const Tr::Cell_handle cell = fit->first;
    const int& index = fit->second;
    if (cell->is_facet_on_surface(index)==true) {
      tri[row]                             = 1 + V[cell->vertex(tr.vertex_triple_index(index, 0))];
      tri[row + c2t3.number_of_facets()]   = 1 + V[cell->vertex(tr.vertex_triple_index(index, 1))];
      tri[row + 2*c2t3.number_of_facets()] = 1 + V[cell->vertex(tr.vertex_triple_index(index, 2))];
    }

  }
  

}

#endif /* CGALMESHSEGMENTATION */
