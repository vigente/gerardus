/*
 * ITK_TRI_RASTERIZATION  Rasterization of triangular mesh to binary
 * segmentation
 *
 * This function provides a Matlab interface to
 * itk::TriangleMeshToBinaryImageFilter.
 *
 * BW = itk_tri_rasterization(TRI, X, RES, SIZE, ORIGIN)
 *
 *   TRI, X describe a triangular mesh.
 *
 *   TRI is a 3-column matrix. Each row represents the indices of the three
 *   vertices that form a triangle. TRI as a whole represents the closed
 *   surface.
 *
 *   X is a 3-column matrix. Each row represents the Cartesian coordinates
 *   of a vertex on the surface, indexed by TRI values.
 *
 *   RES is a 3-vector with the voxel size in (row, column, slice) format.
 *
 *   SIZE is a 3-vector with the number of voxels in (row, column, slice)
 *   format.
 *
 *   ORIGIN is a 3-vector with the coordinates of the centre of the
 *   bottom-left image voxel, in (x, y, z) format.
 *
 *   BW is the output uint8 binary segmentation. Voxels inside the mesh will
 *   be set to 1, and voxels outside to 0. Note that voxels with their
 *   centre exactly on the mesh offer different results depending on whether
 *   they are at the left, right, top or bottom of the mesh. This is a
 *   limitation in itk::TriangleMeshToBinaryImageFilter. See
 *   test_itk_tri_rasterization.m in Gerardus for an example.
 */

/*
 * Important programming note:
 *
 * Matlab saves images first by row and then by column. ITK expects
 * the opposite, but this generally doesn't matter, because what
 * happens is that ITK reads the image with rows and columns transpose
 * with respect to Matlab, processes it, and transposes back the
 * result.
 *
 * When apart from an image we also have to provide a parameter as a
 * (row, colum, slice) vector, there's no problem either.
 *
 * However, when the parameter we provide is in (x, y, z) format,
 * after reading it into a vector of vectors in C++, we need to swap
 * (x, y, z) to (y, x, z), so account for that transposition we
 * mentioned above. In Gerardus, we do that with
 * SwapXYInVectorOfVectors() and SwapXYInVector() from the
 * MatlabImportFilter.
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2013 University of Oxford
  * Version: 0.1.2
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

/* ITK headers */
#include <itkImage.h>
#include <itkMesh.h>
#include <itkTriangleMeshToBinaryImageFilter.h>
#include <itkTriangleCell.h>

// type definitions
static const unsigned int                       Dimension = 3;
typedef uint8_T                                 PixelType;
typedef float                                   MeshDataType;
typedef itk::PointSet<MeshDataType, Dimension>  PointSetType;
typedef PointSetType::CoordRepType              CoordType;
typedef itk::Mesh<MeshDataType, Dimension>      MeshType;
typedef itk::Image<PixelType, Dimension>        ImageType;
typedef PointSetType::PointType                 PointType;
typedef MeshType::CellType                      CellType;
typedef itk::TriangleCell<CellType>             TriangleType;
typedef CellType::CellAutoPointer               CellAutoPointer;
typedef itk::TriangleMeshToBinaryImageFilter<MeshType, ImageType> MeshFilterType;

/*
 * mexFunction(): entry point for the mex function
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // interface to deal with input arguments from Matlab
  enum InputIndexType {IN_TRI, IN_X, IN_RES, IN_SIZE, IN_ORIGIN, InputIndexType_MAX};
  MatlabImportFilter::Pointer matlabImport = MatlabImportFilter::New();
  matlabImport->ConnectToMatlabFunctionInput(nrhs, prhs);

  // check the number of input arguments
  matlabImport->CheckNumberOfArguments(2, InputIndexType_MAX);

  // register the inputs for this function at the import filter
  typedef MatlabImportFilter::MatlabInputPointer MatlabInputPointer;
  MatlabInputPointer inTRI = matlabImport->RegisterInput(IN_TRI, "TRI");
  MatlabInputPointer inX = matlabImport->RegisterInput(IN_X, "X"); // (x, y, z)
  MatlabInputPointer inRES = matlabImport->RegisterInput(IN_RES, "RES"); // (r, c, s)
  MatlabInputPointer inSIZE = matlabImport->RegisterInput(IN_SIZE, "SIZE"); // (r, c, s)
  MatlabInputPointer inORIGIN = matlabImport->RegisterInput(IN_ORIGIN, "ORIGIN"); // (x, y, z)

  // interface to deal with outputs to Matlab
  enum OutputIndexType {OUT_IM, OutputIndexType_MAX};
  MatlabExportFilter::Pointer matlabExport = MatlabExportFilter::New();
  matlabExport->ConnectToMatlabFunctionOutput(nlhs, plhs);
  
  // check that the number of outputs the user is asking for is valid
  matlabExport->CheckNumberOfArguments(0, OutputIndexType_MAX);

  // register the outputs for this function at the export filter
  typedef MatlabExportFilter::MatlabOutputPointer MatlabOutputPointer;
  MatlabOutputPointer outIM = matlabExport->RegisterOutput(OUT_IM, "IM");

  // if any input point set is empty, the outputs are empty too
  if (mxIsEmpty(inTRI->pm) || mxIsEmpty(inX->pm)) {
    matlabExport->CopyEmptyArrayToMatlab(outIM);
    return;
  }

  // get number of rows in inputs X and TRI
  mwSize nrowsX = mxGetM(inX->pm);
  mwSize nrowsTRI = mxGetM(inTRI->pm);

  // instantiate mesh
  MeshType::Pointer mesh = MeshType::New();

  // read vertices
  PointSetType::Pointer xDef = PointSetType::New(); // default: empty point set
  PointSetType::Pointer x = PointSetType::New();
  x->GetPoints()->CastToSTLContainer()
    = matlabImport->ReadVectorOfVectorsFromMatlab<PointType::CoordRepType, PointType>
    (inX, xDef->GetPoints()->CastToSTLContainer());

#ifdef DEBUG
  std::cout << "Number of X points read = " << x->GetNumberOfPoints() << std::endl;
#endif

  // assertion check
  if (nrowsX != x->GetNumberOfPoints()) {
    mexErrMsgTxt(("Input " + inX->name 
		  + ": Number of points read different from number of points provided by user").c_str()); 
  }

  // swap XY coordinates to make them compliant with ITK convention
  // (see important programming note at the help header above)
  matlabImport->SwapXYInVectorOfVectors<PointType::CoordRepType, std::vector<PointType> >
    (x->GetPoints()->CastToSTLContainer(), x->GetNumberOfPoints());

  // populate mesh with vertices
  mesh->SetPoints(x->GetPoints());

  // read triangles
  PointType triDef;
  triDef.Fill(mxGetNaN());
  for (mwIndex i = 0; i < nrowsTRI; ++i) {

    PointType triangle = matlabImport->ReadRowVectorFromMatlab<CoordType, PointType>(inTRI, i, triDef);

    // create a triangle cell to read the vertex indices of the current input triangle
    CellAutoPointer cell;
    cell.TakeOwnership(new TriangleType);

    // assign to the 0, 1, 2 elements in the triangle cell the vertex
    // indices that we have just read. Note that we have to substract
    // 1 to convert Matlab's index convention 1, 2, 3, ... to C++
    // convention 0, 1, 2, ...
    cell->SetPointId(0, triangle[0] - 1);
    cell->SetPointId(1, triangle[1] - 1);
    cell->SetPointId(2, triangle[2] - 1);

    // insert cell into the mesh
    mesh->SetCell(i, cell);
  }

#ifdef DEBUG
  std::cout << "Number of triangles read = " << mesh->GetNumberOfCells() << std::endl;
#endif

  // assertion check
  if (nrowsTRI != mesh->GetNumberOfCells()) {
    mexErrMsgTxt(("Input " + inTRI->name 
		  + ": Number of triangles read different from number of triangles provided by user").c_str()); 
  }

  // get user input parameters for the output rasterization
  ImageType::SpacingType spacingDef;
  spacingDef.Fill(1.0);
  ImageType::SpacingType spacing = matlabImport->
    ReadRowVectorFromMatlab<ImageType::SpacingValueType, ImageType::SpacingType>(inRES, spacingDef);

  ImageType::SizeType sizeDef;
  sizeDef.Fill(10);
  ImageType::SizeType size = matlabImport->
    ReadRowVectorFromMatlab<ImageType::SizeValueType, ImageType::SizeType>(inSIZE, sizeDef);

  ImageType::PointType originDef;
  originDef.Fill(0.0);
  ImageType::PointType origin = matlabImport->
    ReadRowVectorFromMatlab<ImageType::PointType::ValueType, ImageType::PointType>(inORIGIN, originDef);
  // (see important programming note at the help header above)
  matlabImport->SwapXYInVector<ImageType::PointType::ValueType, ImageType::PointType>(origin);

  // instantiate rasterization filter
  MeshFilterType::Pointer meshFilter = MeshFilterType::New();

  // smallest voxel side length
  ImageType::SpacingValueType minSpacing = spacing[0];
  for (mwIndex i = 1; i < Dimension; ++i) {
    minSpacing = std::min(minSpacing, spacing[i]);
  }

  // pass input parameters to the filter
  meshFilter->SetInput(mesh);
  meshFilter->SetSpacing(spacing);
  meshFilter->SetSize(size);
  meshFilter->SetOrigin(origin);
  meshFilter->SetTolerance(minSpacing / 10.0);
  meshFilter->SetInsideValue(1);
  meshFilter->SetOutsideValue(0);

  ImageType::IndexType start;
  start.Fill(0);
  meshFilter->SetIndex(start);

  // convert image size from itk::Size format to std::vector<mwSize>
  // so that we can use it in GraftItkImageOntoMatlab
  std::vector<mwSize> sizeStdVector(Dimension);
  for (unsigned int i = 0; i < Dimension; ++i) {
    sizeStdVector[i] = size[i];
  }

  // graft ITK filter output onto Matlab output
  matlabExport->GraftItkImageOntoMatlab<PixelType, Dimension>
    (outIM, meshFilter->GetOutput(), sizeStdVector);

#ifdef DEBUG
  std::cout << "Resolution (spacing) = " << meshFilter->GetSpacing() << std::endl;
  std::cout << "Size = " << meshFilter->GetSize() << std::endl;
  std::cout << "Origin = " << meshFilter->GetOrigin() << std::endl;
#endif
  
  // run rasterization
  meshFilter->Update();

}
