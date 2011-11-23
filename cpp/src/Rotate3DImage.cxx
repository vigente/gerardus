/*
 * Rotate3DImage.cxx
 *
 * Program to rotate a 3D image in space
 * 
 * Example of usage:
 * 
 *  $ ./rotate3DImage image.mha 0.7830 0.6100 0.1180 0.4400 -0.6790 0.5860 -0.4370 0.4070 0.8010 0.0028 0.0094 -0.0061
 * ($ ./rotate3DImage image.mha  a11   a21   a31   a12    a22   a32    a13   a23   a33    tx   ty   tz )
 * 
 * This rotates the 3D image contained in image.mha using the rotation matrix in extended form A
 * 
 *     [ 0.7830  0.4400 -0.4370  0.0028]
 * A = [ 0.6100 -0.6790  0.4070  0.0094]
 *     [ 0.1180  0.5860  0.8010 -0.0061]
 *     [ 0.0000  0.0000  0.0000  1.0000]
 * 
 * Note: The order a11, a21,... has been selected to make it compatible with Matlab's way of linearizing
 * matrices.
 * A is the transformation from input space to output space coordinates.
 * 
 * The rotation in extended form is
 * 
 *   Y = A*X
 * 
 * where Y'=[y' 1]', X'=[x' 1]' are voxel coordinates in output and input space, respectively. 
 * The extended rotation matrix can be expressed as
 * 
 *   A = [a   (I-a)*m]
 *       [0       1  ]
 * 
 * where a is a (3,3) rotation matrix, and m is the centre of rotation, or also as
 * 
 *   A = [a       t  ]
 *       [0       1  ]
 * 
 * where t can be seen as a translation.
 * 
 * The results are saved to file image-rotated.mha by default, although it's possible to specify the output
 * file name with argument -o --outfile.
 * 
 * The lower and upper bounds of the ouput image are computed so that the whole original image is contained
 * within the output one. This will generally result in a larger output image.
 * 
 * But it is also possible to specify the upper and lower bounds of the output image to crop it, e.g. with
 * arguments
 * 
 * --cxf 5.00e-3 --cxt 1.76e-2 --cyf 7.29e-3 --cyt 2.03e-2 --czf 3.76e-3 --czt 2.61e-2
 * 
 * You can provide any number of cropping boundaries, and those will override the internally computed ones.
 * 
 */ 
 
 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2010-2011 University of Oxford
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
 
#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

// C++ functions
#include <iostream>
#include <limits>
#include <algorithm>

// Boost Filesystem library
#include "boost/filesystem/path.hpp"
#include "boost/filesystem/convenience.hpp"
namespace fs = boost::filesystem;

// Command line parser header file
#include <tclap/CmdLine.h>

// ITK files
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
//#include "itkPointSet.h"
#include "itkTransformMeshFilter.h"
#include "itkMesh.h"

// entry point for the program
int main(int argc, char** argv)
{
    
    /*******************************/
    /** Command line parser block **/
    /*******************************/
    
    static const unsigned int   MatlabPrecision = 15; // number of decimal figures after the point in Matlab
    
    // command line input argument types and variables
    fs::path                            imPath;
    bool                                verbose;
    fs::path                            outImPath;
    float                               rotpVal[12]; // rotation around centroid matrix
    float                               cxf, cxt, cyf, cyt, czf, czt; // cropping coordinates
    float                               bg; // background intensity
    std::string                         interpType; // interpolator type
    float                               autoCrop;
    
    TCLAP::ValueArg< float > cropZToArg( "", "czt", "Crop Z-coordinate upper bound (to)", false, 0.0, "float" );
    TCLAP::ValueArg< float > cropZFromArg( "", "czf", "Crop Z-coordinate lower bound (from)", false, 0.0, "float" );
    TCLAP::ValueArg< float > cropYToArg( "", "cyt", "Crop Y-coordinate upper bound (to)", false, 0.0, "float" );
    TCLAP::ValueArg< float > cropYFromArg( "", "cyf", "Crop Y-coordinate lower bound (from)", false, 0.0, "float" );
    TCLAP::ValueArg< float > cropXToArg( "", "cxt", "Crop X-coordinate upper bound (to)", false, 0.0, "float" );
    TCLAP::ValueArg< float > cropXFromArg( "", "cxf", "Crop X-coordinate lower bound (from)", false, 0.0, "float" );

    TCLAP::ValueArg< float > autoCropArg( "a", "autocrop", "Percent of padding space left around the segmentation mask", false, 0.0, "float" );
    
    try {
        
        // Define the command line object, program description message, separator, version
        TCLAP::CmdLine cmd( "rotate3DImage: rotate a 3D image in space", ' ', "0.0" );
    
        // input argument: cropping coordinates
        cmd.add( cropZToArg );
        cmd.add( cropZFromArg );
        cmd.add( cropYToArg );
        cmd.add( cropYFromArg );
        cmd.add( cropXToArg );
        cmd.add( cropXFromArg );

        // input argument: filename of output segmentation mask
        TCLAP::ValueArg< std::string > outImPathArg( "o", "outfile", "Output image filename", false, "", "file" );
        cmd.add( outImPathArg );

        // input argument: filename of output segmentation mask
        TCLAP::ValueArg< float > bgArg( "b", "bkg", "Background intensity", false, 0.0, "bkg" );
        cmd.add( bgArg );

        // input argument: auto cropping
        cmd.add( autoCropArg );
    
        // input argument: interpolating type
        TCLAP::ValueArg< std::string > interpTypeArg( "i", "interp", "Interpolator type: bspline (default), nn", false, "bspline", "string" );
        cmd.add( interpTypeArg );

        // input argument: verbosity
        TCLAP::SwitchArg verboseSwitch( "v", "verbose", "Increase verbosity of program output", false );
        cmd.add( verboseSwitch );
    
        // input argument: rotation matrix
        TCLAP::UnlabeledValueArg< float > a11Arg( "a11", "(1, 1) element of rotation matrix", true, 0.0, "A11" );
        cmd.add( a11Arg );
        TCLAP::UnlabeledValueArg< float > a21Arg( "a21", "(2, 1) element of rotation matrix", true, 0.0, "A21" );
        cmd.add( a21Arg );
        TCLAP::UnlabeledValueArg< float > a31Arg( "a31", "(3, 1) element of rotation matrix", true, 0.0, "A31" );
        cmd.add( a31Arg );

        TCLAP::UnlabeledValueArg< float > a12Arg( "a12", "(1, 2) element of rotation matrix", true, 0.0, "A12" );
        cmd.add( a12Arg );
        TCLAP::UnlabeledValueArg< float > a22Arg( "a22", "(2, 2) element of rotation matrix", true, 0.0, "A22" );
        cmd.add( a22Arg );
        TCLAP::UnlabeledValueArg< float > a32Arg( "a32", "(3, 2) element of rotation matrix", true, 0.0, "A32" );
        cmd.add( a32Arg );

        TCLAP::UnlabeledValueArg< float > a13Arg( "a13", "(1, 3) element of rotation matrix", true, 0.0, "A13" );
        cmd.add( a13Arg );
        TCLAP::UnlabeledValueArg< float > a23Arg( "a23", "(2, 3) element of rotation matrix", true, 0.0, "A23" );
        cmd.add( a23Arg );
        TCLAP::UnlabeledValueArg< float > a33Arg( "a33", "(3, 3) element of rotation matrix", true, 0.0, "A33" );
        cmd.add( a33Arg );

        TCLAP::UnlabeledValueArg< float > txArg( "a14", "(1, 4) element of rotation matrix", true, 0.0, "A14" );
        cmd.add( txArg );
        TCLAP::UnlabeledValueArg< float > tyArg( "a24", "(2, 4) element of rotation matrix", true, 0.0, "A24" );
        cmd.add( tyArg );
        TCLAP::UnlabeledValueArg< float > tzArg( "a34", "(3, 4) element of rotation matrix", true, 0.0, "A34" );
        cmd.add( tzArg );

        // input argument: filename of input file
        TCLAP::UnlabeledValueArg< std::string > imPathArg( "image", "3D image", true, "", "file" );
        cmd.add( imPathArg );

        // Parse the command line arguments
        cmd.parse( argc, argv );

        // Get the value parsed by each argument
        imPath = fs::path( imPathArg.getValue() );
        outImPath = fs::path( outImPathArg.getValue() );
        verbose = verboseSwitch.getValue();
        bg = bgArg.getValue();
        interpType = interpTypeArg.getValue();
        autoCrop = autoCropArg.getValue();
                
        // the matrix is passed to the parameters vectorin row-major order 
        // (where the column index varies the fastest),
        // while the input arguments are in colum-major order, to make it
        // compatible with Matlab's way of linearizing matrices
        rotpVal[0]  = a11Arg.getValue();
        rotpVal[1]  = a12Arg.getValue();
        rotpVal[2]  = a13Arg.getValue();
        rotpVal[3]  = a21Arg.getValue();
        rotpVal[4]  = a22Arg.getValue();
        rotpVal[5]  = a23Arg.getValue();
        rotpVal[6]  = a31Arg.getValue();
        rotpVal[7]  = a32Arg.getValue();
        rotpVal[8]  = a33Arg.getValue();
        rotpVal[9]  = txArg.getValue();
        rotpVal[10] = tyArg.getValue();
        rotpVal[11] = tzArg.getValue();
        
        cxf = cropXFromArg.getValue();
        cxt = cropXToArg.getValue();
        cyf = cropYFromArg.getValue();
        cyt = cropYToArg.getValue();
        czf = cropZFromArg.getValue();
        czt = cropZToArg.getValue();
        
    } catch (const TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "Error parsing command line: " << std::endl 
        << e.error() << " for arg " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    
    /*******************************/
    /** Load input image block    **/
    /*******************************/

    static const unsigned int   Dimension = 3; // volume data dimension (i.e. 3D volumes)
    typedef double              TScalarType; // data type for scalars (e.g. point coordinates)
    
    typedef float                                        FloatPixelType;
    typedef itk::Image< FloatPixelType, 
                        Dimension >                      InputImageType;
    typedef InputImageType::SizeType                     InputSizeType;
    typedef itk::ImageFileReader< InputImageType >       ReaderType;

    // landmark I/O variables
    ReaderType::Pointer                             imReader;
        
    // image variables
    InputSizeType                           sizeIn;
    InputImageType::Pointer                 imIn;
    
    try {
        
        // create file readers
        imReader = ReaderType::New();
        
        // read input 3D images
        imReader->SetFileName( imPath.string() );
        if ( verbose ) {
            std::cout << "# Input image filename: " << imPath.string() << std::endl;
        }
        imReader->Update();
        
        // get input image
        imIn = imReader->GetOutput();
        
        // get image's size
        sizeIn = imIn->GetLargestPossibleRegion().GetSize();

        if ( verbose ) {
            std::cout << "# Input image dimensions: " << sizeIn[0] << "\t" 
                << sizeIn[1] << "\t" << sizeIn[2] << std::endl; 
        }

    } catch( const std::exception &e )  // catch any exceptions
    {
        std::cerr << "Error loading input image: " << std::endl 
        << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    /********************************************************************************/
    /** Rotate vertices of image frame to figure out how big it's going to be      **/
    /********************************************************************************/

    typedef unsigned short                               UShortPixelType;
    typedef itk::Image< UShortPixelType, 
                        Dimension >                      OutputImageType;
    typedef OutputImageType::SizeType                    OutputSizeType;
    typedef itk::AffineTransform< TScalarType, 
                                  Dimension >            TransformType;
    typedef itk::Index< Dimension >                      IndexType;
    typedef itk::Mesh< TScalarType, Dimension >          MeshType;
    typedef MeshType::PointType                          PointType;
    typedef itk::TransformMeshFilter< MeshType, 
                                      MeshType, 
                                      TransformType >    TransformMeshFilter;
    typedef itk::ImageRegionConstIterator< InputImageType > ConstIteratorType;
    
    PointType                         point;
    PointType                         minpoint, maxpoint;
    IndexType                         idx, minidx, maxidx;
    TransformType::Pointer            transform;
    MeshType::Pointer                 vertices;
//    MeshType::Pointer                 spacingOutMesh;
    TransformMeshFilter::Pointer      transformMesh;
//    TransformMeshFilter::Pointer      spacingTransformMesh;
    OutputSizeType                    sizeOut;
    InputImageType::PointType         originOut;
    TransformType::OutputVectorType   centroid;

    InputImageType::SpacingType spacing;  
    spacing = imIn->GetSpacing();
    const InputImageType::PointType&    origin  = imIn->GetOrigin();

    try {
        
        // init objects
        vertices = MeshType::New();
//        spacingOutMesh = MeshType::New();
        transform = TransformType::New();
        transformMesh = TransformMeshFilter::New();
//        spacingTransformMesh = TransformMeshFilter::New();
        
        // build mesh with the 8 vertices of the frame that contains the image
        // [0,0,0]
        vertices->SetPoint( 0, origin );
        // [0,Ny-1,0]
        idx[0] = 0; idx[1] = sizeIn[1]-1; idx[2] = 0;
        imIn->TransformIndexToPhysicalPoint( idx, point );
        vertices->SetPoint( 1, point );
        // [Nx-1,0,0]
        idx[0] = sizeIn[0]-1; idx[1] = 0; idx[2] = 0;
        imIn->TransformIndexToPhysicalPoint( idx, point );
        vertices->SetPoint( 2, point );
        // [Nx-1,Ny-1,0]
        idx[0] = sizeIn[0]-1; idx[1] = sizeIn[1]-1; idx[2] = 0;
        imIn->TransformIndexToPhysicalPoint( idx, point );
        vertices->SetPoint( 3, point );
        // [0,0,Nz-1]
        idx[0] = 0; idx[1] = 0; idx[2] = sizeIn[2]-1;
        imIn->TransformIndexToPhysicalPoint( idx, point );
        vertices->SetPoint( 4, point );
        // [0,Ny-1,Nz-1]
        idx[0] = 0; idx[1] = sizeIn[1]-1; idx[2] = sizeIn[2]-1;
        imIn->TransformIndexToPhysicalPoint( idx, point );
        vertices->SetPoint( 5, point );
        // [Nx-1,0,Nz-1]
        idx[0] = sizeIn[0]-1; idx[1] = 0; idx[2] = sizeIn[2]-1;
        imIn->TransformIndexToPhysicalPoint( idx, point );
        vertices->SetPoint( 6, point );
        // [Nx-1,Ny-1,Nz-1]
        idx[0] = sizeIn[0]-1; idx[1] = sizeIn[1]-1; idx[2] = sizeIn[2]-1;
        imIn->TransformIndexToPhysicalPoint( idx, point );
        vertices->SetPoint( 7, point );

        // initialize the affine transform to identity
        transform->SetIdentity();

        // initialize an array for the affine transform parameters
        transform->SetIdentity();
        TransformType::ParametersType rotp = transform->GetParameters();

        for (size_t i=0; i<12; ++i) {
            rotp[i] = rotpVal[i];
        }

        // resplace identity transform by the affine transform we want to apply
        transform->SetParameters(rotp);
        
        // apply affine transformation to the image frame's vertices
        transformMesh->SetTransform(transform);
        transformMesh->SetInput(vertices);
        
        transformMesh->Update();

        // find a Cartesian frame that encloses the rotated image frame
        for (int i=0; i<3; i++) // init points that delimitate the rotated frame
        {
            minpoint[i] = std::numeric_limits<TScalarType>::max();
            maxpoint[i] = std::numeric_limits<TScalarType>::min();
        }
        
        vertices->PrepareForNewData();
        vertices = transformMesh->GetOutput();

        for ( int i = 0; i <= 7; i++ )
        { 
            vertices->GetPoint( i, &point );
            
            minpoint[0] = std::min(minpoint[0], point[0]);
            minpoint[1] = std::min(minpoint[1], point[1]);
            minpoint[2] = std::min(minpoint[2], point[2]);

            maxpoint[0] = std::max(maxpoint[0], point[0]);
            maxpoint[1] = std::max(maxpoint[1], point[1]);
            maxpoint[2] = std::max(maxpoint[2], point[2]);
        }
        
// COMMENT OUT this block because it produces an incorrect output image        
//        // Find change in scaling
//        
//        // place a (1, 1, 1) on the rotation center
//        PointType centerRotation = transform->GetCenter();
//        PointType vec = centerRotation;
//        vec[0] += 1.0;
//        vec[1] += 1.0;
//        vec[2] += 1.0;
//
//        // make a mesh with the origin and the (1,1,1) point
//        spacingOutMesh->SetPoint(0, centerRotation);  
//        spacingOutMesh->SetPoint(1, vec);
//        
//        // apply the transformation to said mesh
//        spacingTransformMesh->SetTransform(transform);
//        spacingTransformMesh->SetInput(spacingOutMesh);
//        spacingTransformMesh->Update();
//        spacingOutMesh->PrepareForNewData();
//        spacingOutMesh = spacingTransformMesh->GetOutput();
//        
//        // remove the transformed center to see the change in the vector
//        spacingOutMesh->GetPoint(0, &centerRotation);  
//        spacingOutMesh->GetPoint(1, &vec);
//        vec[0] -= centerRotation[0];
//        vec[1] -= centerRotation[1];
//        vec[2] -= centerRotation[2];
//        
//        // change the scaling of the output image according to the transformation
//        spacing[0] /= vec[0];
//        spacing[1] /= vec[1];
//        spacing[2] /= vec[2];
        
        // Autocrop
        
        // if the user has entered an autocrop percentage, then we have 
        // to compute the dimensions of the segmentation mask
        if (autoCropArg.isSet()) {
            
            // swap max and min values, because we know that the segmentation mask has
            // to be within those values
            point = minpoint;
            minpoint = maxpoint;
            maxpoint = point;

            // we are going to reuse the vertices
            vertices->PrepareForNewData();
            
            // loop all voxels, and find a tight frame around the segmentation mask
            // (note, we are looping in the input image, but we have to convert the 
            // coordinates to output space)
            PointType point;
            ConstIteratorType iterator(imIn, imIn->GetLargestPossibleRegion());
            for (iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
            {
                
                // if the voxel has been marked as belonging to a landmark, 
                // see whether we have to change the cropping area
                if ( iterator.Get() )
                {
                    
                    // compute point coordinates in input space
                    imIn->TransformIndexToPhysicalPoint(iterator.GetIndex(), point);
                    
                    // convert coordinates to output space
                    point = transform->TransformPoint(point);
                    
                    // update tight frame, if necessary
                    minpoint[0] = std::min(minpoint[0], point[0]);
                    minpoint[1] = std::min(minpoint[1], point[1]);
                    minpoint[2] = std::min(minpoint[2], point[2]);
        
                    maxpoint[0] = std::max(maxpoint[0], point[0]);
                    maxpoint[1] = std::max(maxpoint[1], point[1]);
                    maxpoint[2] = std::max(maxpoint[2], point[2]);
                }
            } 
            
            // extend (or reduce) the thight frame according to the autocrop parameter
            PointType delta;
            for (size_t i = 0; i <  Dimension; ++i) {
                delta[i] = (maxpoint[i] - minpoint[i]) * autoCrop / 100.0;
                minpoint[i] -= delta[i];
                maxpoint[i] += delta[i];
            }
            
        } // if (autoCropArg.isSet())
        
        // if the user has entered cropping parameters, then they override
        // anything else computed so far
        if ( cropXFromArg.isSet() ) {  minpoint[0] = cxf;  }
        if ( cropYFromArg.isSet() ) {  minpoint[1] = cyf;  }
        if ( cropZFromArg.isSet() ) {  minpoint[2] = czf;  }
        if ( cropXToArg.isSet() ) {  maxpoint[0] = cxt;  }
        if ( cropYToArg.isSet() ) {  maxpoint[1] = cyt;  }
        if ( cropZToArg.isSet() ) {  maxpoint[2] = czt;  }
        
        // output cropping parameters used, in case they are needed for another image
        if (verbose) {
            std::cout.precision( MatlabPrecision );
            std::cout << "# --cxf " << minpoint[0]
                      << " --cxt "  << maxpoint[0]
                      << " --cyf " << minpoint[1]
                      << " --cyt "  << maxpoint[1]
                      << " --czf " << minpoint[2]
                      << " --czt "  << maxpoint[2] << std::endl;
        }
        
        // compute the size of the new frame (we use imIn because if we assume that the resolution
        // doesn't change, the out size should be correct even if the indices are obtained in
        // input space)
        imIn->TransformPhysicalPointToIndex( minpoint, minidx );
        imIn->TransformPhysicalPointToIndex( maxpoint, maxidx );
        sizeOut[0] = maxidx[0] - minidx[0] + 1;
        sizeOut[1] = maxidx[1] - minidx[1] + 1;
        sizeOut[2] = maxidx[2] - minidx[2] + 1;
        
        // move the new larger frame to its correct position
        originOut = minpoint;
        
    } catch( const std::exception &e )  // catch any exceptions
    {
        std::cerr << "Error finding rotated image limits: " << std::endl 
        << e.what() << std::endl;
        return EXIT_FAILURE;
    }


    /*******************************/
    /** Rotate image              **/
    /*******************************/

    typedef itk::ResampleImageFilter<
                  InputImageType, OutputImageType >      ResampleFilterType;
    // cubic spline
    typedef itk::BSplineInterpolateImageFunction< 
                  InputImageType, TScalarType >     BSplineInterpolatorType;
    typedef itk::NearestNeighborInterpolateImageFunction< 
                  InputImageType, TScalarType >     NearestNeighborInterpolatorType;
    typedef itk::InterpolateImageFunction< 
                  InputImageType, TScalarType >     InterpolatorType;

    // image variables
    OutputImageType::Pointer                             imOut;

    // filters
    ResampleFilterType::Pointer                          resampler;
    TransformType::Pointer                               transformInv;
    InterpolatorType::Pointer                            interpolator;

    try {

        // init objects
        transformInv = TransformType::New();
        
        // create objects for rotation
        resampler = ResampleFilterType::New();
        if (interpType == "bspline") {
            interpolator = BSplineInterpolatorType::New();
        } else if (interpType == "nn") {
            interpolator = NearestNeighborInterpolatorType::New();
        } else {
            throw std::string("Invalid interpolator type");
        }
        
        // the way ITK works, when you define a transform A and apply it to:
        //   * an image: it applies A^{-1} to the coordinates of voxels in output space to see
        //               which input coordinates they correspond to (and interpolate). This is 
        //               equivalent to applying A to the input coordinates
        //   * a mesh:   it applies A to the mesh coordinates, i.e. to the points in
        //               input space. Note that this is consistent with the image behaviour
        transform->GetInverse( transformInv );
        
        // set all the bits and pieces that go into the resampler
        resampler->SetDefaultPixelValue( bg );
        resampler->SetInterpolator( interpolator );
        resampler->SetTransform( transformInv );
        resampler->SetOutputOrigin( originOut );
        resampler->SetOutputSpacing( spacing );
        resampler->SetSize( sizeOut );
        resampler->SetInput( imIn ); 
        
        // rotate image
        resampler->Update();
        imOut = resampler->GetOutput();

        if ( verbose ) {
            std::cout << "# Output Image dimensions: " << sizeOut[0] << "\t" 
                << sizeOut[1] << "\t" << sizeOut[2] << std::endl; 
        }
        
    } catch( const std::exception &e )  // catch exceptions
    {
        std::cerr << "Error rotating input image: " << std::endl 
        << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch( const std::string &e )  // catch exceptions
    {
        std::cerr << "Error rotating input image: " << std::endl 
        << e << std::endl;
        return EXIT_FAILURE;
    }

    /*******************************/
    /** Output block              **/
    /*******************************/

    typedef itk::ImageFileWriter< OutputImageType >      WriterType;

    // I/O variables
    WriterType::Pointer                                  writer;
        
    try {     

        // create writer object        
        writer = WriterType::New();
        
        // create a filename for the output image by appending 
        // "rotated" to the input image filename, if none is
        // provided explicitely in the command line
        if ( outImPath.empty() ) {
            outImPath = imPath.branch_path() 
            / fs::path( fs::basename( imPath ) + "-rotated" 
            + fs::extension( imPath ) );
        }

        if ( verbose ) {
            std::cout << "# Output filename: " << outImPath.string() << std::endl;
        }
        
        // write output file
        writer->SetInput( imOut );
        writer->SetFileName( outImPath.string() );
        writer->SetUseCompression( true );
        writer->Update();
           
    } catch( const std::exception &e )  // catch any exceptions
    {
        std::cerr << "Error writing output image: " << std::endl 
        << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    /*******************************/
    /** End of program            **/
    /*******************************/
    
    return EXIT_SUCCESS; 
}
