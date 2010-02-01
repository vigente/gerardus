/*
 * Rotate3DImage.cxx
 *
 * Program to rotate a 3D image in space
 * 
 * Example of usage:
 * 
 *  $ ./rotate3DImage image.mha 0.011 0.013 0.014 0.783 0.610 0.118 0.440 -0.679 0.586 -0.437 0.407 0.801
 * ($ ./rotate3DImage image.mha  mx    my     mz   a11   a21   a31   a12    a22   a32    a13   a23   a33 )
 * 
 * This rotates the 3D image contained in image.mha around the centroid [0.011 0.013 0.014] using the rotation
 * matrix
 * 
 *     [ 0.783   0.440  -0.437 ]
 * A = [ 0.610  -0.679   0.407 ]
 *     [ 0.118   0.586   0.801 ]
 * 
 * That is, matrix values are provided in row-major (column index changes faster) order.
 * 
 * Note that this corresponds to first moving the image so that the centroid is placed at (0, 0, 0).
 * 
 * Then the rotation is applied to each voxel coordinate v as v' = A v. But note that because of the way ITK
 * works, the transformation is applied to the voxel coordinates of the *output* image. That is, 
 * the rotation will be the opposite of what you expect. This is not a big problem, because
 * the opposite rotation of A is A^T, where T means transpose. Thus, all you have to do is
 * provide the transpose of the rotation matrix if you have that issue.
 * 
 * Finally, the image is translated to the original centroid. 
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
  * Copyright © 2010 University of Oxford
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
#include "itkPointSet.h"
#include "itkTransformMeshFilter.h"
#include "itkMesh.h"

// entry point for the program
int main(int argc, char** argv)
{
    
    /*******************************/
    /** Command line parser block **/
    /*******************************/
    
    // command line input argument types and variables
    fs::path                            imPath;
    bool                                verbose;
    fs::path                            outImPath;
    float                               centroidVal[3]; // rotation centroid
    float                               rotpVal[9]; // rotation matrix
    float                               cxf, cxt, cyf, cyt, czf, czt;
    
    TCLAP::ValueArg< float > cropZToArg( "", "czt", "Crop Z-coordinate upper bound (to)", false, 0.0, "float" );
    TCLAP::ValueArg< float > cropZFromArg( "", "czf", "Crop Z-coordinate lower bound (from)", false, 0.0, "float" );
    TCLAP::ValueArg< float > cropYToArg( "", "cyt", "Crop Y-coordinate upper bound (to)", false, 0.0, "float" );
    TCLAP::ValueArg< float > cropYFromArg( "", "cyf", "Crop Y-coordinate lower bound (from)", false, 0.0, "float" );
    TCLAP::ValueArg< float > cropXToArg( "", "cxt", "Crop X-coordinate upper bound (to)", false, 0.0, "float" );
    TCLAP::ValueArg< float > cropXFromArg( "", "cxf", "Crop X-coordinate lower bound (from)", false, 0.0, "float" );
    
    try {
        
        // Define the command line object, program description message, separator, version
        TCLAP::CmdLine cmd( "rotate3DImage: rotate a 3D image in space", ' ', "0.0" );
    
        // input argument: filename of input segmentation mask
        TCLAP::UnlabeledValueArg< std::string > imPathArg( "image", "3D image", true, "", "file" );
        cmd.add( imPathArg );

        // input argument: rotation centroid coordinates
        TCLAP::UnlabeledValueArg< float > cxArg( "mx", "X-coordinate for rotation centroid", true, 0.0, "mx" );
        cmd.add( cxArg );
        TCLAP::UnlabeledValueArg< float > cyArg( "my", "Y-coordinate for rotation centroid", true, 0.0, "my" );
        cmd.add( cyArg );
        TCLAP::UnlabeledValueArg< float > czArg( "mz", "Z-coordinate for rotation centroid", true, 0.0, "mz" );
        cmd.add( czArg );

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

        // input argument: verbosity
        TCLAP::SwitchArg verboseSwitch( "v", "verbose", "Increase verbosity of program output", false );
        cmd.add( verboseSwitch );
    
        // Parse the command line arguments
        cmd.parse( argc, argv );

        // Get the value parsed by each argument
        imPath = fs::path( imPathArg.getValue() );
        outImPath = fs::path( outImPathArg.getValue() );
        verbose = verboseSwitch.getValue();
        
        centroidVal[0] = cxArg.getValue();
        centroidVal[1] = cyArg.getValue();
        centroidVal[2] = czArg.getValue();
        
        rotpVal[0] = a11Arg.getValue();
        rotpVal[1] = a21Arg.getValue();
        rotpVal[2] = a31Arg.getValue();
        rotpVal[3] = a12Arg.getValue();
        rotpVal[4] = a22Arg.getValue();
        rotpVal[5] = a32Arg.getValue();
        rotpVal[6] = a13Arg.getValue();
        rotpVal[7] = a23Arg.getValue();
        rotpVal[8] = a33Arg.getValue();
        
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
        imReader->SetFileName( imPath.file_string() );
        if ( verbose ) {
            std::cout << "# Input image filename: " << imPath.file_string() << std::endl;
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

    /****************************************************************************/
    /** Rotate vertices of image frame to figure out how big it's going to be  **/
    /****************************************************************************/

    typedef itk::PointSet< TScalarType, 
                            Dimension >                  PointSetType;
    typedef PointSetType::PointType                      PointType;
    typedef unsigned short                               UShortPixelType;
    typedef itk::Image< UShortPixelType, 
                        Dimension >                      OutputImageType;
    typedef OutputImageType::SizeType                    OutputSizeType;
    typedef itk::AffineTransform< TScalarType, 
                                  Dimension >            TransformType;
    typedef itk::Index< Dimension >                      IndexType;
    typedef itk::Mesh< TScalarType, Dimension >          MeshType;
    typedef itk::TransformMeshFilter< MeshType, 
                                      MeshType, 
                                      TransformType >    TransformMeshFilter;
    
    PointType                         point;
    IndexType                         idx, minidx, maxidx;
    TransformType::Pointer            transform, transformInv;
    MeshType::Pointer                 vertices;
    TransformMeshFilter::Pointer      transformMesh;
    OutputSizeType                    sizeOut;
    InputImageType::PointType         originOut;
    TransformType::OutputVectorType   centroid;

    const InputImageType::SpacingType&  spacing = imIn->GetSpacing();
    const InputImageType::PointType&    origin  = imIn->GetOrigin();

    try {
        
        // init objects
        vertices = MeshType::New();
        transform = TransformType::New();
        transformInv = TransformType::New();
        transformMesh = TransformMeshFilter::New();
        
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
  
        // compute affine transformation that will rotate the image around the 
        // segmentation mask centroid to make the heart vertical
        transform->SetIdentity();
        
        // set as center of rotation the centroid of the tissue voxels
        centroid[0] = centroidVal[0];
        centroid[1] = centroidVal[1];
        centroid[2] = centroidVal[2];
        transform->Translate( -centroid );
        
        TransformType::Pointer rot = TransformType::New();
        rot->SetIdentity();
        TransformType::ParametersType rotp = rot->GetParameters();

        rotp[0] = rotpVal[0];
        rotp[1] = rotpVal[1];
        rotp[2] = rotpVal[2];
        rotp[3] = rotpVal[3];
        rotp[4] = rotpVal[4];
        rotp[5] = rotpVal[5];
        rotp[6] = rotpVal[6];
        rotp[7] = rotpVal[7];
        rotp[8] = rotpVal[8];

        rot->SetParameters( rotp );

        // post-compose rotation with existing translation
        transform->Compose( rot, false );
        
        // post-compose translation back to original position
        transform->Translate( centroid, false );
        
        // because the affine transformation defined for the image maps coordinates in 
        // the output image to coordinates in the input image, we need the inverse 
        // transformation to figure out where the input frame will go in the output image
        transform->GetInverse( transformInv );
        
        // apply affine transformation to the image frame's vertices
        transformMesh->SetTransform( transformInv );
        transformMesh->SetInput( vertices );
        
        transformMesh->Update();
        
        // find a Cartesian frame that encloses the rotated image frame
        PointType minpoint, maxpoint;
        for (int i=0; i<3; i++) // init points that delimitate the rotated frame
        {
            minpoint[i] = std::numeric_limits<TScalarType>::max();
            maxpoint[i] = std::numeric_limits<TScalarType>::min();
        }
        
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
        
        // if the user has entered cropping parameters, then the corresponding 
        // ones computed for the output frame need to be overriden
        if ( cropXFromArg.isSet() ) {  minpoint[0] = cxf;  }
        if ( cropYFromArg.isSet() ) {  minpoint[1] = cyf;  }
        if ( cropZFromArg.isSet() ) {  minpoint[2] = czf;  }
        if ( cropXToArg.isSet() ) {  maxpoint[0] = cxt;  }
        if ( cropYToArg.isSet() ) {  maxpoint[1] = cyt;  }
        if ( cropZToArg.isSet() ) {  maxpoint[2] = czt;  }
        
        // compute the size of the new frame
        point[0] = maxpoint[0] - minpoint[0];
        point[1] = maxpoint[1] - minpoint[1];
        point[2] = maxpoint[2] - minpoint[2];
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
                       InputImageType, TScalarType >     InterpolatorType;
  

    // image variables
    OutputImageType::Pointer                             imOut;

    // filters
    ResampleFilterType::Pointer                          resampler;
    InterpolatorType::Pointer                            interpolator;

    try {
        
        // create objects for rotation
        resampler = ResampleFilterType::New();
        interpolator = InterpolatorType::New();
        
        // set all the bits and pieces that go into the resampler
        resampler->SetDefaultPixelValue( 50.0 );
        resampler->SetInterpolator( interpolator );
        resampler->SetTransform( transform );
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
        
    } catch( const std::exception &e )  // catch any exceptions
    {
        std::cerr << "Error rotating input image: " << std::endl 
        << e.what() << std::endl;
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
            std::cout << "# Output filename: " << outImPath.file_string() << std::endl;
        }
        
        // write output file
        writer->SetInput( imOut );
        writer->SetFileName( outImPath.file_string() );
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
