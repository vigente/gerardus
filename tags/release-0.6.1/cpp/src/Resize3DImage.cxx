/*
 * Resize3DImage.cxx
 *
 * Program to resize a 3D image 
 * 
 * Example of usage:
 * 
 *  $ ./resize3DImage 512 640 1024 image.mha 
 * 
 * This resizes the 3D image contained in image.mha to an output size
 * of 512 x 640 x 1024 voxels.
 * 
 * The results are saved to file image-resized.mha by default,
 * although it's possible to specify the output file name with
 * argument -o --outfile.
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
#include <cmath>

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
#include "itkBSplineInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkIdentityTransform.h"
#include "itkRecursiveGaussianImageFilter.h"

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
    std::string                         interpType; // interpolator type
    size_t                              sX, sY, sZ; // output size
    float                               sigX, sigY, sigZ; // user-defined Gaussian std
    bool                                sigmaSeg3D; // whether to use a very similar blurring to Seg3D's
    bool                                sigmaInVoxels; // whether sigma units are in voxels or real world coordinates
    
    try {
        
        // Define the command line object, program description message, separator, version
        TCLAP::CmdLine cmd( "resize3DImage: resize a 3D image", ' ', "0.0" );

	// input argument: override automatically computed sigma values for Gaussian
        // filter by user input
        TCLAP::ValueArg< float > sigXArg("", "sigx", "Gaussian std X", false, -1.0, "float");
        TCLAP::ValueArg< float > sigYArg("", "sigy", "Gaussian std Y", false, -1.0, "float");
        TCLAP::ValueArg< float > sigZArg("", "sigz", "Gaussian std Z", false, -1.0, "float");
	cmd.add(sigXArg);
	cmd.add(sigYArg);
	cmd.add(sigZArg);

        // input argument: Seg3D's low-pass blurring
        TCLAP::SwitchArg sigmaSeg3DSwitch("", "sigmaSeg3D", "Use similar low-pass blurring as Seg3D's Resample tool", false);
        cmd.add(sigmaSeg3DSwitch);
        
        // input argument: sigma units in voxels
        TCLAP::SwitchArg sigmaInVoxelsSwitch("", "sigmaInVoxels", "Sigma values provided by user are in voxels instead of real world coordinates", false);
        cmd.add(sigmaInVoxelsSwitch);
        
        // input argument: filename of output image
        TCLAP::ValueArg< std::string > outImPathArg("o", "outfile", "Output image filename", false, "", "file");
        cmd.add(outImPathArg);

        // input argument: interpolating type
        TCLAP::ValueArg< std::string > interpTypeArg("i", "interp", "Interpolator type: bspline (default), nn", false, "bspline", "string");
        cmd.add(interpTypeArg);

        // input argument: verbosity
        TCLAP::SwitchArg verboseSwitch("v", "verbose", "Increase verbosity of program output", false);
        cmd.add(verboseSwitch);
    
        // input argument: output size
        TCLAP::UnlabeledValueArg< size_t > sXArg("sx", "Output size X", true, 0, "sx");
        TCLAP::UnlabeledValueArg< size_t > sYArg("sy", "Output size Y", true, 0, "sy");
        TCLAP::UnlabeledValueArg< size_t > sZArg("sz", "Output size Z", true, 0, "sz");
        cmd.add(sXArg);
        cmd.add(sYArg);
        cmd.add(sZArg);

        // input argument: filename of input file
        TCLAP::UnlabeledValueArg< std::string > imPathArg("image", "3D image", true, "", "file");
        cmd.add(imPathArg);
        
        // Parse the command line arguments
        cmd.parse(argc, argv);

        // Get the value parsed by each argument
        imPath = fs::path(imPathArg.getValue());
        outImPath = fs::path(outImPathArg.getValue());
        verbose = verboseSwitch.getValue();
        interpType = interpTypeArg.getValue();
        sX = sXArg.getValue();
        sY = sYArg.getValue();
        sZ = sZArg.getValue();
        sigX = sigXArg.getValue();
        sigY = sigYArg.getValue();
        sigZ = sigZArg.getValue();
        sigmaSeg3D = sigmaSeg3DSwitch.getValue();
	sigmaInVoxels = sigmaInVoxelsSwitch.getValue();
        
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
    ReaderType::Pointer                                  imReader;
        
    // image variables
    InputSizeType                                        sizeIn;
    InputImageType::Pointer                              imIn;
    
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

    /*******************************/
    /** Smooth image              **/
    /*******************************/

    // [from ITK's /usr/share/doc/insighttoolkit3-examples/examples/Filtering/SubsampleVolume.cxx.gz]

    typedef itk::RecursiveGaussianImageFilter< 
                                  InputImageType,
                                  InputImageType >       GaussianFilterType;

    typedef unsigned short                               UShortPixelType;
    typedef itk::Image< UShortPixelType, 
                        Dimension >                      OutputImageType;
    typedef OutputImageType::SizeType                    OutputSizeType;

    // image variables
    InputImageType::SpacingType                          spacingIn;  
    OutputImageType::Pointer                             imOut;
    OutputSizeType                                       sizeOut;
 

    GaussianFilterType::Pointer                          smootherX;
    GaussianFilterType::Pointer                          smootherY;
    GaussianFilterType::Pointer                          smootherZ;

    // standard deviation for smoother 
    double sigmaX;
    double sigmaY;
    double sigmaZ;

    try {

        // get parameter values
        spacingIn = imIn->GetSpacing();

        // output size (if command line value is 0, then use input image size)
        if (sX == 0) {sizeOut[0] = sizeIn[0];} else {sizeOut[0] = sX;} 
        if (sY == 0) {sizeOut[1] = sizeIn[1];} else {sizeOut[1] = sY;}
        if (sZ == 0) {sizeOut[2] = sizeIn[2];} else {sizeOut[2] = sZ;}
        
        // instantiate smoother
        smootherX = GaussianFilterType::New();
        smootherY = GaussianFilterType::New();
        smootherZ = GaussianFilterType::New();
        
        // set image to smooth
	// this is delayed until 

        // set standard deviation for smoother 
        sigmaX = spacingIn[0] * (double)sizeIn[0] / (double)sizeOut[0];
        sigmaY = spacingIn[1] * (double)sizeIn[1] / (double)sizeOut[1];
        sigmaZ = spacingIn[2] * (double)sizeIn[2] / (double)sizeOut[2];
        if (sigmaSeg3D) {
            sigmaX *= 0.61;
            sigmaY *= 0.61;
            sigmaZ *= 0.61;
        }

	// override automatically computed values of Gaussian standard
	// deviation by user parameters
	if (sigX >= 0.0) {
	  sigmaInVoxels ? sigmaX = spacingIn[0] * sigX : sigmaX = sigX;
	}
	if (sigY >= 0.0) {
	  sigmaInVoxels ? sigmaY = spacingIn[1] * sigY : sigmaY = sigY;
	}
	if (sigZ >= 0.0) {
	  sigmaInVoxels ? sigmaZ = spacingIn[2] * sigZ : sigmaZ = sigZ;
	}

        smootherX->SetSigma(sigmaX);
        smootherY->SetSigma(sigmaY);
        smootherZ->SetSigma(sigmaZ);
        
        // "we instruct each one of the smoothing filters to act along a particular
        // direction of the image, and set them to use normalization across scale space
        // in order to prevent for the reduction of intensity that accompanies the
        // diffusion process associated with the Gaussian smoothing." (ITK's example)
        smootherX->SetDirection(0);
        smootherY->SetDirection(1);
        smootherZ->SetDirection(2);

        smootherX->SetNormalizeAcrossScale(true);
        smootherY->SetNormalizeAcrossScale(true);
        smootherZ->SetNormalizeAcrossScale(true);
        
        
    } catch( const std::exception &e )  // catch exceptions
    {
        std::cerr << "Error smoothing input image: " << std::endl 
        << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    
    /*******************************/
    /** Resize image              **/
    /*******************************/

    typedef itk::IdentityTransform< TScalarType, 
                                  Dimension >            IdentityTransformType;
    typedef itk::ResampleImageFilter<
                  InputImageType, OutputImageType >      ResampleFilterType;
    // cubic spline
    typedef itk::BSplineInterpolateImageFunction< 
                  InputImageType, TScalarType >          BSplineInterpolatorType;
    typedef itk::NearestNeighborInterpolateImageFunction< 
                  InputImageType, TScalarType >          NearestNeighborInterpolatorType;
    typedef itk::InterpolateImageFunction< 
                  InputImageType, TScalarType >          InterpolatorType;
                  

    OutputImageType::SpacingType                         spacingOut;  

    // filters
    IdentityTransformType::Pointer                       transform;
    ResampleFilterType::Pointer                          resampler;
    InterpolatorType::Pointer                            interpolator;

    try {

        // create objects for downsample
        transform = IdentityTransformType::New();
        resampler = ResampleFilterType::New();
        if (interpType == "bspline") {
            interpolator = BSplineInterpolatorType::New();
        } else if (interpType == "nn") {
            interpolator = NearestNeighborInterpolatorType::New();
        } else {
            throw std::string("Invalid interpolator type");
        }
        
        // compute spacing factor in the output image
        for (size_t i = 0; i < Dimension; ++i) { 
            spacingOut[i] = spacingIn[i] * (double)sizeIn[i] / (double)sizeOut[i];
        }

        // set all the bits and pieces that go into the resampler
        resampler->SetInterpolator(interpolator);
        resampler->SetTransform(transform);
        resampler->SetOutputOrigin(imIn->GetOrigin());
        resampler->SetOutputSpacing(spacingOut);
        resampler->SetSize(sizeOut);

	// create a pipeline for the image depending on which Gaussian
	// filters we are going to use
	if (sigmaZ > 0.0) {
	  resampler->SetInput(smootherZ->GetOutput());
	  if (sigmaY > 0.0) {
	    smootherZ->SetInput(smootherY->GetOutput());
	    if (sigmaX > 0.0) {
	      smootherY->SetInput(smootherX->GetOutput());
	      smootherX->SetInput(imIn);
	    } else {
	      smootherY->SetInput(imIn);
	    }
	  } else if (sigmaX > 0.0) {
	    smootherZ->SetInput(smootherX->GetOutput());
	    smootherX->SetInput(imIn);
	  } else {
	    smootherZ->SetInput(imIn);
	  }
	} else if (sigmaY > 0.0) {
	  resampler->SetInput(smootherY->GetOutput());
	  if (sigmaX > 0.0) {
	    smootherY->SetInput(smootherX->GetOutput());
	    smootherX->SetInput(imIn);
	  } else {
	    smootherY->SetInput(imIn);
	  }
	} else if (sigmaX > 0.0) {
	  resampler->SetInput(smootherX->GetOutput());
	  smootherX->SetInput(imIn);
	} else {
	  resampler->SetInput(imIn);
	}
        
        // rotate image
        resampler->Update();
        imOut = resampler->GetOutput();

        if ( verbose ) {
            std::cout << "# Output Image dimensions: " << sizeOut[0] << "\t" 
                << sizeOut[1] << "\t" << sizeOut[2] << std::endl; 
        }
        
    } catch( const std::exception &e )  // catch exceptions
    {
        std::cerr << "Error resizing input image: " << std::endl 
        << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch( const std::string &e )  // catch exceptions
    {
        std::cerr << "Error resizing input image: " << std::endl 
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
            / fs::path(fs::basename(imPath) + "-resized" 
            + fs::extension(imPath));
        }

        if ( verbose ) {
            std::cout << "# Output filename: " << outImPath.string() << std::endl;
        }
        
        // write output file
        writer->SetInput(imOut);
        writer->SetFileName(outImPath.string());
        writer->SetUseCompression(true);
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
        
