/*
 * Vesselness3DImage.cxx
 *
 * Program to compute the vesselness measure for each voxel of a 3D
 * image.
 *
 * Vesselness is a measure of how much does the voxel look like
 * belonging to a long cilindrical structure, as opposed to, e.g. a
 * blob or a plate.
 * 
 * Example of usage:
 * 
 *  $ ./vesselness3DImage -s .61e-4 image.mha 
 * 
 * This will highlight vessels with a diameter of approximately 
 * 4 * 0.61e-4 m = 0.24 mm.
 *
 * This program is basically a tuned wrapper of ITK's
 * HessianRecursiveGaussianImageFilter and
 * Hessian3DToVesselnessMeasureImageFilter.
 *
 * The vesselnes measure is corrected so that you will obtain the same
 * result independently from a linear scaling of intensities, or voxel
 * resolution. 
 *
 * Note that the measure is not necessarily between 0 and
 * 1, it could be smaller, so you need to take that into account.
 * 
 * The results are saved to file image-vesselness.mha by default,
 * although it's possible to specify the output file name with
 * argument -o --outfile.
 * 
 */ 
 
 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version 0.1
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
#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"

// entry point for the program
int main(int argc, char** argv) {
  
  /*******************************/
  /** Command line parser block **/
  /*******************************/
  
  // command line input argument types and variables
  fs::path                            imPath;
  bool                                verbose;
  fs::path                            outImPath;
  float                               sigma, alpha1, alpha2;
    
  try {
        
    // Define the command line object, program description message, separator, version
    TCLAP::CmdLine cmd( "vesselness3DImage: compute vesselness measure from a 3D image", ' ', "0.0" );

    // vesselness filters parameters
    TCLAP::ValueArg< float > sigmaArg("s", "sigma", "Gaussian std to compute Hessian (make it 1/4 of vessel diameter)", false, 1.0, "float");
    TCLAP::ValueArg< float > alpha1Arg("a", "alpha1", "Alpha_1 parameter (default 0.5)", false, 0.5, "float");
    TCLAP::ValueArg< float > alpha2Arg("b", "alpha2", "Alpha_2 parameter (default 2.0)", false, 2.0, "float");
    cmd.add(alpha2Arg);
    cmd.add(alpha1Arg);
    cmd.add(sigmaArg);

    // input argument: filename of output image
    TCLAP::ValueArg< std::string > outImPathArg("o", "outfile", "Output image filename", false, "", "file");
    cmd.add(outImPathArg);

    // input argument: verbosity
    TCLAP::SwitchArg verboseSwitch("v", "verbose", "Increase verbosity of program output", false);
    cmd.add(verboseSwitch);

    // input argument: filename of input file
    TCLAP::UnlabeledValueArg< std::string > imPathArg("image", "3D image", true, "", "file");
    cmd.add(imPathArg);
        
    // Parse the command line arguments
    cmd.parse(argc, argv);

    // Get the value parsed by each argument
    imPath = fs::path(imPathArg.getValue());
    outImPath = fs::path(outImPathArg.getValue());
    verbose = verboseSwitch.getValue();
    sigma = sigmaArg.getValue();
    alpha1 = alpha1Arg.getValue();
    alpha2 = alpha2Arg.getValue();
	
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
    
  typedef float                                        PixelType;
  typedef itk::Image<PixelType, Dimension>             ImageType;
  typedef ImageType::SizeType                          SizeType;
  typedef itk::ImageFileReader<ImageType>              ReaderType;

  // landmark I/O variables
  ReaderType::Pointer                                  imReader;
        
  // image variables
  SizeType                                             sizeIn;
  ImageType::Pointer                                   imIn;
    
  try {
        
    // create file readers
    imReader = ReaderType::New();
        
    // read input 3D images
    imReader->SetFileName(imPath.string());
    if ( verbose ) {
      std::cout << "# Input image filename: " 
		<< imPath.string() << std::endl;
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

  /********************************/
  /** Compute vesselness measure **/
  /********************************/

  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> 
    RescaleIntensityImageFilterType;
  typedef itk::HessianRecursiveGaussianImageFilter<ImageType>              
    HessianFilterType;
  typedef itk::Hessian3DToVesselnessMeasureImageFilter<PixelType> 
    VesselnessMeasureFilterType;
  typedef itk::MultiplyByConstantImageFilter<ImageType, TScalarType, ImageType> 
    MultiplyFilterType;
  
  MultiplyFilterType::Pointer
    multiplyFilter = MultiplyFilterType::New();

  try {

    HessianFilterType::Pointer 
      hessianFilter = HessianFilterType::New();
    VesselnessMeasureFilterType::Pointer 
      vesselnessFilter = VesselnessMeasureFilterType::New();
    RescaleIntensityImageFilterType::Pointer 
      intensityRescaleFilter = RescaleIntensityImageFilterType::New();
    // rescale intensities
    intensityRescaleFilter->SetInput(imIn);
    intensityRescaleFilter->SetOutputMaximum(1);
    intensityRescaleFilter->SetOutputMinimum(0);

    // Hessian filter
    hessianFilter->SetInput(intensityRescaleFilter->GetOutput());
    hessianFilter->SetSigma(sigma);

    // vesselness filter
    vesselnessFilter->SetInput(hessianFilter->GetOutput());
    vesselnessFilter->SetAlpha1(alpha1);
    vesselnessFilter->SetAlpha2(alpha2);

    // get spacing to normalize output of vesselness filter, because
    // this changes the scale of the vesselness values
    double dx = imIn->GetSpacing()[0];
    double dy = imIn->GetSpacing()[1];

    multiplyFilter->SetInput(vesselnessFilter->GetOutput());
    multiplyFilter->SetConstant(dx * dy);

    multiplyFilter->Update();

  } catch(const std::exception &e) { // catch any exceptions
    
    std::cerr << "Error computing vesselness measure: " << std::endl 
	      << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  
  /*******************************/
  /** Output block              **/
  /*******************************/
  
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
  // I/O variables
  WriterType::Pointer writer;
  
  try {     
    
    // create writer object        
    writer = WriterType::New();
    
    // create a filename for the output image by appending 
    // "skeleton" to the input image filename, if none is
    // provided explicitely in the command line
    if ( outImPath.empty() ) {
      outImPath = imPath.branch_path() 
	/ fs::path(fs::basename(imPath) + "-vesselness" 
		   + fs::extension(imPath));
    }
    
    if ( verbose ) {
      std::cout << "# Output filename: " 
		<< outImPath.string() << std::endl;
    }
    
    // write output file
    writer->SetInput(multiplyFilter->GetOutput());
    writer->SetFileName(outImPath.string());
    writer->SetUseCompression(true);
    writer->Update();
    
  } catch( const std::exception &e ) {  // catch any exceptions
    
    std::cerr << "Error writing output image: " << std::endl 
	      << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  
}
