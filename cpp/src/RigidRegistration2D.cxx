/*
 * RigidRegistration2D.cxx
 *
 * Program to align one 2D image to another using a rigid
 * transformation (rotation and translation)
 *
 * Example of usage:
 *
 * $ ./rigidRegistration2D source.bmp target.bmp -i 1000 -m 1.5
 *
 * This loads source.bmp and target.bmp assuming they are RGB images,
 * converts them to their luminance value (grayscale), computes a
 * rigid registration that stops either when the rotation step is
 * smaller than 1.5º, or after 1000 iterations, and creates a file
 * source-reg.bmp with the RGB solution.
 *
 * USAGE: 
 *
 *   cpp/src/rigidRegistration2D  [-v] [-o <file>] [-i <uint>] [-m <deg>] [-M
 *                                <deg>] [--] [--version] [-h] <source>
 *                                <target>
 *
 *
 * Where: 
 *
 *   -v,  --verbose
 *     Increase verbosity of program output
 *
 *   -o <file>,  --outfile <file>
 *     Output image filename
 *
 *   -i <uint>,  --maxiter <uint>
 *     Maximum number of iterations (default 200)
 *
 *   -m <deg>,  --minstep <deg>
 *     Minimum step length (default rotation 0.5º)
 *
 *   -M <deg>,  --maxstep <deg>
 *     Maximum step length (default rotation 10º)
 *
 *   --,  --ignore_rest
 *     Ignores the rest of the labeled arguments following this flag.
 *
 *   --version
 *     Displays version information and exits.
 *
 *   -h,  --help
 *     Displays usage information and exits.
 *
 *   <source>
 *     (required)  source 2D image
 *
 *   <target>
 *     (required)  target 2D image
 *
 *
 *   rigidRegistration2D:  rigid registration of two 2D images
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
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
#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkCenteredRigid2DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkVectorResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkMath.h"
#include "itkRGBPixel.h"
#include "itkRGBToLuminanceImageFilter.h"

// entry point for the program
int main(int argc, char** argv)
{
    
  /*******************************/
  /** Command line parser block **/
  /*******************************/
  
  // command line input argument types and variables
  fs::path                            imsPath, imtPath;
  bool                                verbose;
  fs::path                            outImPath;
  double                              minimumStepLength, maximumStepLength;
  unsigned int                        maximumNumberOfIterations;
  
  try {
    
    // Define the command line object, program description message, separator, version
    TCLAP::CmdLine cmd( "rigidRegistration2D:  rigid registration of two 2D images", ' ', "0.0" );
    
    // input argument: optimizer parameters
    TCLAP::ValueArg< double > maximumStepLengthArg("M", "maxstep", "Maximum step length (default rotation 10º)", false, 
						   10.0, "deg");
    TCLAP::ValueArg< double > minimumStepLengthArg("m", "minstep", "Minimum step length (default rotation 0.5º)", false, 
						   0.5, "deg");
    TCLAP::ValueArg< unsigned int > maximumNumberOfIterationsArg("i", "maxiter", "Maximum number of iterations (default 200)", false, 
						   200, "uint");
    cmd.add(maximumStepLengthArg);
    cmd.add(minimumStepLengthArg);
    cmd.add(maximumNumberOfIterationsArg);
    
    // input argument: filename of output image
    TCLAP::ValueArg< std::string > outImPathArg("o", "outfile", "Output image filename", false, "", "file");
    cmd.add(outImPathArg);
    
    // input argument: verbosity
    TCLAP::SwitchArg verboseSwitch("v", "verbose", "Increase verbosity of program output", false);
    cmd.add(verboseSwitch);
    
    // input argument: filename of input files, source and target
    TCLAP::UnlabeledValueArg< std::string > imsPathArg("source", "source 2D image", true, "", "source");
    cmd.add(imsPathArg);
    TCLAP::UnlabeledValueArg< std::string > imtPathArg("target", "target 2D image", true, "", "target");
    cmd.add(imtPathArg);
    
    // Parse the command line arguments
    cmd.parse(argc, argv);
    
    // Get the value parsed by each argument
    imsPath = fs::path(imsPathArg.getValue());
    imtPath = fs::path(imtPathArg.getValue());
    maximumStepLength = maximumStepLengthArg.getValue();
    minimumStepLength = minimumStepLengthArg.getValue();
    maximumNumberOfIterations = maximumNumberOfIterationsArg.getValue();
    outImPath = fs::path(outImPathArg.getValue());
    verbose = verboseSwitch.getValue();
  
  } catch (const TCLAP::ArgException &e) { // catch any exceptions
    
    std::cerr << "Error parsing command line: " << std::endl 
	      << e.error() << " for arg " << e.argId() << std::endl;
    return EXIT_FAILURE;
  }

  /*******************************/
  /** Load input images         **/
  /*******************************/
  
  static const unsigned int   Dimension = 2; // data dimension (i.e. 2D images)
  typedef double              TScalarType; // data type for scalars (e.g. point coordinates)
  typedef itk::RGBPixel<unsigned char> RGBPixelType; // pixel type (intensity values)
  typedef unsigned char GrayPixelType; // pixel type (intensity values)
  
  typedef itk::Image<RGBPixelType, Dimension>        InputImageType;
  typedef itk::Image<GrayPixelType, Dimension>       RegistrationImageType;
  typedef InputImageType::SizeType                   InputSizeType;
  typedef itk::ImageFileReader<InputImageType>       ReaderType;

  // landmark I/O variables
  ReaderType::Pointer                                sourceImageReader, targetImageReader;
  
  // image variables
  InputSizeType                                      sourceSize, targetSize;
  InputImageType::Pointer                            sourceImage, targetImage;
  
  try {
    
    // create file readers
    sourceImageReader = ReaderType::New();
    targetImageReader = ReaderType::New();
    
    // read input images
    sourceImageReader->SetFileName(imsPath.string());
    targetImageReader->SetFileName(imtPath.string());
    if (verbose) {
      std::cout << "# Source image filename: " << imsPath.string() << std::endl;
      std::cout << "# Target image filename: " << imtPath.string() << std::endl;
    }
    sourceImageReader->Update();
    targetImageReader->Update();
    
    // get input image
    sourceImage = sourceImageReader->GetOutput();
    targetImage = targetImageReader->GetOutput();
    
    // get image's size
    sourceSize = sourceImage->GetLargestPossibleRegion().GetSize();
    targetSize = targetImage->GetLargestPossibleRegion().GetSize();
    
    if (verbose) {
      std::cout << "# Source image dimensions: " << sourceSize << std::endl; 
      std::cout << "# Target image dimensions: " << targetSize << std::endl; 
      std::cout << "# Source image spacing: " << sourceImage->GetSpacing() << std::endl; 
      std::cout << "# Target image spacing: " << targetImage->GetSpacing() << std::endl; 
    }
    
  } catch( const std::exception &e ) { // catch any exceptions
    
    std::cerr << "Error loading input images: " << std::endl 
	      << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  
  /*******************************/
  /** Register images           **/
  /*******************************/

  typedef itk::RGBToLuminanceImageFilter<InputImageType,
					 RegistrationImageType > RGBToLuminanceFilterType;

  typedef itk::RegularStepGradientDescentOptimizer     OptimizerType;
  typedef itk::MeanSquaresImageToImageMetric< RegistrationImageType,
					      RegistrationImageType > MetricType;
  typedef itk:: LinearInterpolateImageFunction< RegistrationImageType, 
						TScalarType> InterpolatorType;
  typedef itk::ImageRegistrationMethod< RegistrationImageType,
					RegistrationImageType > RegistrationType;
  typedef itk::CenteredRigid2DTransform< TScalarType > TransformType;
  typedef itk::CenteredTransformInitializer<TransformType,
  					    RegistrationImageType,
  					    RegistrationImageType > TransformInitializerType;

  typedef OptimizerType::ScalesType OptimizerScalesType;

  // cast input image to a luminance image
  RGBToLuminanceFilterType::Pointer sourceCaster = RGBToLuminanceFilterType::New();
  RGBToLuminanceFilterType::Pointer targetCaster = RGBToLuminanceFilterType::New();
  sourceCaster->SetInput(sourceImage);
  targetCaster->SetInput(targetImage);
  sourceCaster->Update();
  targetCaster->Update();

  // instantiate registration components
  MetricType::Pointer metric = MetricType::New();
  TransformType::Pointer transform = TransformType::New();
  OptimizerType::Pointer optimizer = OptimizerType::New();
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  RegistrationType::Pointer registration = RegistrationType::New();
  
  // connect components to registration method
  registration->SetMetric(metric);
  registration->SetOptimizer(optimizer);
  registration->SetTransform(transform);
  registration->SetInterpolator(interpolator);

  // connect input images to registration method
  registration->SetFixedImage(targetCaster->GetOutput());
  registration->SetMovingImage(sourceCaster->GetOutput());

  // use whole image for registration
  registration->SetFixedImageRegion(targetCaster->GetOutput()->GetBufferedRegion());

  // initial parameters of the transformation
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();
  initializer->SetTransform(transform );
  initializer->SetFixedImage(targetCaster->GetOutput());
  initializer->SetMovingImage(sourceCaster->GetOutput());
  initializer->GeometryOn();
  initializer->InitializeTransform();
  // transform->SetAngle(0.0);

  if ( verbose ) {
    std::cout << "# Number of parameters: " 
	      << transform->GetNumberOfParameters() << std::endl;
    std::cout << "# Initial Rotation angle: " 
	      << transform->GetParameters()[0] / 3.14159265 * 180.0
	      << "º" << std::endl;
    std::cout << "# Initial Center of Rotation: " << transform->GetParameters()[1] 
	      << ", " << transform->GetParameters()[2] << std::endl;
    std::cout << "# Initial Translation: " << transform->GetParameters()[3] 
	      << ", " << transform->GetParameters()[4] << std::endl;
  }

  registration->SetInitialTransformParameters(transform->GetParameters());

  // optimizer parameters
  // From Luis Ibanez: "A typical rule of thumb is to assume that
  // rotating by 0.57 radians, (45 degrees) is as dramatic as
  // translating by half the image length.  Therefore, you want to
  // compute the length of your image in millimeters and compute its
  // ratio to the 0.57 radians".
  // http://itk-insight-users.2283740.n2.nabble.com/Confused-abour-Optimizer-Scales-td4010857.html
  // http://www.itk.org/pipermail/insight-users/2007-March/021435.html
  OptimizerScalesType optimizerScales(transform->GetNumberOfParameters());
  optimizerScales[0] = 1.0; // rotation
  optimizerScales[1] = (itk::Math::pi / 180.0 * 45.0) / (sourceSize[0] / 2.0 * sourceImage->GetSpacing()[0]); // center of rotation x
  optimizerScales[2] = (itk::Math::pi / 180.0 * 45.0) / (sourceSize[1] / 2.0 * sourceImage->GetSpacing()[1]); // center of rotation y
  optimizerScales[3] = (itk::Math::pi / 180.0 * 45.0) / (sourceSize[0] / 2.0 * sourceImage->GetSpacing()[0]); // translation x
  optimizerScales[4] = (itk::Math::pi / 180.0 * 45.0) / (sourceSize[1] / 2.0 * sourceImage->GetSpacing()[1]); // translation y
  optimizer->SetScales(optimizerScales);

  // for RegularStepGradientDescentOptimizer
  optimizer->SetMaximumStepLength(itk::Math::pi / 180.0 * maximumStepLength); // initial step change
  optimizer->SetMinimumStepLength(itk::Math::pi / 180.0 * minimumStepLength); // don't take steps smaller than this
  optimizer->SetNumberOfIterations(maximumNumberOfIterations);

  try {
    registration->Update();
    if (verbose) {
      std::cout << "# Final Rotation angle: " 
		<< transform->GetParameters()[0] / 3.14159265 * 180.0
		<< "º" << std::endl;
      std::cout << "# Final Center of Rotation: " << transform->GetParameters()[1] 
		<< ", " << transform->GetParameters()[2] << std::endl;
      std::cout << "# Final Translation: " << transform->GetParameters()[3] 
		<< ", " << transform->GetParameters()[4] << std::endl;
      std::cout << "# Stop condition: " 
      		<< optimizer->GetStopConditionDescription() << std::endl;
    }
  } catch (const std::exception &e) { // catch any exceptions
    
    std::cerr << "Error with registration: " << std::endl 
	      << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  /*******************************/
  /** Output block              **/
  /*******************************/
  
  typedef itk::VectorResampleImageFilter< InputImageType,
					  InputImageType >   ResampleFilterType;
  typedef itk::Image< RGBPixelType, Dimension > OutputImageType;
  typedef itk::ImageFileWriter< OutputImageType > WriterType;

  // I/O variables
  WriterType::Pointer                                  writer;
  
  try {     
    
    // create a filename for the output image by appending 
    // "reg" to the input image filename, if none is
    // provided explicitely in the command line
    if ( outImPath.empty() ) {
      outImPath = imsPath.branch_path() 
	/ fs::path(fs::basename(imsPath) + "-reg" 
		   + fs::extension(imsPath));
    }
    
    if ( verbose ) {
      std::cout << "# Output filename: " << outImPath.string() << std::endl;
    }
    
    // resampler filter
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetInput(sourceImage);
    resampler->SetTransform(registration->GetOutput()->Get());
    resampler->SetSize(targetImage->GetLargestPossibleRegion().GetSize());
    resampler->SetOutputOrigin(targetImage->GetOrigin());
    resampler->SetOutputSpacing(targetImage->GetSpacing());
    RGBPixelType background;
    background[0] = 0;
    background[1] = 0;
    background[2] = 0;
    resampler->SetDefaultPixelValue(background);

    // create writer object        
    writer = WriterType::New();
    
    // write output file
    writer->SetInput(resampler->GetOutput());
    writer->SetFileName(outImPath.string());
    writer->SetUseCompression(true);
    writer->Update();
    
  } catch( const std::exception &e ) { // catch any exceptions
    
    std::cerr << "Error writing output image: " << std::endl 
	      << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  
  /*******************************/
  /** End of program            **/
  /*******************************/
  
  return EXIT_SUCCESS; 
  
}
