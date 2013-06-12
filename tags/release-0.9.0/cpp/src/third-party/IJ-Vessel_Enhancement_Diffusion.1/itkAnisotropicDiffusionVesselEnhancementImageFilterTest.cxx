/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAnisotropicDiffusionVesselEnhancementImageFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2007/06/20 16:03:23 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


#include "itkAnisotropicDiffusionVesselEnhancementImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main(int argc, char* argv [] )
{
  if ( argc < 3 )
    {
    std::cerr << "Missing Parameters: " 
              << argv[0]
              << " Input_Image"
              << " Vessel_Enhanced_Output_Image [SigmaMin SigmaMax NumberOfScales NumberOfIteration]" << std::endl; 
    return EXIT_FAILURE;
    }
 
 
  // Define the dimension of the images
  const unsigned int Dimension = 3;
  typedef double      InputPixelType;
  typedef double      OutputPixelType;

  // Declare the types of the images
  typedef itk::Image< InputPixelType, Dimension>           InputImageType;
  typedef itk::Image< InputPixelType, Dimension>           OutputImageType;

  typedef itk::ImageFileReader< InputImageType  >      ImageReaderType;

  ImageReaderType::Pointer   reader = ImageReaderType::New();
  reader->SetFileName ( argv[1] ); 

  std::cout << "Reading input image : " << argv[1] << std::endl;
  try
    {
    reader->Update();
    }
  catch ( itk::ExceptionObject &err )
    {
    std::cerr << "Exception thrown: " << err << std::endl;
    return EXIT_FAILURE;
    }


  // Declare the anisotropic diffusion vesselness filter
  typedef itk::AnisotropicDiffusionVesselEnhancementImageFilter< InputImageType,
                                            OutputImageType>  VesselnessFilterType;

  // Create a vesselness Filter
  VesselnessFilterType::Pointer VesselnessFilter = 
                                      VesselnessFilterType::New();
  
//  VesselnessFilter->DebugOn();

  VesselnessFilter->SetInput( reader->GetOutput() );

  if ( argc >= 4 ) 
    { 
    VesselnessFilter->SetSigmaMin( atof(argv[3])  ); 
    }
 
  if ( argc >= 5 )
    {
    VesselnessFilter->SetSigmaMax( atof(argv[4]) ); 
    }

  if ( argc >= 6 )
    {
    VesselnessFilter->SetNumberOfSigmaSteps( atoi(argv[5]) ); 
    }

  if ( argc >= 7 )
    {
    VesselnessFilter->SetNumberOfIterations( atoi(argv[6]) ); 
    }

  //Test Set/Get VED parameters
 
  VesselnessFilter->SetSensitivity( 4.0 );
  VesselnessFilter->SetWStrength( 24.0 );
  VesselnessFilter->SetEpsilon( 0.01 );

  if ( fabs( VesselnessFilter->GetSensitivity() - 4.0 ) > 0.01 )
    {
    std::cerr << "Error Set/Get Sensitivity" << std::endl;
    return EXIT_FAILURE;
    }

  if ( fabs( VesselnessFilter->GetWStrength() - 24.0 ) > 0.01 )
    {
    std::cerr << "Error Set/Get Sensitivity" << std::endl;
    return EXIT_FAILURE;
    }

  if ( fabs( VesselnessFilter->GetEpsilon() - 0.01 ) > 0.01 )
    {
    std::cerr << "Error Set/Get Sensitivity" << std::endl;
    return EXIT_FAILURE;
    }

  VesselnessFilter->SetSensitivity( 5.0 );
  VesselnessFilter->SetWStrength( 25.0 );
  VesselnessFilter->SetEpsilon( 10e-2 );

  std::cout << "Enhancing vessels.........: " << argv[1] << std::endl;

  try
    {
    VesselnessFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Writing out the enhanced image to " <<  argv[2] << std::endl;

  typedef itk::ImageFileWriter< OutputImageType  >      ImageWriterType;
  ImageWriterType::Pointer writer = ImageWriterType::New();

  writer->SetFileName( argv[2] );
  writer->SetInput ( VesselnessFilter->GetOutput() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;

}

