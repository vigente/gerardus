/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2007/04/01 21:19:46 $
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


#include "itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkRescaleIntensityImageFilter.h"

int main(int argc, char* argv [] )
{
  if ( argc < 3 )
    {
    std::cerr << "Missing Parameters: " 
              << argv[0]
              << " Input_Image"
              << " Vessel_Enhanced_Output_Image [SigmaMin SigmaMax NumberOfScales]" << std::endl; 
    return EXIT_FAILURE;
    }
 
 
  // Define the dimension of the images
  const unsigned int Dimension = 3;
  typedef short      InputPixelType;
  typedef double     OutputVesselnessPixelType;

  // Declare the types of the images
  typedef itk::Image< InputPixelType, Dimension>            InputImageType;

  typedef itk::Image< OutputVesselnessPixelType, Dimension> VesselnessOutputImageType;

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


  // Declare the type of multiscale vesselness filter
  typedef itk::MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter<
                                            InputImageType,
                                            VesselnessOutputImageType>  
                                            MultiScaleVesselnessFilterType;

  // Create a vesselness Filter
  MultiScaleVesselnessFilterType::Pointer MultiScaleVesselnessFilter = 
                                      MultiScaleVesselnessFilterType::New();

  MultiScaleVesselnessFilter->SetInput( reader->GetOutput() );

  if ( argc >= 4 ) 
    { 
    MultiScaleVesselnessFilter->SetSigmaMin( atof(argv[3])  ); 
    }
 
  if ( argc >= 5 )
    {
    MultiScaleVesselnessFilter->SetSigmaMax( atof(argv[4]) ); 
    }

  if ( argc >= 6 )
    {
    MultiScaleVesselnessFilter->SetNumberOfSigmaSteps( atoi(argv[5]) ); 
    }

  try
    {
    MultiScaleVesselnessFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Writing out the enhanced image to " <<  argv[2] << std::endl;

  //Rescale the output of the vesslness image
  typedef itk::Image<unsigned char, 3>              OutputImageType; 
  typedef itk::RescaleIntensityImageFilter< VesselnessOutputImageType,
                                            OutputImageType> 
                                            RescaleFilterType;

  RescaleFilterType::Pointer rescale = RescaleFilterType::New();
  rescale->SetInput( MultiScaleVesselnessFilter->GetOutput() );
  rescale->SetOutputMinimum(   0 );
  rescale->SetOutputMaximum( 255 );
  rescale->Update();

  typedef itk::ImageFileWriter< OutputImageType  >      ImageWriterType;
  ImageWriterType::Pointer writer = ImageWriterType::New();

  writer->SetFileName( argv[2] );
  writer->SetInput ( rescale->GetOutput() );

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

