/*
 * Skeletonize3DSegmentation.cxx
 *
 * Program to reduce a segmentation to its skeleton.
 * 
 * Example of usage:
 * 
 *  $ ./skeletonize3DSegmentation seg.mha 
 * 
 * This extracts the skeleton from the segmentation in seg.mha, using
 * Hanno Homann's implementation of a 3D thinning algorithm [1].
 *
 * [1] T.C. Lee, R.L. Kashyap, and C.N. Chu. Building skeleton models
 * via 3-D medial surface/axis thinning algorithms. Computer Vision,
 * Graphics, and Image Processing, 56(6):462–478, 1994.
 * 
 * The results are saved to file seg-skeleton.mha by default,
 * although it's possible to specify the output file name with
 * argument -o --outfile.
 * 
 */ 
 
 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
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
#include "itkBinaryThinningImageFilter3D.h"

// entry point for the program
int main(int argc, char** argv) {
  /*******************************/
  /** Command line parser block **/
  /*******************************/
  
  // command line input argument types and variables
  fs::path                            maskPath;
  bool                                verbose;
  fs::path                            outMaskPath;
  
  try {
    
    // Define the command line object, program description message,
    // separator, version
    TCLAP::CmdLine cmd( "skeletonize3DSegmentation: Reduce segmentation to its skeleton", ' ', "0.0" );
    
    // input argument: filename of output image
    TCLAP::ValueArg< std::string > 
      outMaskPathArg("o", "outfile", "Output image filename", false, "", "file");
    cmd.add(outMaskPathArg);
    
    // input argument: verbosity
    TCLAP::SwitchArg verboseSwitch("v", "verbose", "Increase verbosity of program output", false);
    cmd.add(verboseSwitch);
    
    // input argument: filename of input file
    TCLAP::UnlabeledValueArg< std::string > maskPathArg("image", "3D image", true, "", "file");
    cmd.add(maskPathArg);
    
    // Parse the command line arguments
    cmd.parse(argc, argv);
    
    // Get the value parsed by each argument
    maskPath = fs::path(maskPathArg.getValue());
    outMaskPath = fs::path(outMaskPathArg.getValue());
    verbose = verboseSwitch.getValue();
    
  } catch (const TCLAP::ArgException &e) { // catch any exceptions
    std::cerr << "Error parsing command line: " << std::endl 
	      << e.error() << " for arg " << e.argId() << std::endl;
    return EXIT_FAILURE;
  }
  

  /*******************************/
  /** Load input image block    **/
  /*******************************/
  
  static const unsigned int Dimension = 3; // volume data dimension (i.e. 3D volumes)
  typedef double TScalarType; // data type for scalars
  //  typedef bool BinaryPixelType;
  typedef unsigned short PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef ImageType::SizeType BinarySizeType;
  typedef itk::ImageFileReader< ImageType > BinaryReaderType;
  typedef itk::ImageRegionConstIterator< ImageType > ConstBinaryIteratorType;
  
  // landmark I/O variables
  BinaryReaderType::Pointer maskReader;
  
  try {
    
    // create file readers
    maskReader = BinaryReaderType::New();
    
    // read input 3D images
    maskReader->SetFileName( maskPath.file_string() );
    if ( verbose ) {
      std::cout << "# Segmentation mask filename: " 
		<< maskPath.file_string() << std::endl;
    }
    maskReader->Update();
    
    
  } catch( const std::exception &e ) { // catch any exceptions
    std::cerr << "Error loading input binary image: " << std::endl 
	      << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  
  /*******************************/
  /** Skeletonize image         **/
  /*******************************/

  typedef itk::BinaryThinningImageFilter3D< ImageType, 
    ImageType > ThinningFilterType;
  ThinningFilterType::Pointer 
    thinningFilter = ThinningFilterType::New();
  thinningFilter->SetInput(maskReader->GetOutput());
  thinningFilter->Update();

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
    if ( outMaskPath.empty() ) {
      outMaskPath = maskPath.branch_path() 
	/ fs::path(fs::basename(maskPath) + "-skeleton" 
		   + fs::extension(maskPath));
    }
    
    if ( verbose ) {
      std::cout << "# Output filename: " 
		<< outMaskPath.file_string() << std::endl;
    }
    
    // write output file
    writer->SetInput(thinningFilter->GetOutput());
    writer->SetFileName(outMaskPath.file_string());
    writer->SetUseCompression(true);
    writer->Update();
    
  } catch( const std::exception &e ) {  // catch any exceptions
    
    std::cerr << "Error writing output image: " << std::endl 
	      << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  
  /*******************************/
  /** End of program            **/
  /*******************************/
  
  return EXIT_SUCCESS; 
}
