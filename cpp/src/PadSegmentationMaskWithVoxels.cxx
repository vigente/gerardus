/*
 * PadSegmentationMaskWithVoxels.cxx
 *
 * Program to pad the sides of a 3D segmentation mask (e.g. produced by Seg3D) with new voxels
 * 
 * Example of usage:
 * 
 * $ ./padSegmentationMaskWithVoxels mask.mha 339 225 244 365 39 111
 * 
 * This will pad the volume in mask.mha with 339 voxeks to the left and 225 to the right, in the 
 * X-axis. Also, 244 voxels to the front and 365 to the back of the volume (Y-axis). Finally,
 * 39 voxels to the bottom and 111 voxels to the top of the volume (Z-axis).
 * 
 */ 
 
 /*
  * Author: Ramón Casero <rcasero@gmail.com>
  * Copyright © 2009 University of Oxford
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
#include "itkPasteImageFilter.h"
#include "itkImageRegionIterator.h"

// entry point for the program
int main(int argc, char** argv)
{
	
	/*************************************/
	/** Types and variables definitions **/
	/*************************************/

	static const unsigned int	Dimension = 3; // volume data dimension (i.e. 3D volumes)
	typedef double 				TScalarType; // data type for scalars
	typedef itk::Index< Dimension >										IndexType;
	
	
	typedef unsigned char												UCharPixelType;
	typedef itk::Image< UCharPixelType, 
						Dimension >										UCharImageType;
	typedef UCharImageType::SizeType									UCharSizeType;
	typedef itk::ImageFileReader< UCharImageType >						UCharReaderType;

	typedef itk::PasteImageFilter< UCharImageType, 
									UCharImageType, 
									UCharImageType >					PasteImageFilterType;

	typedef itk::ImageFileWriter< UCharImageType > 					UCharWriterType;
	
	typedef itk::ImageRegionIterator< UCharImageType > 				IteratorType;

	// command line input argument types and variables
	fs::path maskPath;
	bool verbose;
	fs::path outMaskPath;
	
	// landmark I/O variables
	UCharReaderType::Pointer 				maskReader;
	
	// image variables
	UCharSizeType							sizeIn, sizeOut; // size of input and output images
	size_t 									padXa, padXb, padYa, padYb, padZa, padZb; // number of voxels to pad on each side
	UCharImageType::Pointer					imIn, imOut; // input and padded images
	UCharImageType::RegionType				regionIn, regionOut; // intermediate variable to create output image
	UCharImageType::IndexType				startIn, startOut; // start corner of output image
	PasteImageFilterType::Pointer			pasteFilter; // filter to past input image onto output padded image
	
	/*******************************/
	/** Command line parser block **/
	/*******************************/
	
	try {
		
		// Define the command line object, program description message, separator, version
		TCLAP::CmdLine cmd( "padSegmentationMaskWithVoxels: Pad the sides of a segmentation mask with new voxels", ' ', "0.0" );
	
		// input argument: filename of input segmentation mask
		TCLAP::UnlabeledValueArg< std::string > maskPathArg( "mask", "Segmentation mask filename (binary image volume)", true, "", "file" );
		cmd.add( maskPathArg );

		// input argument: number of padding voxels. Each axis has 2 values, for both sides of the volume
		TCLAP::UnlabeledValueArg< size_t > padXaArg( "xa", "Number of padding voxels, X axis, before volume", true, 0, "xa" );
		TCLAP::UnlabeledValueArg< size_t > padXbArg( "xb", "Number of padding voxels, X axis, after volume", true, 0, "xb" );
		cmd.add( padXaArg );		
		cmd.add( padXbArg );		
		TCLAP::UnlabeledValueArg< size_t > padYaArg( "ya", "Number of padding voxels, Y axis, before volume", true, 0, "ya" );
		TCLAP::UnlabeledValueArg< size_t > padYbArg( "yb", "Number of padding voxels, Y axis, after volume", true, 0, "yb" );
		cmd.add( padYaArg );		
		cmd.add( padYbArg );		
		TCLAP::UnlabeledValueArg< size_t > padZaArg( "za", "Number of padding voxels, Z axis, before volume", true, 0, "za" );
		TCLAP::UnlabeledValueArg< size_t > padZbArg( "zb", "Number of padding voxels, Z axis, after volume", true, 0, "zb" );
		cmd.add( padZaArg );		
		cmd.add( padZbArg );		
	
		// input argument: filename of output segmentation mask
		TCLAP::ValueArg< std::string > outMaskPathArg( "o", "outfile", "Output mask filename (binary image volume)", false, "", "file" );
		cmd.add( outMaskPathArg );

		// input argument: verbosity
		TCLAP::SwitchArg verboseSwitch( "v", "verbose", "Increase verbosity of program output", false );
    	cmd.add( verboseSwitch );
	
		// Parse the command line arguments
		cmd.parse( argc, argv );

		// Get the value parsed by each argument
		maskPath = fs::path( maskPathArg.getValue() );
		outMaskPath = fs::path( outMaskPathArg.getValue() );
		verbose = verboseSwitch.getValue();
		padXa = padXaArg.getValue();
		padXb = padXbArg.getValue();
		padYa = padYaArg.getValue();
		padYb = padYbArg.getValue();
		padZa = padZaArg.getValue();
		padZb = padZbArg.getValue();
		
	} catch (const TCLAP::ArgException &e)  // catch any exceptions
	{
		std::cerr << "Error parsing command line: " << std::endl 
		<< e.error() << " for arg " << e.argId() << std::endl;
		return EXIT_FAILURE;
	}
	
	/*******************************/
	/** Load input image block    **/
	/*******************************/
		
	try {
		
		// create file readers
		maskReader = UCharReaderType::New();
		
		// read input 3D images
		maskReader->SetFileName( maskPath.file_string() );
		if ( verbose ) {
			std::cout << "# Segmentation mask filename: " << maskPath.file_string() << std::endl;
		}
		maskReader->Update();
		
		// get input image
		imIn = maskReader->GetOutput();
		
		// get mask's size
		sizeIn = imIn->GetLargestPossibleRegion().GetSize();

		if ( verbose ) {
			std::cout << "# Mask dimensions: " << sizeIn[0] << "\t" << sizeIn[1] << "\t" << sizeIn[2] << std::endl; 
		}

	} catch( const std::exception &e )  // catch any exceptions
	{
		std::cerr << "Error loading input segmentation mask: " << std::endl 
		<< e.what() << std::endl;
		return EXIT_FAILURE;
	}

	/*******************************/
	/** Pad input image block     **/
	/*******************************/
		
	try {
		
		// create an output image that will contain the padded image
		imOut = UCharImageType::New();
		
		// make output image big enough to contain the padding
		sizeOut[0] = sizeIn[0] + padXa + padXb;
		sizeOut[1] = sizeIn[1] + padYa + padYb; 
		sizeOut[2] = sizeIn[2] + padZa + padZb;
		regionOut.SetSize( sizeOut );

		// set starting corner of the output image
		startIn = startOut = imIn->GetLargestPossibleRegion().GetIndex();
		startOut[0] -= padXa;
		startOut[1] -= padYa;
		startOut[2] -= padZa;
		
		regionOut.SetIndex( startOut );
		
		// copy the metainformation from the input image
		imOut->CopyInformation( imIn );

		// finish creating the output image 		
		imOut->SetRegions( regionOut );
		imOut->Allocate();
		imOut->Update();
		
		// set output image region the same as input image
		regionOut.SetSize( sizeIn );
		regionOut.SetIndex( startIn );
		imOut->SetRequestedRegion( regionOut );
		
		// iterate input image and copy values to selected output image region
		IteratorType itIn( imIn, imIn->GetLargestPossibleRegion() );
		IteratorType itOut( imOut, imOut->GetRequestedRegion() );
	
  		
  		for ( itIn = itIn.Begin(), itOut = itOut.Begin(); !itIn.IsAtEnd(); ++itIn, ++itOut ) 
  		{
  			itOut.Set( itIn.Get() );
  		}
		
		imOut->Update();
		
		
	} catch( const std::exception &e )  // catch any exceptions
	{
		std::cerr << "Error padding input segmentation mask: " << std::endl 
		<< e.what() << std::endl;
		return EXIT_FAILURE;
	}
		
	
	/*******************************/
	/** Output block              **/
	/*******************************/
		
	try {
		
		// create writer for the results file
		UCharWriterType::Pointer writer = UCharWriterType::New();
		
		// create a filename for the registered image by appending 
		// "padded" to the input image filename, if none is
		// provided explicitely in the command line
		if ( outMaskPath.empty() ) {
			outMaskPath = maskPath.branch_path() 
			/ fs::path( fs::basename( maskPath ) + "-padded" 
			+ fs::extension( maskPath ) );
		}

		if ( verbose ) {
			std::cout << "# Output filename: " << outMaskPath.file_string() << std::endl;
		}

		// write output file
		writer->SetInput( imOut );
		writer->SetFileName( outMaskPath.file_string() );
		writer->SetUseCompression( true );
		writer->Update();
		
	} catch (const std::exception &e)  // catch any exceptions
	{
		std::cerr << "Error saving results: " << std::endl 
		<< e.what() << std::endl;
		return EXIT_FAILURE;
	}
		
	/*******************************/
	/** End of program            **/
	/*******************************/
	
	return EXIT_SUCCESS; 
	
}
