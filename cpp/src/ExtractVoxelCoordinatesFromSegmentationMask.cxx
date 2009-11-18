/*
 * ExtractVoxelCoordinatesFromSegmentationMask.cxx
 *
 * Program that finds the voxels set to true in a segmentation mask,
 * and outputs a matrix with their coordinates
 * 
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
#include "itkPointSet.h"
#include "itkImageFileWriter.h"

// entry point for the program
int main(int argc, char** argv)
{
	
	/*************************************/
	/** Types and variables definitions **/
	/*************************************/

	static const unsigned int	Dimension = 3; // volume data dimension (i.e. 3D volumes)
	static const unsigned int   MatlabPrecision = 15; // number of decimal figures after the point in Matlab
	typedef double 				TScalarType; // data type for scalars
	typedef bool														BinaryPixelType;
	typedef itk::Image< BinaryPixelType, 
						Dimension >										BinaryImageType;
	typedef BinaryImageType::SizeType									BinarySizeType;
	typedef itk::ImageFileReader< BinaryImageType >						BinaryReaderType;
	typedef itk::ImageRegionConstIterator< BinaryImageType >			ConstBinaryIteratorType;
	typedef itk::PointSet< TScalarType, 
							Dimension >									PointSetType;
	typedef PointSetType::PointsContainer								PointsContainer;
	typedef PointSetType::PointType										PointType;
	typedef itk::Index< Dimension >										IndexType;

	// command line input argument types and variables
	fs::path maskPath;
	bool verbose;
	bool coordsAsIndex;
	
	// landmark I/O variables
	BinaryReaderType::Pointer 				maskReader;
	
	/*******************************/
	/** Command line parser block **/
	/*******************************/
	
	try {
		
		// Define the command line object, program description message, separator, version
		TCLAP::CmdLine cmd( "extractVoxelCoordinatesFromSegmentationMask: Extract the coordinates of voxels selected in a segmentation mask", ' ', "0.0" );
	
		// input argument: filename of input segmentation mask
		TCLAP::UnlabeledValueArg< std::string > maskPathArg( "mask", "Segmentation mask filename (binary image volume)", true, "", "file name" );
		cmd.add( maskPathArg );
		
		// input argument: coordinates format
		TCLAP::SwitchArg coordsAsIndexSwitch( "i", "index", "Format output coordinates as indices (as opposed to real world coordinates)", false );
    	cmd.add( coordsAsIndexSwitch );
	
		// input argument: verbosity
		TCLAP::SwitchArg verboseSwitch( "v", "verbose", "Increase verbosity of program output", false );
    	cmd.add( verboseSwitch );
	
		// Parse the command line arguments
		cmd.parse( argc, argv );

		// Get the value parsed by each argument
		maskPath = fs::path( maskPathArg.getValue() );
		coordsAsIndex = coordsAsIndexSwitch.getValue();
		verbose = verboseSwitch.getValue();
		
	} catch (const TCLAP::ArgException &e)  // catch any exceptions
	{
		std::cerr << "Error parsing command line: " << std::endl 
		<< e.error() << " for arg " << e.argId() << std::endl;
		return EXIT_FAILURE;
	}

	/*******************************/
	/** Load input images block   **/
	/*******************************/
		
	try {
		
		// create file readers
		maskReader = BinaryReaderType::New();
		
		// read input 3D images
		maskReader->SetFileName( maskPath.file_string() );
		if ( verbose ) {
			std::cout << "# Segmentation mask filename: " << maskPath.file_string() << std::endl;
		}
		maskReader->Update();
		

	} catch( const std::exception &e )  // catch any exceptions
	{
		std::cerr << "Error loading input landmarks masks: " << std::endl 
		<< e.what() << std::endl;
		return EXIT_FAILURE;
	}

	/******************************************************/
	/** Extract landmark coordinates from landmark mask  **/
	/******************************************************/
		
	try {
		
		// instantiate configurations of landmarks and related objects
		PointType point;
		IndexType index;
		unsigned int pointId;

		// instantiate iterators to read voxels in the binary masks
		ConstBinaryIteratorType iterator( maskReader->GetOutput(), 
										  maskReader->GetOutput()->GetLargestPossibleRegion() );
														
		// travel the binary mask
		for ( iterator.GoToBegin(), pointId = 0; !iterator.IsAtEnd(); ++iterator )
		{
			
			// if the voxel has been marked as belonging to a landmark, 
			// add its coordinates to the landmark configuration
			if ( iterator.Get() )
			{
				if ( coordsAsIndex )
				{
					index = iterator.GetIndex();
					std::cout << index[0] << ",\t" << index[1] << ",\t" << index[2] << std::endl;
				} else {
					maskReader->GetOutput()->TransformIndexToPhysicalPoint( iterator.GetIndex(), point );
					std::cout.precision( MatlabPrecision );
					std::cout << point[0] << ",\t" << point[1] << ",\t" << point[2] << std::endl;
				}
			pointId++;
			}
			
		} 
		
		if ( verbose )
		{
			BinarySizeType size;
			
			size = maskReader->GetOutput()->GetLargestPossibleRegion().GetSize();
			
			std::cout.precision( 2 );
			std::cout << "# Voxels in segmentation mask" << std::endl;
			std::cout << "#    Total:    " << size[0] * size[1] * size[2] << std::endl;
			std::cout << "#    Selected: " << pointId <<  "(" 
						<< (float)pointId 
						/ (float)( size[0] * size[1] * size[2] ) * 100.0
						<< " %)"	<< std::endl;
		}
		
	} catch( const std::exception &e )  // catch any exceptions
	{
		std::cerr << "Error extracting coordinates from segmentation mask: " << std::endl 
		<< e.what() << std::endl;
		return EXIT_FAILURE;
	}
	
	/*******************************/
	/** End of program            **/
	/*******************************/
		
	return EXIT_SUCCESS; 
}