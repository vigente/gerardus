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
#include "itkScaleTransform.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

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
    
    try {
        
        // Define the command line object, program description message, separator, version
        TCLAP::CmdLine cmd( "resize3DImage: resize a 3D image", ' ', "0.0" );
    
        // input argument: filename of output image
        TCLAP::ValueArg< std::string > outImPathArg("o", "outfile", "Output image filename", false, "", "file");
        cmd.add(outImPathArg);

        // input argument: interpolating type
        TCLAP::ValueArg< std::string > interpTypeArg("i", "interp", "Interpolator type: bspline (default), nn", false, "bspline", "string");
        cmd.add(interpTypeArg);

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
        interpType = interpTypeArg.getValue();
        
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

    /*******************************/
    /** Resize image              **/
    /*******************************/

    typedef unsigned short                               UShortPixelType;
    typedef itk::Image< UShortPixelType, 
                        Dimension >                      OutputImageType;
    typedef OutputImageType::SizeType                    OutputSizeType;

    typedef itk::ScaleTransform< TScalarType, 
                                  Dimension >            TransformType;
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
    OutputSizeType                                       sizeOut;
 
    InputImageType::SpacingType spacing;  
    spacing = imIn->GetSpacing();
    InputImageType::PointType    origin  = imIn->GetOrigin();

    TransformType::ScaleType                            scaling;

    // filters
    TransformType::Pointer                               transform;
    ResampleFilterType::Pointer                          resampler;
    InterpolatorType::Pointer                            interpolator;

    try {

        // create objects for rotation
        transform = TransformType::New();
        resampler = ResampleFilterType::New();
        if (interpType == "bspline") {
            interpolator = BSplineInterpolatorType::New();
        } else if (interpType == "nn") {
            interpolator = NearestNeighborInterpolatorType::New();
        } else {
            throw std::string("Invalid interpolator type");
        }
        
        // compute output size rounding to integer number of pixels
        for (size_t i = 0; i < Dimension; ++i) {
            sizeOut[i] = round((double)sizeIn[i] * 0.5); 
        }        
        
        // recompute scaling factor to take into account the rounding of size
        for (size_t i = 0; i < Dimension; ++i) { 
//            scaling[i] = (double)sizeOut[i] / (double)sizeIn[i];
//            spacing[i] /= scaling[i]; 
            scaling[i] = (double)sizeIn[i] / (double)sizeOut[i];
//            spacing[i] *= scaling[i];
            origin[i] /= scaling[i];
        }
        transform->SetScale(scaling);
        
        // set all the bits and pieces that go into the resampler
        resampler->SetInterpolator(interpolator);
        resampler->SetTransform(transform);
        resampler->SetOutputOrigin(origin);
        resampler->SetOutputSpacing(spacing);
        resampler->SetSize(sizeOut);
        resampler->SetInput(imIn); 
        
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
            std::cout << "# Output filename: " << outImPath.file_string() << std::endl;
        }
        
        // write output file
        writer->SetInput(imOut);
        writer->SetFileName(outImPath.file_string());
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
        