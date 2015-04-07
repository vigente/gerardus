

# Important note #

This page is obsolete. The architecture of itk\_imfilter() has completely changed since it was written.

# Introduction #

This page provides step-by-step instructions on how to add a new ITK filter to the Gerardus Matlab MEX function itk\_imfilter().

If the automatic script wants to be used, then it's necessary to use bash on Linux or cygwin on Windows.

All filters must inherit from class [itk::ImageToImageFilter](http://www.itk.org/Doxygen320/html/classitk_1_1ImageToImageFilter.html).

# Common steps to all filters #

  1. Make the filter's C++ code available to Gerardus
    * If the filter is already part of ITK (whether the standard or Review packages), then Gerardus can already see the code, and you can jump to the next step
    * If the filter comes from an external source, e.g. the [itk::AnisotropicDiffusionVesselEnhancementImageFilter](http://hdl.handle.net/1926/558) from the Insight Journal, then
      1. Put the filter's C++ files in a new subdirectory within `gerardus/cpp/src/third-party`, e.g.
```
gerardus/cpp/src/third-party/IJ-Vessel_Enhancement_Diffusion.1
```
      1. Add the new subdirectory to the subversion repository, add a reference in the `ChangeLog` file to where the code was downloaded from, and commit it, e.g.
```
$ cd gerardus
$ svn add cpp/src/third-party/IJ-Vessel_Enhancement_Diffusion.1
$ svn propset svn:eol-style native cpp/src/third-party/IJ-Vessel_Enhancement_Diffusion.1/*
$ emacs ChangeLog
$ svn ci
```
  1. Run bash script [add\_filter\_template.sh](http://code.google.com/p/gerardus/source/browse/trunk/matlab/ItkToolbox/add_filter_template.sh) as described below to get the new filter into `itk_imfilter()`:
```
$ cd gerardus/matlab/ItkToolbox
$ ./add_filter_template.sh AnisotropicDiffusionVesselEnhancementImageFilter advess
```
> > where `AnisotropicDiffusionVesselEnhancementImageFilter` is the name of the filter, that is assumed to correspond to an `itk::AnisotropicDiffusionVesselEnhancementImageFilter` class that inherits from `itk::ImageToImageFilter`; and `advess` is the short name of the filter. That is, from Matlab the filter can be invoked running
```
>> im2 = itk_imfilter('advess', im);
```

The script `add_filter_template.sh` automates the steps for "vanilla" filters that:
  * return the same output image type as the input image type
  * don't require user-provided parameters (e.g. radius, threhold, etc.)
  * don't require specific setup steps

Below we describe the manual changes that are necessary for more complicated filters.

## If filter accepts user-provided parameters ##

For example, filter `AnisotropicDiffusionVesselEnhancement` accepts, amongst others, the following user-provided parameters according to the paper and the code:
  * Minimum sigma (`double`, default 0.2)
  * Maximum sigma (`double`, default 2.0)
  * Number of sigma steps (`int`, default 10)

  1. Edit the .hpp file automatically created for the new filter, [MexAnisotropicDiffusionVesselEnhancementImageFilter.hpp](http://code.google.com/p/gerardus/source/browse/trunk/matlab/ItkToolbox/MexAnisotropicDiffusionVesselEnhancementImageFilter.hpp). Add accepted user-provided parameters as member variables
```
protected:

  // user-provided input parameters
  double sigmaMin;
  double sigmaMax;
  int    numSigmaSteps;
```
  1. Edit the .cpp file automatically created, [MexAnisotropicDiffusionVesselEnhancementImageFilter.cpp](http://code.google.com/p/gerardus/source/browse/trunk/matlab/ItkToolbox/MexAnisotropicDiffusionVesselEnhancementImageFilter.cpp).
    1. Correct the number of input parameters that are accepted
```
  if (this->nparam > 3) {
    mexErrMsgTxt("Too many input arguments");
  }
```
    1. read the parameters from the input
```
  // get user-provided parameters: 
  //    parameter name
  //    index (0 = first parameter)
  //    default value
  this->sigmaMin = this->template
    GetScalarParamValue<double>("SIGMAMIN",    0, 0.2);
  this->sigmaMax = this->template
    GetScalarParamValue<double>("SIGMAMAX",    1, 2.0);
  this->numSigmaSteps = this->template
    GetScalarParamValue<int>("NUMSIGMASTEPS",  2, 10);
```
  1. Go back to the .hpp file, and uncomment the declaration of `FilterAdvancedSetup()`
```
  // if this particular filter needs to redefine one or more MexBaseFilter
  // virtual methods, the corresponding declarations go here
  // void CheckNumberOfOutputs();
  void FilterAdvancedSetup();
  // void ExportOtherFilterOutputsToMatlab();
```
  1. Go back again to the .cpp file, uncomment the template definition of `FilterAdvancedSetup()`, and add the necessary code to pass the parameters to the filter. Note that because these methods are specific to each filter, it's not possible to call them using `filter->SetXXX()`. Instead, it's necessary to use the provided local pointer `localFilter->SetXXX()`
```
// add code here if you need to pass user-provided parameters to the
// filter, or perform any other kind of filter setup
template <class InVoxelType, class OutVoxelType>
void MexAnisotropicDiffusionVesselEnhancementImageFilter<InVoxelType, 
			    OutVoxelType>::FilterAdvancedSetup() {
  
  // create a local pointer to the filter so that we can use
  // methods that are not part of the MexBaseFilter
  typename FilterType::Pointer localFilter = 
    dynamic_cast<typename MexAnisotropicDiffusionVesselEnhancementImageFilter<InVoxelType,
				 OutVoxelType>::FilterType *>(this->filter.GetPointer());

  // set user-provided parameters
  localFilter->SetSigmaMin(this->sigmaMin);
  localFilter->SetSigmaMax(this->sigmaMax);
  localFilter->SetNumberOfSigmaSteps(this->numSigmaSteps);

}
```