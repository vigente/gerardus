# This project is designed to be built outside the Insight source tree.
PROJECT(VED)

INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

ENABLE_TESTING()
INCLUDE(Dart)

SET(Required_Packages
ITK
)

FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

SET(Libraries
ITKIO
ITKBasicFilters
ITKCommon
)

SET(CurrentExe "ImageCompare")
ADD_EXECUTABLE(${CurrentExe} ImageCompare.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

ADD_EXECUTABLE( itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilterTest
                itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilterTest.cxx 
              )

TARGET_LINK_LIBRARIES( itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilterTest
  ITKIO
  ITKCommon
  ITKBasicFilters
  )

ADD_EXECUTABLE( itkAnisotropicDiffusionVesselEnhancementImageFilterTest
                itkAnisotropicDiffusionVesselEnhancementImageFilterTest.cxx 
              )

TARGET_LINK_LIBRARIES( itkAnisotropicDiffusionVesselEnhancementImageFilterTest
  ITKIO
  ITKCommon
  ITKBasicFilters
  )

ADD_TEST(MultiScaleHessianEnhancedTest
itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilterTest
${CMAKE_SOURCE_DIR}/CroppedWholeLungCTScan.mhd
${CMAKE_SOURCE_DIR}/HessianEnhanced0.5-4.0.mha 0.4 5.0 10 )

ADD_TEST(MultiScaleHessianEnhancedRegressionTest ImageCompare
${CMAKE_SOURCE_DIR}/HessianEnhanced0.5-4.0.mha
${CMAKE_SOURCE_DIR}/HessianEnhanced0.5-4.0_baseline.mha)

ADD_TEST(VEDEnhancedsTest
itkAnisotropicDiffusionVesselEnhancementImageFilterTest
${CMAKE_SOURCE_DIR}/CroppedWholeLungCTScan.mhd
${CMAKE_SOURCE_DIR}/VEDEnhanced0.5-4.0.mha 0.4 5.0 10 10)

ADD_TEST(VEDEnhancedRegressionTest ImageCompare
${CMAKE_SOURCE_DIR}/VEDEnhanced0.5-4.0.mha
${CMAKE_SOURCE_DIR}/VEDEnhanced0.5-4.0_baseline.mha)

