ADD_EXECUTABLE(DeconvolveRawImages DeconvolveRawImages.cxx)

TARGET_LINK_LIBRARIES(DeconvolveRawImages
  Clarity
  ${FOUND_CUFFT}
)
