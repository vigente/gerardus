ADD_EXECUTABLE(BlindMaximumLikelihoodTest BlindMaximumLikelihoodTest.cxx)

TARGET_LINK_LIBRARIES(BlindMaximumLikelihoodTest
  Clarity
  ${FOUND_CUFFT}
)
