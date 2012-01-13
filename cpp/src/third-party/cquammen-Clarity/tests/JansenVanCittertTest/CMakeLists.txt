ADD_EXECUTABLE(JansenVanCittertTest JansenVanCittertTest.cxx)

TARGET_LINK_LIBRARIES(JansenVanCittertTest
  Clarity
  ${FOUND_CUFFT}
)
