PROJECT(Clarity)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
IF(COMMAND cmake_policy)
  CMAKE_POLICY(SET CMP0003 OLD)
ENDIF(COMMAND cmake_policy)

OPTION(BUILD_SHARED_LIBS "Build Clarity with shared libraries." ON)
OPTION(BUILD_WITH_OPENMP "Build with OpenMP." ON)
OPTION(BUILD_WITH_CUDA "Build with CUDA FFT acceleration." OFF)
OPTION(TIME_MAP "Time mapping operations." OFF)
OPTION(TIME_REDUCE "Time reduction operations." OFF)
OPTION(TIME_TRANSFER "Time data transfer to/from GPU." OFF)
OPTION(TIME_FFT "Time FFT forward and inverse operations." OFF)
OPTION(TIME_DECONVOLVE "Time deconvolution operation only." OFF)
OPTION(TIME_TOTAL "Time total deconvolution operation (includes image padding and clipping." OFF)

# ramc: if CMAKE_MODULE_PATH has been set from the root CMakeLists.txt
# of Gerardus, do not overwrite it
IF(NOT CMAKE_MODULE_PATH)
  SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/")
ENDIF(NOT CMAKE_MODULE_PATH)

# ramc: operations when compiling CUDA variant of Clarity
IF(BUILD_WITH_CUDA)
  # Load CUDA macros
  FIND_PACKAGE(CUDA)
  IF(CMAKE_COMPILER_IS_GNUCC)
    # get version of compiler
    EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER} -dumpversion
      OUTPUT_VARIABLE GCC_VERSION)
    IF(GCC_VERSION VERSION_LESS 4.5)
      # ramc: CUDA v4.0 RC2 is incompatible with gcc >= 4.5. Gcc
      # 4.4. is selected in the root gerardus CMakeLists.txt, but this
      # flag needs to be set by hand, because otherwise nvcc will read
      # the version of the default gcc in the system, that is v4.5
      SET(CUDA_NVCC_FLAGS "-D__GNUC_MINOR__=4")
    ENDIF(GCC_VERSION VERSION_LESS 4.5)
  ENDIF(CMAKE_COMPILER_IS_GNUCC)
ENDIF(BUILD_WITH_CUDA)

# add FFTW library
FIND_PACKAGE(FFTW)

SET(Clarity_MAJOR_VERSION "1")
SET(Clarity_MINOR_VERSION "1")
SET(Clarity_BUILD_VERSION "0")

SET( Clarity_VERSION_FLAT "${Clarity_MAJOR_VERSION}${Clarity_MINOR_VERSION}${Clarity_BUILD_VERSION}" )
SET( Clarity_VERSION "${Clarity_MAJOR_VERSION}.${Clarity_MINOR_VERSION}.${Clarity_BUILD_VERSION}" )

INCLUDE_DIRECTORIES(${Clarity_BINARY_DIR})
INCLUDE_DIRECTORIES(${Clarity_SOURCE_DIR}
  ${Clarity_SOURCE_DIR}/include
  ${Clarity_SOURCE_DIR}/include/internal
  ${FFTW_INCLUDE_DIR}
  ${CUDA_INCLUDE}
  )

SET( Clarity_SRCS
  src/Clarity.cxx
  src/ComputePrimitives.cxx
  src/Convolve.cxx
  src/FFT.cxx
  src/ImageClip.cxx
  src/ImagePadSpatialShift.cxx
  src/JansenVanCittertDeconvolve.cxx
  src/Memory.cxx
  src/MaximumLikelihoodDeconvolve.cxx
  src/Stopwatch.cxx
  src/WienerDeconvolve.cxx
  )

CONFIGURE_FILE(
  ${Clarity_SOURCE_DIR}/Clarity.h.in
  ${Clarity_BINARY_DIR}/Clarity.h
)

IF(BUILD_WITH_CUDA)
  SET (Clarity_CUDA_SRCS
    src/ComputePrimitivesGPU.cu
    src/ModulateGPU.cu
    src/JansenVanCittertDeconvolveGPU.cu
    src/WienerDeconvolveGPU.cu
    )
  CUDA_ADD_LIBRARY( Clarity ${Clarity_SRCS} ${Clarity_CUDA_SRCS})
  CUDA_ADD_CUFFT_TO_TARGET( Clarity )
ELSE(BUILD_WITH_CUDA)
  INCLUDE_DIRECTORIES(${Clarity_SOURCE_DIR}/include)
  ADD_LIBRARY( Clarity ${Clarity_SRCS})
ENDIF(BUILD_WITH_CUDA)

TARGET_LINK_LIBRARIES(Clarity ${FFTW_LIBRARIES} ${FFTW_EXECUTABLE_LIBRARIES})

# Enable linking with OpenMP when compiling with gcc
SET (USE_GOMP ${CMAKE_COMPILER_IS_GNUCXX}
     AND ${BUILD_WITH_OPENMP} STREQUAL "ON")
IF (USE_GOMP)
  FIND_PACKAGE(Threads)
  TARGET_LINK_LIBRARIES(Clarity gomp ${CMAKE_THREAD_LIBS_INIT})
ENDIF (USE_GOMP)


SET_TARGET_PROPERTIES(Clarity PROPERTIES
  VERSION ${Clarity_VERSION}
  SOVERSION ${Clarity_MAJOR_VERSION}.${Clarity_MINOR_VERSION}
)

ADD_SUBDIRECTORY(apps)
ADD_SUBDIRECTORY(tests)

################################################################
## installation of library
################################################################
INSTALL(TARGETS Clarity
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
