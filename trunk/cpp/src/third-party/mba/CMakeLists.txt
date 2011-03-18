CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# MBA library
#FIND_PACKAGE(OpenGL REQUIRED)
ADD_LIBRARY(MBA SHARED
src/filterRestrict2D.cpp
src/MBA.cpp
src/MBAdata.cpp
src/PointAccessUtils.cpp
src/UCBsplines.cpp
src/UCBsplineSurface.cpp
src/UCButils.cpp)
ADD_DEFINITIONS("-O2 -Wall")
INCLUDE_DIRECTORIES(include)
#TARGET_LINK_LIBRARIES(MBA glut GLU GL)
# this installs in /usr/lib/.
#SET(CMAKE_INSTALL_PREFIX .)
#INSTALL(TARGETS MBA 
#LIBRARY DESTINATION .)
