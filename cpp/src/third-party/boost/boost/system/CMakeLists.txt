# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

## Changes by Ramon Casero <rcasero@gmail.com> for project Gerardus

cmake_minimum_required(VERSION 2.8.8)
project(BoostSystem CXX)

include(Ryppl)

ryppl_find_and_use_package(BoostConfig)
ryppl_find_and_use_package(BoostInteger) # cstdint.hpp
ryppl_find_and_use_package(BoostUtility) # assert.hpp, noncopyable.hpp, enable_if.hpp, operators.hpp

include_directories(include)

add_library(boost_system
  src/error_code.cpp
  src/local_free_on_destruction.hpp
  )

if(NOT BUILD_SHARED_LIBS)
  set_target_properties(boost_system PROPERTIES
    PREFIX "lib"
    )
endif(NOT BUILD_SHARED_LIBS)

#add_subdirectory(doc)


# ryppl_export(
#   TARGETS
#     boost_system
#   DEPENDS
#     BoostConfig
#     BoostInteger
#     BoostUtility
#   INCLUDE_DIRECTORIES
#     include
#   )

# install compiled libraries
if(WIN32)
  install(TARGETS
    boost_system
    RUNTIME ARCHIVE
    DESTINATION ${GERARDUS_SOURCE_DIR}/lib)
else(WIN32)
  install(TARGETS
    boost_system
    LIBRARY ARCHIVE
    DESTINATION ${GERARDUS_SOURCE_DIR}/lib)
endif(WIN32)
