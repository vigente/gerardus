# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostCore)

include(Ryppl)

set(Boost_VERSION_MAJOR  1)
set(Boost_VERSION_MINOR 50)
set(Boost_VERSION_PATCH  0)

math(EXPR BOOST_VERSION
  "${Boost_VERSION_MAJOR} * 100000 + ${Boost_VERSION_MINOR} * 100 + ${Boost_VERSION_PATCH}"
  )
string(REPLACE "." "_" BOOST_LIB_VERSION "${Boost_VERSION}")

configure_file(include/boost/version.hpp.in include/boost/version.hpp)

#ryppl_add_test_subdirectory(test)

ryppl_export(
  DEPENDS
    BoostConfig
#    BoostInteger
#    BoostIterator
#    BoostMPL
#    BoostPreprocessor
#    BoostSmartPtr
#    BoostStaticAssert
#    BoostTypeTraits
#    BoostUtility
  INCLUDE_DIRECTORIES
    "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_BINARY_DIR}/include"
  VERSION "1.50"
  )
