# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostCompatibility)

#include(BoostDocumentation)
include(Ryppl)

#boost_documentation(index.html)

ryppl_export(
  INCLUDE_DIRECTORIES
    include
  )
