# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostDisjointSets)

include(Ryppl)

#ryppl_test_suite(
#  RUN
#    disjoint_set_test.cpp
#  )


ryppl_export(
  DEPENDS
    BoostGraph
  INCLUDE_DIRECTORIES
    include
  )
