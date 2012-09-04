# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostInspect CXX)

include(Ryppl)

ryppl_find_and_use_package(BoostAlgorithm)
ryppl_find_and_use_package(BoostConfig)
ryppl_find_and_use_package(BoostConversion)
ryppl_find_and_use_package(BoostFilesystem)
ryppl_find_and_use_package(BoostRegex)
ryppl_find_and_use_package(BoostSmartPtr)
ryppl_find_and_use_package(BoostTest)
ryppl_find_and_use_package(BoostUtility)

add_executable(inspect
  apple_macro_check.cpp
  ascii_check.cpp
  assert_macro_check.cpp
  copyright_check.cpp
  crlf_check.cpp
  end_check.cpp
  inspect.cpp
  license_check.cpp
  link_check.cpp
  minmax_check.cpp
  path_name_check.cpp
  tab_check.cpp
  unnamed_namespace_check.cpp
  )

target_link_libraries(inspect
  boost_filesystem
  boost_system
  boost_regex
  )

# ryppl_export(TARGETS inspect)
