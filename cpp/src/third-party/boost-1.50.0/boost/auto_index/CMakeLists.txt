# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostAutoIndex CXX)

include(Ryppl)

ryppl_find_and_use_package(BoostArray)
ryppl_find_and_use_package(BoostFilesystem)
ryppl_find_and_use_package(BoostProgramOptions)

add_executable(auto_index
  src/auto_index.cpp
  src/file_scanning.cpp
  src/index_generator.cpp
  src/tiny_xml.cpp
  )

target_link_libraries(auto_index
  boost_program_options
  boost_filesystem
  boost_system
  boost_regex
  )

# ryppl_export(TARGETS auto_index)
