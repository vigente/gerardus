# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostProgramOptions CXX)

include(Ryppl)

ryppl_find_and_use_package(BoostAny)
ryppl_find_and_use_package(BoostBind)
ryppl_find_and_use_package(BoostConfig)
ryppl_find_and_use_package(BoostConversion)
ryppl_find_and_use_package(BoostCore)
ryppl_find_and_use_package(BoostException)
ryppl_find_and_use_package(BoostFunction)
ryppl_find_and_use_package(BoostIterator)
ryppl_find_and_use_package(BoostSmartPtr)
ryppl_find_and_use_package(BoostStaticAssert)
ryppl_find_and_use_package(BoostTokenizer)
ryppl_find_and_use_package(BoostTypeTraits)
ryppl_find_and_use_package(BoostUtility)

include_directories(include)

add_library(boost_program_options
  src/cmdline
  src/config_file
  src/options_description
  src/parsers
  src/variables_map
  src/value_semantic
  src/positional_options
  src/utf8_codecvt_facet
  src/convert
  src/winmain
  src/split
  )

# ryppl_export(
#   TARGETS
#     boost_program_options
#   DEPENDS
#     BoostAny
#     BoostConfig
#     BoostConversion
#     BoostCore
#     BoostException
#     BoostFunction
#     BoostIterator
#     BoostSmartPtr
#     BoostStaticAssert
#     BoostTypeTraits
#     BoostUtility
#   INCLUDE_DIRECTORIES
#     include
#   )
