# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(Quickbook CXX)

include(Ryppl)

ryppl_find_and_use_package(BoostFilesystem)
ryppl_find_and_use_package(BoostProgramOptions)
ryppl_find_and_use_package(BoostSpirit)

if(MSVC)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS
    /wd4244 /wd4267 /wd4355 /wd4511 /wd4512 /wd4701 /wd4702 /wd4800)
endif()

add_executable(quickbook
  src/quickbook.cpp
  src/actions.cpp
  src/doc_info_actions.cpp
  src/state.cpp
  src/utils.cpp
  src/files.cpp
  src/string_ref.cpp
  src/input_path.cpp
  src/values.cpp
  src/id_manager.cpp
  src/post_process.cpp
  src/collector.cpp
  src/template_stack.cpp
  src/code_snippet.cpp
  src/markups.cpp
  src/syntax_highlight.cpp
  src/grammar.cpp
  src/main_grammar.cpp
  src/block_element_grammar.cpp
  src/phrase_element_grammar.cpp
  src/doc_info_grammar.cpp
  )

target_link_libraries(quickbook
  boost_program_options
  boost_filesystem
  boost_system
  )

ryppl_export(TARGETS quickbook)
