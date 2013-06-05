# Add a subdirectory for tests, documentation or examples.
#
#   ryppl_add_test_subdirectory(source_dir [binary_dir] [EXCLUDE_FROM_ALL])
#   ryppl_add_doc_subdirectory(source_dir [binary_dir] [EXCLUDE_FROM_ALL])
#   ryppl_add_example_subdirectory(source_dir [binary_dir] [EXCLUDE_FROM_ALL])
#
# Tests, documentation, and examples may have additional build dependencies.
# Ryppl needs a standardized way to disable them explicitely.
#
# ryppl_add_test_subdirectory() will call add_subdirectory() unless either
# RYPPL_DISABLE_TESTS or <project>_DISABLE_TESTS is true.
#
# ryppl_add_doc_subdirectory() will call add_subdirectory() unless either
# RYPPL_DISABLE_DOCS or <project>_DISABLE_DOCS is true.
#
# ryppl_add_example_subdirectory() will call add_subdirectory() unless either
# RYPPL_DISABLE_EXAMPLES or <project>_DISABLE_EXAMPLES is true.

#=============================================================================
# Copyright (C) 2012 Daniel Pfeifer <daniel@pfeifer-mail.de>
#
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at
#   http://www.boost.org/LICENSE_1_0.txt
#=============================================================================

macro(ryppl_add_test_subdirectory)
  if(NOT RYPPL_DISABLE_TESTS AND NOT ${PROJECT_NAME}_DISABLE_TESTS)
    add_subdirectory(${ARGV})
  endif()
endmacro(ryppl_add_test_subdirectory)

macro(ryppl_add_doc_subdirectory)
  if(NOT RYPPL_DISABLE_DOCS AND NOT ${PROJECT_NAME}_DISABLE_DOCS)
    add_subdirectory(${ARGV})
  endif()
endmacro(ryppl_add_doc_subdirectory)

macro(ryppl_add_example_subdirectory)
  if(NOT RYPPL_DISABLE_EXAMPLES AND NOT ${PROJECT_NAME}_DISABLE_EXAMPLES)
    add_subdirectory(${ARGV})
  endif()
endmacro(ryppl_add_example_subdirectory)
