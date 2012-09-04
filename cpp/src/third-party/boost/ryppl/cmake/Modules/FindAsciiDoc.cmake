################################################################################
# Copyright (C) 2010-2011 Daniel Pfeifer <daniel@pfeifer-mail.de>              #
#                                                                              #
# Distributed under the Boost Software License, Version 1.0.                   #
# See accompanying file LICENSE_1_0.txt or copy at                             #
#   http://www.boost.org/LICENSE_1_0.txt                                       #
################################################################################

if(ASCIIDOC_EXECUTABLE)
  return()
endif(ASCIIDOC_EXECUTABLE)

find_program(ASCIIDOC_EXECUTABLE
  NAMES
    asciidoc
  DOC
    "Text based document generation"
  )

if(ASCIIDOC_EXECUTABLE)
  execute_process(COMMAND ${ASCIIDOC_EXECUTABLE} --version
    OUTPUT_VARIABLE ASCIIDOC_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  string(REGEX REPLACE "^asciidoc (.+)$" "\\1"
    ASCIIDOC_VERSION "${ASCIIDOC_VERSION}"
    )
endif(ASCIIDOC_EXECUTABLE)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ASCIIDOC
  REQUIRED_VARS ASCIIDOC_EXECUTABLE
  VERSION_VAR ASCIIDOC_VERSION
  )

mark_as_advanced(ASCIIDOC_EXECUTABLE)
