################################################################################
# Copyright (C) 2010-2011 Daniel Pfeifer <daniel@pfeifer-mail.de>              #
#                                                                              #
# Distributed under the Boost Software License, Version 1.0.                   #
# See accompanying file LICENSE_1_0.txt or copy at                             #
#   http://www.boost.org/LICENSE_1_0.txt                                       #
################################################################################

find_program(PANDOC_EXECUTABLE
  NAMES
    pandoc
  PATHS
    $ENV{PROGRAMFILES}/Pandoc/bin
  DOC
    "a universal document converter"
  )

if(PANDOC_EXECUTABLE)
  execute_process(COMMAND ${PANDOC_EXECUTABLE} --version
    OUTPUT_VARIABLE PANDOC_VERSION
    )
  string(REGEX REPLACE "^pandoc ([.0-9]+).*" "\\1"
    PANDOC_VERSION "${PANDOC_VERSION}"
    )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PANDOC
  REQUIRED_VARS PANDOC_EXECUTABLE
  VERSION_VAR PANDOC_VERSION
  )

mark_as_advanced(PANDOC_EXECUTABLE)
