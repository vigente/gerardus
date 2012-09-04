################################################################################
# Copyright (C) 2010-2011 Daniel Pfeifer <daniel@pfeifer-mail.de>              #
#                                                                              #
# Distributed under the Boost Software License, Version 1.0.                   #
# See accompanying file LICENSE_1_0.txt or copy at                             #
#   http://www.boost.org/LICENSE_1_0.txt                                       #
################################################################################

find_program(FO_PROCESSOR
  NAMES
    fop
    xep.bat
  PATHS
    $ENV{PROGRAMFILES}/RenderX/XEP
  DOC
    "An XSL-FO processor"
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FOProcessor DEFAULT_MSG FO_PROCESSOR)
