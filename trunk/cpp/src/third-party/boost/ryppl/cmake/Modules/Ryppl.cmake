#=============================================================================
# Copyright (C) 2012 Daniel Pfeifer <daniel@pfeifer-mail.de>
#
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at
#   http://www.boost.org/LICENSE_1_0.txt
#=============================================================================

# silence warnings about unused variables
set(_ryppl_disable_warnings
    RYPPL_DISABLE_TESTS 
    RYPPL_DISABLE_DOCS 
    RYPPL_DISABLE_EXAMPLES
    BUILD_SHARED_LIBS 
    RYPPL_INITIAL_PASS
    )
foreach(varname ${_ryppl_disable_warnings})
  if(${varname})
  endif()
endforeach()

include(RypplAddSubdirectory)
include(RypplExport)
include(RypplFindPackage)
