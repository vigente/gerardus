##########################################################################
# Copyright (C) 2012 Daniel Pfeifer <daniel@pfeifer-mail.de>             #
#                                                                        #
# Distributed under the Boost Software License, Version 1.0.             #
# See accompanying file LICENSE_1_0.txt or copy at                       #
#   http://www.boost.org/LICENSE_1_0.txt                                 #
##########################################################################

set(failed)
separate_arguments(TESTS)

foreach(test ${TESTS})
  file(READ ${test}_ok.txt ok LIMIT 1)
  if(NOT ok)
    list(APPEND failed ${test})  
  endif(NOT ok)
endforeach(test)

if(failed)
  string(REPLACE ";" "\n  " failed "${failed}")
  message(FATAL_ERROR "The following tests failed:\n  ${failed}")
endif(failed)
