##########################################################################
# Copyright (C) 2012 Daniel Pfeifer <daniel@pfeifer-mail.de>             #
#                                                                        #
# Distributed under the Boost Software License, Version 1.0.             #
# See accompanying file LICENSE_1_0.txt or copy at                       #
#   http://www.boost.org/LICENSE_1_0.txt                                 #
##########################################################################

if(PYTHONPATH)
  set(ENV{PYTHONPATH} ${PYTHONPATH})
endif(PYTHONPATH)

set(success 1)

separate_arguments(COMMANDS)

foreach(cmd ${COMMANDS})
  set(command "${${cmd}}")
  set(fail "${${cmd}_FAIL}")

  separate_arguments(command)
  execute_process(COMMAND ${command}
    RESULT_VARIABLE result
    OUTPUT_QUIET ERROR_QUIET # TODO: write to file
    )

  if((fail AND result EQUAL 0) OR (NOT fail AND NOT result EQUAL 0))
    set(success 0)
    break()
  endif()
endforeach(cmd)

file(WRITE ${OUTPUT} ${success})
