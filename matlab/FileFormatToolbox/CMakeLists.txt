# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2014 University of Oxford
# Version: 0.1.1
# $Rev$
# $Date$
#
# University of Oxford means the Chancellor, Masters and Scholars of
# the University of Oxford, having an administrative office at
# Wellington Square, Oxford OX1 2JD, UK. 
#
# This file is part of Gerardus.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details. The offer of this
# program under the terms of the License is subject to the License
# being interpreted in accordance with English Law and subject to any
# action against the University of Oxford being under the jurisdiction
# of the English Courts.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 2.8)

# compilation of Matlab-code functions as executables. 
if(NOT WIN32)

  # This approach works, but makes CMake build the function every time
  # "make" is invoked, regardless of whether the source code has
  # changed
  #
  # # vmu2png
  # add_custom_command(
  #   OUTPUT  run_vmu2png.sh vmu2png
  #   COMMAND ${MATLAB_ROOT}/bin/mcc
  #   ARGS    -m vmu2png.m
  #   DEPENDS vmu2png.m
  #   COMMENT "Compiling vmu2png.m into a binary"
  #   )
  # add_custom_target(do_vmu2png ALL
  #   DEPENDS run_vmu2png.sh vmu2png
  #   SOURCES vmu2png.m
  #   )

  # This approach compiles the Matlab function at configuration rather
  # than build time, but it's a simple way to avoid building every
  # time we run "make" (see comment above). On the other hand, if
  # vmu2png.m changes, we'll need to manually delete run_vmu2png.sh to
  # trigger another compilation
  if (NOT EXISTS "${GERARDUS_MATLAB_SOURCE_DIR}/FileFormatToolbox/run_vmu2png.sh")
    message(STATUS "Compiling Matlab function vmu2png.m...")
    execute_process(
      COMMAND "${MATLAB_ROOT}/bin/mcc" -m vmu2png.m
      RESULT_VARIABLE result
      WORKING_DIRECTORY "${GERARDUS_MATLAB_SOURCE_DIR}/FileFormatToolbox"
      )
    message(STATUS "... compilation result: ${result}")
  else()
    message(STATUS "Matlab function vmu2png.m already compiled. Skipping compilation.")
    endif()

endif(NOT WIN32)
