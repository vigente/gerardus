# http://www.boost.org/doc/libs/release/more/getting_started/windows.html#library-naming

#=============================================================================
# Copyright (C) 2012 Daniel Pfeifer <daniel@pfeifer-mail.de>
#
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at
#   http://www.boost.org/LICENSE_1_0.txt
#=============================================================================

function(boost_library_naming)
  # TODO: make this an option
  set(BUILD_MULTI_THREADED ON)

  # Toolset detection.
  if (MSVC60)
    set(BOOST_TOOLSET "vc6")
  elseif(MSVC70)
    set(BOOST_TOOLSET "vc7")
  elseif(MSVC71)
    set(BOOST_TOOLSET "vc71")
  elseif(MSVC80)
    set(BOOST_TOOLSET "vc80")
  elseif(MSVC90)
    set(BOOST_TOOLSET "vc90")
  elseif(MSVC10)
    set(BOOST_TOOLSET "vc100")
  elseif(MSVC)
    set(BOOST_TOOLSET "vc")
  elseif(BORLAND)
    set(BOOST_TOOLSET "bcb")
  else()
    string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" BOOST_TOOLSET)
    if(CMAKE_CXX_COMPILER_VERSION MATCHES "([0-9]+)[.]([0-9]+)")
      set(BOOST_TOOLSET "${BOOST_TOOLSET}${CMAKE_MATCH_1}${CMAKE_MATCH_2}")
    endif()
  endif()

  # Append the Boost version number to the versioned name
  string(REPLACE "." "_" boost_version "${BoostCore_VERSION}")

  # The versioned name starts with the full Boost toolset
  if(WIN32)
    set(tag_toolset "-${BOOST_TOOLSET}")
    set(tag_version "-${boost_version}")
  else(WIN32)
    set(tag_toolset "")
    set(tag_version "")
  endif(WIN32)

  # Add -mt for multi-threaded libraries
  if(BUILD_MULTI_THREADED)
    set(tag_mt "-mt")
  else(BUILD_MULTI_THREADED)
    set(tag_mt "")
  endif(BUILD_MULTI_THREADED)

  # Using the debug version of the runtime library.
  # With Visual C++, this comes automatically with debug
  if(MSVC)
    set(tag_rtdebug "g")
  else(MSVC)
    set(tag_rtdebug "")
  endif(MSVC)

  # CMAKE_<CONFIG>_POSTFIX
  set_target_properties(${ARGN} PROPERTIES
    DEBUG_POSTFIX   "${tag_toolset}${tag_mt}-${tag_rtdebug}d${tag_version}"
    RELEASE_POSTFIX "${tag_toolset}${tag_mt}${tag_version}"
    )
endfunction(boost_library_naming)
