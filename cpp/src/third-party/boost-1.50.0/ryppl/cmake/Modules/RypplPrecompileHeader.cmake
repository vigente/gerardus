##########################################################################
# Copyright (C) 2010-2011 Daniel Pfeifer <daniel@pfeifer-mail.de>        #
#                                                                        #
# Distributed under the Boost Software License, Version 1.0.             #
# See accompanying file LICENSE_1_0.txt or copy at                       #
#   http://www.boost.org/LICENSE_1_0.txt                                 #
##########################################################################


##
function(ryppl_add_pch name source_list)
  set(pch_header "${CMAKE_CURRENT_BINARY_DIR}/${name}_pch.hpp")
  set(pch_source "${CMAKE_CURRENT_BINARY_DIR}/${name}_pch.cpp")

  file(WRITE ${pch_header}.in "/* ${name} precompiled header file */\n\n")
  foreach(header ${ARGN})
    if(header MATCHES "^<.*>$")
      file(APPEND ${pch_header}.in "#include ${header}\n")
    else()
      get_filename_component(header ${header} ABSOLUTE)
      file(APPEND ${pch_header}.in "#include \"${header}\"\n")
    endif()
  endforeach(header)
  configure_file(${pch_header}.in ${pch_header} COPYONLY)

  file(WRITE ${pch_source}.in "#include \"${pch_header}\"\n")
  configure_file(${pch_source}.in ${pch_source} COPYONLY)

  if(MSVC_IDE)
    set(pch_binary "$(IntDir)/${name}.pch")
    set_source_files_properties(${pch_source} PROPERTIES
      COMPILE_FLAGS "/Yc\"${pch_header}\" /Fp\"${pch_binary}\""
      OBJECT_OUTPUTS "${pch_binary}"
      )
    set_source_files_properties(${${source_list}} PROPERTIES
      COMPILE_FLAGS "/Yu\"${pch_header}\" /FI\"${pch_header}\" /Fp\"${pch_binary}\""
      OBJECT_DEPENDS "${pch_binary}"
      )
    set(${source_list} ${pch_source} ${${source_list}} PARENT_SCOPE)
  endif(MSVC_IDE)
  
  set(PCH_HEADER ${pch_header} PARENT_SCOPE)
endfunction(ryppl_add_pch)


##
function(ryppl_add_pch_to_target target header)
  if(NOT header OR BOOST_DISABLE_PCH)
    return()
  endif(NOT header OR BOOST_DISABLE_PCH)

  if(XCODE_VERSION)
    set_target_properties(${target} PROPERTIES
      XCODE_ATTRIBUTE_GCC_PREFIX_HEADER "${header}"
      XCODE_ATTRIBUTE_GCC_PRECOMPILE_PREFIX_HEADER "YES"
      )
    return()
  endif(XCODE_VERSION)

  if(CMAKE_COMPILER_IS_GNUCXX)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" build_type)
    set(compile_flags ${CMAKE_CXX_FLAGS_${build_type}})

    get_target_property(target_type ${target} TYPE)
    if(${target_type} STREQUAL SHARED_LIBRARY)
      list(APPEND compile_flags "-fPIC")
    endif(${target_type} STREQUAL SHARED_LIBRARY)

    get_directory_property(include_directories INCLUDE_DIRECTORIES)
    foreach(dir ${include_directories})
      list(APPEND compile_flags "-I${dir}")
    endforeach(dir)

    get_directory_property(definitions DEFINITIONS)
    list(APPEND compile_flags ${definitions})

    separate_arguments(compile_flags)

    set(pch_header "${CMAKE_CURRENT_BINARY_DIR}/${target}.hpp")
    file(WRITE "${pch_header}" "#include \"${header}\"\n")

    set(pch_binary "${pch_header}.gch")

    add_custom_command(OUTPUT ${pch_binary} 	
      COMMAND ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1}
      ${compile_flags} -x c++-header -o ${pch_binary} ${pch_header}
      )
    set_property(TARGET ${target} APPEND PROPERTY
      COMPILE_FLAGS "-include ${pch_header} -Winvalid-pch"
      )

    # I wish CMake would support targets depending on files...
    #add_dependencies(${target} ${pch_binary})
    add_custom_target(${target}-pch DEPENDS ${pch_binary})
    add_dependencies(${target} ${target}-pch)
  endif(CMAKE_COMPILER_IS_GNUCXX)
endfunction(ryppl_add_pch_to_target)
