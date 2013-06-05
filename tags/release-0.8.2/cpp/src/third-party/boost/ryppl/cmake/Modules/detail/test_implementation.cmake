################################################################################
# Copyright (C) 2012 Daniel Pfeifer <daniel@pfeifer-mail.de>                   #
#                                                                              #
# Distributed under the Boost Software License, Version 1.0.                   #
# See accompanying file LICENSE_1_0.txt or copy at                             #
#   http://www.boost.org/LICENSE_1_0.txt                                       #
################################################################################


set(__boost_test_run "${CMAKE_CURRENT_LIST_DIR}/test_run.cmake")
set(__boost_test_summary "${CMAKE_CURRENT_LIST_DIR}/test_summary.cmake")


macro(__boost_add_test_compile fail)
  get_filename_component(name ${FILE} NAME_WE)
  set(output ${name}_ok.txt)

  get_filename_component(SOURCE ${FILE} ABSOLUTE)
  set(OBJECT ${name}.o)

  string(TOUPPER "${CMAKE_BUILD_TYPE}" build_type)
  set(FLAGS ${CMAKE_CXX_FLAGS_${build_type}})

  get_directory_property(include_directories INCLUDE_DIRECTORIES)
  foreach(dir ${include_directories})
    list(APPEND FLAGS "-I${dir}")
  endforeach(dir)

  string(REGEX REPLACE "<([A-Z_]+)>" "@\\1@" compile
    "${CMAKE_CXX_COMPILE_OBJECT}"
    )
  string(CONFIGURE "${compile}" compile @ONLY)
  string(REPLACE ";" " " compile "${compile}")

  add_custom_command(OUTPUT ${output}
    COMMAND ${CMAKE_COMMAND}
      -D "OUTPUT=${output}"
      -D "COMMANDS=COMPILE"
      -D "COMPILE=${compile}"
      -D "COMPILE_FAIL=${fail}"
      -P "${__boost_test_run}"
    DEPENDS ${FILE}
    COMMENT "compile test: ${name}"
    )

  list(APPEND TEST_NAMES ${name})
  list(APPEND TEST_FILES ${output})
endmacro(__boost_add_test_compile)


macro(__boost_add_test_link link_rule fail)
  get_filename_component(name ${FILE} NAME_WE)
  set(output ${name}_ok.txt)

  get_filename_component(SOURCE ${FILE} ABSOLUTE)
  set(TARGET ${name}_ok)
  set(OBJECT ${name}.o)
  set(OBJECTS ${OBJECT})

  string(REGEX REPLACE "<([A-Z_]+)>" "@\\1@" compile
    "${CMAKE_CXX_COMPILE_OBJECT}"
    )
  string(REGEX REPLACE "<([A-Z_]+)>" "@\\1@" link
    "${link_rule}"
    )
  string(CONFIGURE "${compile}" compile @ONLY)
  string(CONFIGURE "${link}" link @ONLY)

  add_custom_command(OUTPUT ${output}
    COMMAND ${CMAKE_COMMAND}
      -D "OUTPUT=${output}"
      -D "COMMANDS=COMPILE LINK"
      -D "COMPILE=${compile}"
      -D "LINK=${link}"
      -D "LINK_FAIL=${fail}"
      -P "${__boost_test_run}"
    DEPENDS ${FILE}
    COMMENT "link test: ${name}"
    )

  list(APPEND TEST_NAMES ${name})
  list(APPEND TEST_FILES ${output})
endmacro(__boost_add_test_link)


macro(__boost_add_test_run driver fail)
  get_filename_component(name ${FILE} NAME_WE)
  set(output ${name}_ok.txt)

  add_custom_command(OUTPUT ${output}
    COMMAND ${CMAKE_COMMAND}
      -D "OUTPUT=${output}"
      -D "COMMANDS=RUN"
      -D "RUN=$<TARGET_FILE:${driver}> ${name}"
      -D "RUN_FAIL=${fail}"
      -P "${__boost_test_run}"
    DEPENDS ${FILE}
    COMMENT "Running test: ${name}"
    )

  list(APPEND TEST_NAMES ${name})
  list(APPEND TEST_FILES ${output})
endmacro(__boost_add_test_run)


macro(__boost_add_test_run_deprecated driver fail)
  get_filename_component(name ${FILE} NAME_WE)
  set(name ${target}-${name})
  set(output ${name}_ok.txt)

  set(testdriver ${driver}${SUFFIX})
  math(EXPR SUFFIX "${SUFFIX} + 1")
  add_executable(${testdriver} EXCLUDE_FROM_ALL
    ${FILE}
    ${TEST_ADDITIONAL_SOURCES}
    )
  target_link_libraries(${testdriver}
    ${TEST_LINK_LIBRARIES}
    )

  add_custom_command(OUTPUT ${output}
    COMMAND ${CMAKE_COMMAND}
      -D "OUTPUT=${output}"
      -D "COMMANDS=RUN"
      -D "RUN=$<TARGET_FILE:${testdriver}>"
      -D "RUN_FAIL=${fail}"
      -P "${__boost_test_run}"
    DEPENDS ${FILE}
    COMMENT "Running test: ${name}"
    )

  list(APPEND TEST_NAMES ${name})
  list(APPEND TEST_FILES ${output})
endmacro(__boost_add_test_run_deprecated)


macro(__boost_add_test_python fail)
  get_filename_component(name ${FILE} NAME_WE)
  set(output ${name}_ok.txt)

  set(module "${PROJECT_NAME}-test-${name}-ext")
  add_library(${module} MODULE EXCLUDE_FROM_ALL ${FILE})
  target_link_libraries(${module} ${TEST_LINK_LIBRARIES})
  set_target_properties(${module} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    OUTPUT_NAME "${name}_ext"
    PREFIX ""
    )

  add_custom_command(OUTPUT ${output}
    COMMAND ${CMAKE_COMMAND}
      -D "OUTPUT=${output}"
      -D "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
      -D "COMMANDS=RUN"
      -D "RUN=${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.py"
      -D "RUN_FAIL=${fail}"
      -P "${__boost_test_run}"
    DEPENDS ${FILE} ${module}
    COMMENT "Running test: ${name}"
    )

  list(APPEND TEST_NAMES ${name})
  list(APPEND TEST_FILES ${output})
endmacro(__boost_add_test_python)
