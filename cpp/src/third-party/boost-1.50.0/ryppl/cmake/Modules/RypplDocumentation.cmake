##########################################################################
# Copyright (C) 2008 Douglas Gregor <doug.gregor@gmail.com>              #
# Copyright (C) 2011-2012 Daniel Pfeifer <daniel@pfeifer-mail.de>        #
#                                                                        #
# Distributed under the Boost Software License, Version 1.0.             #
# See accompanying file LICENSE_1_0.txt or copy at                       #
#   http://www.boost.org/LICENSE_1_0.txt                                 #
##########################################################################

if(NOT TARGET documentation)
  add_custom_target(documentation)
endif(NOT TARGET documentation)  

if(CMAKE_HOST_WIN32)
  set(dev_null NUL)
else()
  set(dev_null /dev/null)
endif()

find_package(Boostbook QUIET)
find_package(DBLATEX QUIET)
find_package(FOProcessor QUIET)
find_package(HTMLHelp QUIET)
find_package(XSLTPROC REQUIRED)

get_filename_component(Ryppl_RESOURCE_PATH
  "${CMAKE_CURRENT_LIST_DIR}/../Resources" ABSOLUTE CACHE
  )

#include(CMakeParseArguments)

# Adds documentation for the current library or tool project
#
#   ryppl_documentation(<docbook-file>)
#
function(ryppl_documentation input)
  if(RYPPL_DISABLE_DOCS)
    return()
  endif()

  set(doc_targets)
  set(html_dir "${CMAKE_CURRENT_BINARY_DIR}/html")

  file(COPY
      "${Ryppl_RESOURCE_PATH}/images"
      "${Ryppl_RESOURCE_PATH}/ryppl.css"
    DESTINATION
      "${html_dir}"
    )

  get_filename_component(ext ${input} EXT)
  get_filename_component(name ${input} NAME_WE)
  get_filename_component(input ${input} ABSOLUTE)

  if(HTML_HELP_COMPILER)
    set(hhp_output "${html_dir}/htmlhelp.hhp")
    set(chm_output "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.chm")
    xsltproc(
      INPUT      ${input}
      OUTPUT     ${hhp_output}
      CATALOG    ${BOOSTBOOK_CATALOG}
      STYLESHEET ${Ryppl_RESOURCE_PATH}/docbook-xsl/htmlhelp.xsl
      PARAMETERS "htmlhelp.chm=../${PROJECT_NAME}.chm"
      )
    set(hhc_cmake ${CMAKE_CURRENT_BINARY_DIR}/hhc.cmake)
    file(WRITE ${hhc_cmake}
      "execute_process(COMMAND \"${HTML_HELP_COMPILER}\" htmlhelp.hhp"
      " WORKING_DIRECTORY \"${html_dir}\" OUTPUT_QUIET)"
      )
    add_custom_command(OUTPUT ${chm_output}
      COMMAND "${CMAKE_COMMAND}" -P "${hhc_cmake}"
      DEPENDS ${hhp_output}
      )
    list(APPEND doc_targets ${chm_output})
    install(FILES    "${chm_output}"
      DESTINATION    "doc"
      COMPONENT      "doc"
      CONFIGURATIONS "Release"
      )
  else() # generate HTML and manpages
    set(output_html "${html_dir}/index.html")
    xsltproc(
      INPUT      ${input}
      OUTPUT     ${output_html}
      CATALOG    ${BOOSTBOOK_CATALOG}
      STYLESHEET ${Ryppl_RESOURCE_PATH}/docbook-xsl/xhtml.xsl
      )
    list(APPEND doc_targets ${output_html})
#   set(output_man  ${CMAKE_CURRENT_BINARY_DIR}/man/man.manifest)
#   xsltproc(${output_man} ${BOOSTBOOK_XSL_DIR}/manpages.xsl ${input})
#   list(APPEND doc_targets ${output_man})
    install(DIRECTORY "${html_dir}/"
      DESTINATION     "share/doc/${PROJECT_NAME}"
      COMPONENT       "doc"
      CONFIGURATIONS  "Release"
      )
  endif()

  set(target "${PROJECT_NAME}-doc")
  add_custom_target(${target} DEPENDS ${doc_targets})
  set_target_properties(${target} PROPERTIES
    FOLDER "${PROJECT_NAME}"
    PROJECT_LABEL "${PROJECT_NAME} (documentation)"
    )
  add_dependencies(documentation ${target})

  # build documentation as pdf
  if(DBLATEX_FOUND OR FOPROCESSOR_FOUND)
    set(pdf_dir ${CMAKE_BINARY_DIR}/pdf)
    set(pdf_file ${pdf_dir}/${PROJECT_NAME}.pdf)
    file(MAKE_DIRECTORY ${pdf_dir})

    if(FOPROCESSOR_FOUND)
      set(fop_file ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.fo)
      xsltproc(
        INPUT      ${input}
        OUTPUT     ${fop_file}
        CATALOG    ${BOOSTBOOK_CATALOG}
        STYLESHEET ${BOOSTBOOK_XSL_DIR}/fo.xsl
        PARAMETERS "img.src.path=${CMAKE_CURRENT_BINARY_DIR}/images/"
        )
      add_custom_command(OUTPUT ${pdf_file}
        COMMAND ${FO_PROCESSOR} ${fop_file} ${pdf_file} 2>${dev_null}
        DEPENDS ${fop_file}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
    elseif(DBLATEX_FOUND)
      add_custom_command(OUTPUT ${pdf_file}
        COMMAND ${DBLATEX_EXECUTABLE} -o ${pdf_file} ${input} 2>${dev_null}
        DEPENDS ${input}
        )
    endif()

    set(target "${PROJECT_NAME}-pdf")
    add_custom_target(${target} DEPENDS ${pdf_file})
    set_target_properties(${target} PROPERTIES
      FOLDER "${PROJECT_NAME}"
      PROJECT_LABEL "${PROJECT_NAME} (pdf)"
      )
  endif(DBLATEX_FOUND OR FOPROCESSOR_FOUND)
endfunction(ryppl_documentation)
