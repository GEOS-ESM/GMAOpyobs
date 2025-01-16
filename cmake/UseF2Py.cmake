# Copyright (c) 2024, Henry Schreiner
# Developed under NSF AWARD OAC-2209877 and by the respective contributors.
# All rights reserved.
#
# SPDX-License-Identifier: Apache-2.0

if(CMAKE_VERSION VERSION_LESS 3.17)
  message(FATAL_ERROR "CMake 3.17+ required")
endif()

include_guard(GLOBAL)

if(TARGET Python::NumPy)
  set(_Python Python CACHE INTERNAL "" FORCE)
elseif(TARGET Python3::NumPy)
  set(_Python Python3 CACHE INTERNAL "" FORCE)
else()
  message(FATAL_ERROR "You must find Python or Python3 with the NumPy component before including F2PY!")
endif()

execute_process(
  COMMAND "${${_Python}_EXECUTABLE}" -c
          "import numpy.f2py; print(numpy.f2py.get_include())"
                  OUTPUT_VARIABLE F2PY_inc_output
                  ERROR_VARIABLE F2PY_inc_error
                  RESULT_VARIABLE F2PY_inc_result
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)

if(NOT F2PY_inc_result EQUAL 0)
  message(FATAL_ERROR "Can't find f2py, got ${F2PY_inc_output} ${F2PY_inc_error}")
endif()

set(F2PY_INCLUDE_DIR "${F2PY_inc_output}" CACHE STRING "" FORCE)
set(F2PY_OBJECT_FILES "${F2PY_inc_output}/fortranobject.c;${F2PY_inc_output}/fortranobject.h" CACHE STRING "" FORCE)
mark_as_advanced(F2PY_INCLUDE_DIR F2PY_OBJECT_FILES)

add_library(F2Py::Headers IMPORTED GLOBAL INTERFACE)
target_include_directories(F2Py::Headers INTERFACE "${F2PY_INCLUDE_DIR}")

function(f2py_object_library NAME TYPE)
  add_library(${NAME} ${TYPE} "${F2PY_INCLUDE_DIR}/fortranobject.c")
  target_link_libraries(${NAME} PUBLIC ${_Python}::NumPy F2Py::Headers)
  if("${TYPE}" STREQUAL "OBJECT")
    set_property(TARGET ${NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
  endif()
endfunction()

function(f2py_generate_module NAME)
  cmake_parse_arguments(
    PARSE_ARGV 1
    F2PY
    "NOLOWER;F77;F90"
    "OUTPUT_DIR;OUTPUT_VARIABLE"
    "F2PY_ARGS"
  )
  set(ALL_FILES ${F2PY_UNPARSED_ARGUMENTS})

  if(NOT ALL_FILES)
    message(FATAL_ERROR "One or more input files must be specified")
  endif()

  if(NOT F2PY_OUTPUT_DIR)
    set(F2PY_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
  endif()

  if(NAME MATCHES "\\.pyf$")
    set(_file_arg "${NAME}")
    get_filename_component(NAME "${NAME}" NAME_WE)
  else()
    set(_file_arg -m ${NAME})
  endif()

  if(F2PY_F77 AND F2PY_F90)
    message(FATAL_ERROR "Can't specify F77 and F90")
  elseif(NOT F2PY_F77 AND NOT F2PY_F90)
    set(HAS_F90_FILE FALSE)

    foreach(file IN LISTS ALL_FILES)
        if("${file}" MATCHES "\\.f90$")
            set(HAS_F90_FILE TRUE)
            break()
        endif()
    endforeach()

    if(HAS_F90_FILE)
      set(F2PY_F90 ON)
    else()
      set(F2PY_F77 ON)
    endif()
  endif()

  if(F2PY_F77)
    set(wrapper_files ${NAME}-f2pywrappers.f)
  else()
    set(wrapper_files ${NAME}-f2pywrappers.f ${NAME}-f2pywrappers2.f90)
  endif()

  if(F2PY_NOLOWER)
    set(lower "--no-lower")
  else()
    set(lower "--lower")
  endif()

  set(abs_all_files)
  foreach(file IN LISTS ALL_FILES)
    if(IS_ABSOLUTE "${file}")
      list(APPEND abs_all_files "${file}")
    else()
      list(APPEND abs_all_files "${CMAKE_CURRENT_SOURCE_DIR}/${file}")
    endif()
  endforeach()

  add_custom_command(
    OUTPUT ${NAME}module.c ${wrapper_files}
    DEPENDS ${ALL_FILES}
    VERBATIM
    COMMAND
      "${${_Python}_EXECUTABLE}" -m numpy.f2py
      "${abs_all_files}" ${_file_arg} ${lower} ${F2PY_F2PY_ARGS}
    COMMAND
      "${CMAKE_COMMAND}" -E touch ${wrapper_files}
    WORKING_DIRECTORY "${F2PY_OUTPUT_DIR}"
    COMMENT
      "F2PY making ${NAME} wrappers"
  )

  if(F2PY_OUTPUT_VARIABLE)
    set(${F2PY_OUTPUT_VARIABLE} ${NAME}module.c ${wrapper_files} PARENT_SCOPE)
  endif()
endfunction()
