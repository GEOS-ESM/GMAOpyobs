#
# Top cmake rules for GMAOpyobs
#

# Boiler plate preamble
# ---------------------
cmake_minimum_required (VERSION 3.13)
cmake_policy (SET CMP0053 NEW)
cmake_policy (SET CMP0054 NEW)

project (
  GMAOpyobs
  VERSION 1.0.0
  LANGUAGES Fortran CXX C)  # Note - CXX is required for ESMF

# Enforce out of source directory builds
# --------------------------------------
  if ("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
      message(SEND_ERROR "In-source builds are disabled. Please
              issue cmake command in separate build directory.")
  endif ("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# Set the default build type to release
# -------------------------------------
if (NOT CMAKE_BUILD_TYPE)
  
     message (STATUS "Setting build type to 'Release' as none was specified.")
     set (CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
     
     # Set the possible values of build type for cmake-gui
     # ---------------------------------------------------
     set_property (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                   "Debug" "Release" "Aggressive") 
  endif ()

# Where to find relevant cmake macros
# -----------------------------------
  if (NOT COMMAND esma)
    foreach (dir cmake @cmake cmake@)
       if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/${dir})
         list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/${dir}")
         set (ESMA_CMAKE_PATH "${CMAKE_CURRENT_LIST_DIR}/${dir}"
              CACHE PATH "Path to ESMA_cmake code")
       endif ()
    endforeach ()
  
    include (esma)

    set (PYOBS_STANDALONE TRUE)
  endif ()
   ecbuild_declare_project()

# Generic DFLAGS
# These should be relocated and/or eliminated.
# --------------------------------------------
  add_definitions(-Dsys${CMAKE_SYSTEM_NAME} -DESMA64)
  add_definitions(${MPI_Fortran_COMPILE_FLAGS})
  include_directories(${MPI_Fortran_INCLUDE_PATH})

# Recursively build source tree
# -----------------------------
  esma_add_subdirectory(env)
  esma_add_subdirectory(src)

# https://www.scivision.dev/cmake-auto-gitignore-build-dir/
# --- auto-ignore build directory
  if(NOT EXISTS ${PROJECT_BINARY_DIR}/.gitignore)
    file(WRITE ${PROJECT_BINARY_DIR}/.gitignore "*")
  endif()

# Piggyback that file into install
# --------------------------------
  install(
     FILES        ${PROJECT_BINARY_DIR}/.gitignore
     DESTINATION  ${CMAKE_INSTALL_PREFIX}
   )

# Adds ability to tar source
# --------------------------
  if (PYOBS_STANDALONE)
    include(esma_cpack OPTIONAL)
  endif ()  
