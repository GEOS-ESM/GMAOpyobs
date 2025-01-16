#
# Top cmake rules for GMAOpyobs
#

# Boiler plate preamble
# ---------------------
cmake_minimum_required (VERSION 3.24)
cmake_policy (SET CMP0053 NEW)
cmake_policy (SET CMP0054 NEW)
cmake_policy (SET CMP0132 NEW)

project (
  GMAOpyobs
  VERSION 1.1.0
  LANGUAGES Fortran C
)

# Enforce out of source directory builds
# --------------------------------------
if ("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are disabled. Please
          issue cmake command in separate build directory.")
endif ("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Recursively build source tree
# -----------------------------

# mepo can now clone subrepos in three styles
set (ESMA_ENV_DIRS
  env
  @env
  env@
)
foreach (dir IN LISTS ESMA_ENV_DIRS)
  if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/${dir})
    add_subdirectory (${dir})
  endif ()
endforeach ()

add_subdirectory(src)

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
