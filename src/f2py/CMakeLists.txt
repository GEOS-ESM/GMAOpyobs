#
# Cmake rules for creating python modules with f2py and 
#

# cmake requirements
# ------------------

# NOTE: For CI purposes, f2py is hard to support. Until
#       a solution can be found, we add a flag to allow
#       the user to disable f2py in these circumstances
option(USE_F2PY "Use f2py to build python modules" ON)
if (USE_F2PY)

  # This setting is needed on macOS to ensure
  # that the desired Python is found (i.e., if
  # a user has installed their own Python, it
  # should be favored over the system Python or
  # one from Homebrew, etc.)
  set (CMAKE_FIND_FRAMEWORK LAST)

  find_package(
    Python
    COMPONENTS Interpreter Development.Module NumPy
    REQUIRED)

  include(UseF2Py)

  # Library: GMAOpyobs
  # ------------------
  add_library (GMAOpyobs
    VegType_Mod.F90
    VegType_io.c
    glint_mod.F90
    sgp4_mod.F90
    TLE_mod.F90
    csUtils.F90
  )
  install(TARGETS GMAOpyobs DESTINATION lib)

  # Module: binObs
  # --------------

  f2py_object_library(binObs_object OBJECT)
  f2py_generate_module(binObs_ binObs_py.F OUTPUT_VARIABLE binObs_files)
  python_add_library(binObs_ MODULE "${binObs_files}" WITH_SOABI)
  target_link_libraries(binObs_ PRIVATE binObs_object)
  target_link_libraries(binObs_ PRIVATE GMAOpyobs)
  install(TARGETS binObs_ DESTINATION lib/Python)


  # Module: IGBP
  # ------------

  f2py_object_library(IGBP_object OBJECT)
  f2py_generate_module(IGBP_ IGBP_py.F90 OUTPUT_VARIABLE IGBP_files)
  python_add_library(IGBP_ MODULE "${IGBP_files}" WITH_SOABI)
  target_link_libraries(IGBP_ PRIVATE IGBP_object)
  target_link_libraries(IGBP_ PRIVATE GMAOpyobs)
  install(TARGETS IGBP_ DESTINATION lib/Python/pyobs)

    #esma_add_f2py3_module(IGBP_
        #SOURCES       IGBP_py.F90
        #DESTINATION   lib/Python/pyobs
        #ONLY          getsimpleveg getdetailedveg
        #LIBRARIES     GMAOpyobs
        #INCLUDEDIRS   ${CMAKE_CURRENT_BINARY_DIR}
                      #${CMAKE_BINARY_DIR}/lib
                      #${include_GMAOpyobs}
      #)
      #add_dependencies(IGBP_ GMAOpyobs)

  # Module: sgp4
  # ------------

  f2py_object_library(sgp4_object OBJECT)
  f2py_generate_module(sgp4_ sgp4_py.F90 OUTPUT_VARIABLE sgp4_files)
  python_add_library(sgp4_ MODULE "${sgp4_files}" WITH_SOABI)
  target_link_libraries(sgp4_ PRIVATE sgp4_object)
  target_link_libraries(sgp4_ PRIVATE GMAOpyobs)
  install(TARGETS sgp4_ DESTINATION lib/Python/pyobs)


    #esma_add_f2py3_module(sgp4_
      #SOURCES        sgp4_py.F90
      #DESTINATION    lib/Python/pyobs
      #LIBRARIES      GMAOpyobs
      #INCLUDEDIRS    ${CMAKE_CURRENT_BINARY_DIR}
                      #${CMAKE_BINARY_DIR}/lib
                      #${include_GMAOpyobs}
    #)
    #add_dependencies(sgp4_ GMAOpyobs)

endif ()
