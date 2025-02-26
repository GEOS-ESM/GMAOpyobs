#
# Cmake rules for creating python modules with f2py and 
#

# The OVERRIDE tells esma_set_this() to name the library as
# what we set here rather than the name of the directory
esma_set_this ( OVERRIDE GMAOpyobs )

# cmake requirements
# ------------------
  find_package(F2PY3 REQUIRED)


# Libray
# ------
  esma_add_library (${this}
     SRCS           VegType_Mod.F90
                    VegType_io.c
                    glint_mod.F90  
                    sgp4_mod.F90  
                    TLE_mod.F90
                    csUtils.F90
  )

# Module: binObs
# --------------
  esma_add_f2py3_module(binObs_
     SOURCES        binObs_py.F
     DESTINATION    lib/Python
     INCLUDEDIRS    ${CMAKE_CURRENT_BINARY_DIR}
                    ${CMAKE_BINARY_DIR}/lib
                    ${include_${this}}
   )
   add_dependencies(binObs_ ${this})

 # Module: IGBP
 # ------------
   esma_add_f2py3_module(IGBP_
      SOURCES       IGBP_py.F90
      DESTINATION   lib/Python/pyobs
      ONLY          getsimpleveg getdetailedveg
      LIBRARIES     GMAOpyobs
      INCLUDEDIRS   ${CMAKE_CURRENT_BINARY_DIR}
                    ${CMAKE_BINARY_DIR}/lib
                    ${include_${this}}
    )
    add_dependencies(IGBP_ ${this})

# Module: sgp4
# ------------
  esma_add_f2py3_module(sgp4_
     SOURCES        sgp4_py.F90
     DESTINATION    lib/Python/pyobs
     LIBRARIES      GMAOpyobs
     INCLUDEDIRS    ${CMAKE_CURRENT_BINARY_DIR}
                    ${CMAKE_BINARY_DIR}/lib
                    ${include_${this}}
   )  
   add_dependencies(sgp4_ ${this})

