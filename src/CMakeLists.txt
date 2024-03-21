#
# Cmake rules for pyobs, generic packages for accessing NASA observations.
#

# f2py modules
# ------------
  add_subdirectory(f2py)

# Install 'pyobs' package
# -----------------------
  file(GLOB pyobs_files pyobs/*.py)
  install(FILES ${pyobs_files} DESTINATION lib/Python/pyobs)

# Python executables
# ------------------  
  install(PROGRAMS active_aeronet csBinner wp2traj trj_sampler stn_sampler aop DESTINATION bin)

# TO DO: Convert these to pure python
#        ext_sampler.py  stn_sampler.py stn_xsect.py  

