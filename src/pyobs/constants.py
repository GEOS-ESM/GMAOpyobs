"""
Python version of MAPL Constants.
"""

try:
    from MAPL.constants import *   # use MAPL if available

except:  # hardwire otherwise for portability

    MAPL_PI     = 3.14159265358979323846
    MAPL_DEGREES_TO_RADIANS = MAPL_PI / 180.0
    MAPL_RADIANS_TO_DEGREES = 180.0 / MAPL_PI

    MAPL_UNDEF   = 1.0e15

    MAPL_PSDRY   = 98305.0                # dry surface pressure [Pa]
    MAPL_SECONDS_PER_SIDEREAL_DAY = 86164.0 #s
    MAPL_GRAV    = 9.80665                # m^2/s
    MAPL_RADIUS  = 6371.0E3               # m
    MAPL_OMEGA   = 2.0*MAPL_PI/MAPL_SECONDS_PER_SIDEREAL_DAY    # 1/s
    MAPL_RUNIV   = 8314.47                # J/(Kmole K)
    MAPL_H2OMW   = 18.015                 # kg/Kmole
    MAPL_EARTH_ECCENTRICITY = 8.1819190842622E-2 # --
    MAPL_EARTH_SEMIMAJOR_AXIS = 6378137 # m
    MAPL_KM_PER_DEG = (1.0/(MAPL_RADIUS/1000.)) * MAPL_RADIANS_TO_DEGREES
    MAPL_DEG_PER_KM = (MAPL_RADIUS/1000.) * MAPL_DEGREES_TO_RADIANS

    MAPL_AIRMW   = 28.965                 # kg/Kmole
    MAPL_RDRY    = MAPL_RUNIV/MAPL_AIRMW  # J/(kg K)
    MAPL_CPDRY   = 3.5*MAPL_RDRY          # J/(kg K)
    MAPL_CVDRY   = MAPL_CPDRY-MAPL_RDRY   # J/(kg K)
    MAPL_RVAP    = MAPL_RUNIV/MAPL_H2OMW  # J/(kg K)
    MAPL_CPVAP   = 4.*MAPL_RVAP           # J/(kg K)
    MAPL_CVVAP   = MAPL_CPVAP-MAPL_RVAP   # J/(kg K)
    MAPL_KAPPA   = MAPL_RDRY/MAPL_CPDRY   # (2.0/7.0)
    MAPL_EPSILON = MAPL_H2OMW/MAPL_AIRMW  # --
    MAPL_DELTAP  = MAPL_CPVAP/MAPL_CPDRY  # --
    MAPL_DELTAV  = MAPL_CVVAP/MAPL_CVDRY  # --
    MAPL_GAMMAD  = MAPL_CPDRY/MAPL_CVDRY  # --
    MAPL_RGAS    = MAPL_RDRY              # J/(kg K) (DEPRECATED)
    MAPL_CP      = MAPL_RGAS/MAPL_KAPPA   # J/(kg K) (DEPRECATED)
    MAPL_VIREPS  = 1.0/MAPL_EPSILON-1.0   #          (DEPRECATED)
    MAPL_P00     = 100000.0               # Pa
    MAPL_CAPICE  = 2000.                  # J/(K kg)
    MAPL_CAPWTR  = 4218.                  # J/(K kg)
    MAPL_RHOWTR  = 1000.                  # kg/m^3
    MAPL_NUAIR   = 1.533E-5               # m^2/S (@ 18C)
    MAPL_TICE    = 273.16                 # K
    MAPL_SRFPRS  = 98470                  # Pa
    MAPL_KARMAN  = 0.40                   # --
    MAPL_USMIN   = 1.00                   # m/s
    MAPL_RHO_SEAWATER  = 1026.0          # sea water density [kg/m^3]. SA: should it be = 1026 kg/m^3?
    MAPL_RHO_SEAICE    = 917.0           # sea ice   density [kg/m^3]. SA: should it be = 917  kg/m^3?
    MAPL_RHO_SNOW      = 330.0           # snow density      [kg/m^3]. SA: should it be = 330  kg/m^3?
    MAPL_CELSIUS_TO_KELVIN = 273.15      # K

    MAPL_STFBOL  = 5.6734E-8              # W/(m^2 K^4)
    MAPL_AVOGAD  = 6.023E26               # 1/kmol

    MAPL_O3MW    = 47.9982                # kg/Kmole
    MAPL_ALHL    = 2.4665E6               # J/kg @15C
    MAPL_ALHF    = 3.3370E5               # J/kg
    MAPL_ALHS    = MAPL_ALHL+MAPL_ALHF    # J/kg

