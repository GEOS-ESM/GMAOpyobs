#!/bin/env python
"""
   Implements Python interface to several OMPS datasets.
"""

import os
import sys
from types    import *
from glob     import glob

import h5py
from   numpy    import ones, concatenate, savez, load, array, tile
from   datetime import date, datetime, timedelta
from .omps import OMPS_L2, orbits, MISSING

MISSING = -1.267651e+30
kxOMPS = 326
ktAI   = 81


SDS = dict (

      OMPS =  {'ANC_DATA':
                  ('CloudPressure','TerrainPressure'),          
               'GEOLOCATION_DATA':
                  ('GroundPixelQualityFlags',
                   'Latitude','Longitude','RelativeAzimuthAngle','SecondsInDay',
                   'SolarAzimuthAngle','SolarZenithAngle','Time',
                   'UTC_CCSDS_A','ViewingAzimuthAngle','ViewingZenithAngle'),
               'CALIBRATION_DATA':
                  ('Wavelength',),  
               'SCIENCE_DATA':
                  ('UVAerosolIndex',)  
               },


             )

#....................................................................

class OMPS_AI(OMPS_L2):


    """
    Implements Python interface to the OMPS UV aerosol product.
    """


    def __init__ (self,Path,SDS=SDS['OMPS'],keep=None,Verbose=0,only_good=True):


        """
        Creates an OMPS object defining the attributes corresponding
        to the SDS's on input.

        The optional parameter *keep* is used to specify the number of scan
        lines (from the left of the swath) to keep. This is needed for
        coping with the row anomaly problem.

        """

        OMPS_L2.__init__ (self,Path,SDS,keep,Verbose,only_good)

#....................................................................

if __name__ == "__main__":

    if len(sys.argv) <= 1:
      Files = sorted(glob('/nobackup/OMPS/Level2/Y2013/M08/*2013m0819*.h5'))
      q = OMPS_AI(Files,SDS['OMPS'],Verbose=1)
    
#   ----------------------------------------------------------------------

    else:
      """
        To extract synoptic Aerosol index.
      """

      nymd = sys.argv[1]
      hh   = sys.argv[2]
      if len(sys.argv) == 4:
        res = sys.argv[3]
      else:
       res = 'e'

      nymd_ = int(nymd)
      hh_   = int(hh)
      print("nymd: ",nymd)
      print("hour: ",hh)
      year, month, day = (nymd_/10000, (nymd_%10000)/100, nymd_%100)
      syn_time = datetime(year,month,day,hh_,0,0)
      files = orbits("/nobackup/OMPS/Level2","OMPS",syn_time,Verbose=0)
      print(files)

      q = OMPS_AI(files,SDS['OMPS'],Verbose=1)
      print(q.nymd[5],q.nhms[6],q.UVAerosolIndex[0][5])
#   -------------------------------------------------------------------------
