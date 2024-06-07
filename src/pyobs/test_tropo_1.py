# -*- coding: utf-8 -*-
"""
Created on Thu May 16 16:50:19 2024

@author: sgasso
"""

#   time:units = "seconds since 2010-01-01 00:00:00" ;
#_FillValue = -1.267651e+30f

import os
import sys

from numpy    import zeros, ones, sqrt, std, mean, unique,\
                     concatenate, where, array, linspace,\
                     shape, arange, interp
from datetime import date, datetime, timedelta
from glob     import glob
import h5py

# from .npz import NPZ

# try:
#     from pyods import ODS # must be inported before pyhdf
# except:
#     pass

# from pyhdf.SD import SD, HDF4Error

from bits import BITS

SDS = dict (
      GEODATA = ('latitude','longitude','viewing_zenith_angle','RelativeAzimuthAngle','solar_zenith_angle',
                  'TerrainPressure','delta_time'),
      SCIDATA = ('FinalAerosolLayerHeight','FinalAerosolOpticalDepth','FinalAerosolSingleScattAlb',
                 'ImaRefractiveIndex','FinalAlgorithmFlags','UVAerosolIndex','FinalAlgorithmFlags','MeasurementQualityFlags'),
      MISC    = ('SurfaceAlbedo', 'Reflectivity','SurfaceType','NormRadiance','AerosolType', 'Irradiance','Residue')
            )

ALIAS = dict ( latitude            ='lat',
               longitude           ='lon',
               viewing_zenith_angle='vza',
               RelativeAzimuthAngle='rza',
               solar_zenith_angle  ='sza',
               TerrainPressure     ='psurf')

# KEYS=[1,2] # SDS to be selected in [1,2]='GEODATA','SCIDATA', note this is 1-order
SELGROUP=['GEODATA','SCIDATA'] # SDS to be selected

CHANNELS=(354.0 ,388.0 ,500.0)
DATE_START = datetime(2010,1,1,0,0,0)
#...........................................................................

class TROPO(object):
    """
    This class implements the MODIS Level 2 AEROSOL products, usually
    referred to as MOD04 (TERRA satellite) and MYD04 (AQUA satellite).
    """

    def __init__ (self,Path,SDS,Verbose=0,only_good=True):
        """
        Creates an AURA object defining the attributes corresponding
        to the SDS's on input.

        The optional parameter *keep* is used to specify the number of scan
        lines (from the left of the swath) to keep. This is needed for
        coping with the row anomaly problem.

        """

        # Initially are lists of numpy arrays for each granule
        # ----------------------------------------------------
        self.verb = Verbose
        self.SDS = SDS

        # Variable names
        # --------------
        self.Names = []
        for group in list(self.SDS.keys()):
            for name in self.SDS[group]:
                self.Names.append(name)
        self.Names += ['nymd','nhms']

        # Create empty lists for SDS to be read from orbit file;
        #  each element of the list contains data for one orbit
        # ------------------------------------------------------
        for name in self.Names:
            self.__dict__[name] = []
        self.time = [] # to hold datetime objects

        # Read each orbit, appending them to the list
        # -------------------------------------------
        if type(Path) is list:
            if len(Path) == 0:
                self.nobs = 0
                print("WARNING: Empty AURA object created")
                return
        else:
            Path = [Path, ]

        self._readList(Path) ##### note this read the actual file, see example in mxd04.py

    	# Make each attribute a single numpy array
    # ----------------------------------------
        # for name in self.SDS:
        #      try:
        #          self.__dict__[name] = concatenate(self.__dict__[name])
        #      except:
        #          print("Failed concatenating "+name)

        # for group in list(self.SDS.keys()):
        #     for name in self.SDS[group]:
        #         try :
        #              self.__dict__[name] = concatenate(self.__dict__[name])
        #         except:
        #             print("Failed concatenating,name=  "+name)
        #             print("Failed concatenating,group=  "+group)

        for igroup in SELGROUP:
            key_content_names=tuple(self.SDS.get(igroup))
            for sds_name in key_content_names:
                try :
                     self.__dict__[sds_name] = concatenate(self.__dict__[sds_name])
                except:
                    print("Failed concatenating,name=  "+sds_name)
                    print("Failed concatenating,group=  "+igroup)


       # Create corresponding python time, from mod04.py
       # --------------------------------
        # self.Time = array([DATE_START+timedelta(seconds=s) for s in self.Scan_Start_Time])

    # ODS friendly attributes
    # -----------------------
    pass

    ##-------------------- from Aura.py
    def _readList(self,List):
        """
        Recursively, look for files in list; list items can
        be files or directories.
        """
        for item in List:
            # print('hello 0= ',item)
            if os.path.isdir(item):      self._readDir(item)
            elif os.path.isfile(item):   self._readOrbit(item)
            else:
                print("%s is not a valid file or directory, ignoring it"%item)

    #-------------- from Aura.py, called by  _readList(self,List)
    def _readDir(self,dir):
        """Recursively, look for files in directory."""
        for item in os.listdir(dir):
            path = dir + os.sep + item
        if os.path.isdir(path):      self._readDir(path)
        elif os.path.isfile(path):   self._readOrbit(path)
        else:
            print("%s is not a valid file or directory, ignoring it"%item)

    #-------------- from Aura.py called _readList(self,List)
    def _readOrbit(self,filename):
        """Reads one OMRUV/OMSO2 granule with Level 2 aerosol data."""

    # Reference time
    # --------------
        # REF_DATE = = datetime(2010,1,1,0,0,0)

    # Open the AURA file and loop over the datasets,
    # extracting GEOLOCATION and Data fields
    # ----------------------------------------------
        if self.verb:
            print("... working with file <%s>"%filename)
            f = h5py.File(filename,mode='r')
        print(' ')

        #### added make sure  to select the group name (which are ht key names in the dictionary)
        #### then once the group is selected get the respective elelent sin those keys and get the data
        # keys_to_select=[1,2] # select GEODATA amd SCIDATA, using index  start in 1
        # keys_to_select=KEYS
        # for igroup in keys_to_select:
        #     key_name = list(self.SDS.keys())[igroup - 1]
        #     key_content_names=tuple(self.SDS.get(key_name,"Group name not found"))

        for igroup in SELGROUP:
            key_content_names=tuple(self.SDS.get(igroup))
            # print(f"The name of the key at position {igroup} is: {key_name}")
            print(f"The name of the selected group is {igroup} and it contains the following names: {key_content_names}")
            print('')
            ### now open corresponding group
            # g=f.get(key_name)
            g=f.get(igroup)
            # print(key_content_names[0])
            # print('')

            ### now loop over the SDS names stores in key_content_name
        # Read select variables (reshape to allow concatenation later)
        # ------------------------------------------------------------
            for ig in key_content_names:
                # print('hello2 :',ig)
                if ig == 'delta_time':
                   delta_time=g.get(ig)
                   delta_time=delta_time.astype('float')
                   Time      = g.get('time')
                   # print(dir(delta_time))
                   # print(len(delta_time))
                   # print(dir(delta_time))
                   nobs = len(delta_time)
                   nymd  = ones(nobs).astype('int')
                   nhms  = ones(nobs).astype('int')
                   # print(nobs,nymd,nhms)
                   print('')
                   for i in range(nobs):
                     t_secs = delta_time[i]
                     n_secs = timedelta(seconds=t_secs)
                     t = DATE_START + n_secs

                     yy, mm, dd = (t.year,t.month,t.day)
                     h, m, s = (t.hour,t.minute,t.second)

                     nymd[i] = 10000 * yy + 100 * mm + dd
                     nhms[i] = 10000 * h  + 100 * m  + s
                     self.time.append(t)

                   self.delta_time.append(delta_time[:]+Time) # time as on file
                   self.nymd.append(nymd)
                   self.nhms.append(nhms)

                else:
                     data = g.get(ig)
                     self.__dict__[ig].append(data)
                     print('Read .... ',self.__dict__[ig])
            # print('End of _readOrbit \n')
        print('self=',dir(self))

###### -------- test area
#............................................................................

if __name__ == "__main__":

    current_os=sys.platform
    if current_os=='win32':
        pthin="D:\Satellite\Tropomi\Level2/"
    else:
        pthin='/Users/sgasso/'
#         pthin='/nobackup/TROPOMAER/2019/255/'
        pthin='/nobackup/TROPOMAER/2019/255/'

#     filename='TROPOMI-Sentinel-5P_L2-TROPOMAER_2020m0616t140632-o13864_v01-2021m1113t104402.nc'
    filename='TROPOMI-Sentinel-5P_L2-TROPOMAER_2019m0912t105038-o09918_v01-2021m1016t191711.nc'
#    syn_time = datetime(2008,6,30,0,0,0)
    syn_time = datetime(2019,9,12,10,50,0)
    Files = pthin+filename
    # Files = granules('/discover/nobackup/dao_ops/intermediate/flk/modis','MOD04',syn_time,coll='006')
    # select group of data from dictionary defined begining of code
    keys=[1,2] # select GEODATA amd SCIDATA, using index  start in 1


    data = TROPO(Files,SDS,Verbose=1,only_good=True)