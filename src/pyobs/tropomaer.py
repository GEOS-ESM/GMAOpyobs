# -*- coding: utf-8 -*-
"""
July/15/2024 - _v2 : built from _v1, this adds a plotting of TROPO data,
good for plotting along CALIOP data

Jul/08/2024 Former code test_tropo_caliop_3.py now renamed to Tropo_Caliop_1.py

Jul/05/2024 _v3 : based on v2 start to implement suggested cahnges in v2

This code works

Explore stretuctures
TROP=TROPO()
vars(TROPO)
TROP.__dict__.items()
TROP.__dict__.keys() or print(CALI.keys())
help(TROP)

Jul/04/2024 - This code collocates TROPOMAER LEve2 data with
CALIPSO

Only spatial collocation (do temporal check in a later version)
assume the two closest Caliop and Tropomi orbit have been found.

This code
1) read user defined input TROPO file. Calipso data has been preproceesed and
saved in output.nc
2) indexes for CALIPSO pixels that are contained in TROPO pixels are found
3) Also, an array with  with the indexes of CALIPSO profile that are in each TROPO
are found

Because size of array some maniuplation must be done.
As a result code is very fast (~ 1 sec per orbit)


Things left to do:
    Save TROPO data (UVA, Z, AOD) along the CALIOP profile
    Do a plot verification to check if things make sense (see http://meteothink.org/examples/meteoinfolab/satellite/calipso.html)
    Start to optimize code:
        Make sure of indexes are mapped to original file.
        (see if I do need to save the intermediate CALIOP file)
        Think about the temporal matching
        Generalize code over all orbits in day and then over all days


@author: sgassoumd
"""


import os, sys, time
import numpy as np
from datetime import date, datetime, timedelta

import h5py
import netCDF4 as nc

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm  # For custom colormaps
import cartopy.crs as ccrs  # For cartographic projections
import cartopy.feature as cfeature  # For adding features to the map
# from glob     import glob



# SDS = dict (
#       GEODATA = ('latitude','longitude','viewing_zenith_angle','RelativeAzimuthAngle','solar_zenith_angle',
#                   'TerrainPressure','delta_time'),
#       SCIDATA = ('FinalAerosolLayerHeight','FinalAerosolOpticalDepth','FinalAerosolSingleScattAlb',
#                  'ImaRefractiveIndex','FinalAlgorithmFlags','UVAerosolIndex','FinalAlgorithmFlags','MeasurementQualityFlags'),
#       MISC    = ('SurfaceAlbedo', 'Reflectivity','SurfaceType','NormRadiance','AerosolType', 'Irradiance','Residue')
#             )

# ALIAS = dict ( latitude            ='lat',
#                longitude           ='lon',
#                viewing_zenith_angle='vza',
#                RelativeAzimuthAngle='rza',
#                solar_zenith_angle  ='sza',
#                TerrainPressure     ='psurf')


SDS = dict (
      GEODATA = ('latitude','longitude','latitude_bounds','longitude_bounds','delta_time'),
      SCIDATA = ('UVAerosolIndex','FinalAlgorithmFlags','FinalAerosolLayerHeight',\
                 'SurfaceType','Reflectivity')
            )

ALIAS = dict ( latitude            ='lat',
               longitude           ='lon',
               latitude_bounds     ='lat4',
               longitude_bounds    ='lon4',
               delta_time          ='dtime', # seconds since DATE_START = datetime(2010,1,1,0,0,0)
               UVAerosolIndex      ='uvai',
               FinalAerosolLayerHeight='ztro',
               FinalAlgorithmFlags ='algQA',
               SurfaceType         ='surftype',
               Reflectivity        ='ref')

# KEYS=[1,2] # SDS to be selected in [1,2]='GEODATA','SCIDATA', note this is 1-order
SELGROUP=['GEODATA','SCIDATA'] # SDS to be selected

CHANNELS=(354.0 ,388.0 ,500.0)
DATE_START = datetime(2010,1,1,0,0,0)
#...........................................................................

class TROPOAER_L2(object):
    """
    This class loads TROPOMI data, for variable input types: single files, several directories.
    Based on existing codes in GEOSPYobs aura.py and mxd04.py
    """

    def __init__ (self,Path,SDS=SDS,Verbose=0,only_good=True,alias=None):
        """
        Put General Ancilliary info here

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

       # Add/Substitute some aliases if given
       # ------------------------------------
        self.ALIAS = ALIAS.copy()
        if alias is not None:
           for a in alias: self.ALIAS[a] = alias[a]
        # print('self.ALIAS ', self.ALIAS)

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
                print("WARNING: Empty TROPO object created")
                return
        else:
            Path = [Path, ]
        # print('List of Files to read:')
        # print(Path)
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
                     self.__dict__[sds_name] = np.concatenate(self.__dict__[sds_name])
                except:
                    print("Failed concatenating,name=  "+sds_name)
                    print("Failed concatenating,group=  "+igroup)

       # USE ALIAS dic to shorten the variable names
        Alias = list(self.ALIAS.keys())
        for group in list(self.SDS.keys()):
            for name in self.SDS[group]:
                if name in Alias:
                   self.__dict__[self.ALIAS[name]] = self.__dict__[name]
                   ###  now remove the longer name  attribute
                   delattr(self, name)
                else:                   # this maybe removed if i want to keep variables
                    delattr(self, name) # that are no ALIAS for

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
                sys.exit('Exit _readList')

    #-------------- from Aura.py, called by  _readList(self,List)
    def _readDir(self,dir):
        """Recursively, look for files in directory."""
        for item in os.listdir(dir):
            path = dir + os.sep + item
        if os.path.isdir(path):      self._readDir(path)
        elif os.path.isfile(path):   self._readOrbit(path)
        else:
            print("%s is not a valid file or directory, ignoring it"%item)
            sys.exit('Exit _readDir')
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
            print("\n ... working with file <%s>"%filename)
        f = h5py.File(filename,mode='r')
        print(' ')

        #### added make sure  to select the group name (which are ht key names in the dictionary)
        #### then once the group is selected get the respective elelent sin those keys and get the data
        #### check the class type stored in the key. It has to be a tuple or a list. Also if it is a single
        #### item (astring with a name) , there is no need to turn it into tuple

        for igroup in SELGROUP:
            key_content_names=self.SDS[igroup]

            print(f"The name of the selected group is {igroup} and it contains the following names: {key_content_names}")

            ### now open corresponding group
            g=f.get(igroup)

            ### now loop over the SDS names stores in key_content_name
        # Read select variables (reshape to allow concatenation later)
        # ------------------------------------------------------------
            for ig in key_content_names:
                # print('hello2 :',ig)
                if ig == 'delta_time': #  "milliseconds since 2019-08-13 00:00:00"
                   delta_time=g.get(ig)
                   delta_time=delta_time.astype('float')
                   Time      = g.get('time')
                   nobs = len(delta_time)
                   nymd  = np.ones(nobs).astype('int')
                   nhms  = np.ones(nobs).astype('int')
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
                # print('')
            # print('End of _readOrbit \n')
        print('   ')


##### Now the modules
#............................................................................
###### -------- test area
if __name__ == "__main__":

   current_os=sys.platform
   if current_os=='win32': #in PC office
      pthin="D:\Satellite\Tropomi\Level2\/226\/"
   elif current_os=='darwin': #in Laptop
      # pthin='/Users/sgasso/'
      pthin='/Volumes/ExtData1/SatData/Tropomi/Level2/'
   elif current_os == 'linux': # in Calculon
      pthin='/nobackup/TROPOMAER/2019/226/'
   else:
       print('Current operating system no recognized.')
       print('Cannot set path to Level1 files. Terminate Run')
       sys.exit()


   filename_tro='TROPOMI-Sentinel-5P_L2-TROPOMAER_2019m0814t131912-o09508_v01-2021m1015t022508.nc'
   syn_time = datetime(2019,8,14,13,19,12)
   Files = pthin+filename_tro
   # print(Files)
   # Files = granules('/discover/nobackup/dao_ops/intermediate/flk/modis','MOD04',syn_time,coll='006')
   # select group of data from dictionary defined begining of code
   keys=[1,2,3] # select GEODATA amd SCIDATA, using index  start in 1


   TROP = TROPOAER_L2(Files,Verbose=1,only_good=True)
