# -*- coding: utf-8 -*-
"""

Generic script to load and plot aerosol data from TROPOAER_L2 product. 

Sept/2024 - by Santiago Gassó , sgasso@umd.edu

Internal notes: 

@author: sgassoumd
"""


import os, sys, time
import numpy as np
from datetime import date, datetime, timedelta

import h5py
import netCDF4 as nc

### plotting modules 
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm  # For custom colormaps
import cartopy.crs as ccrs  # For cartographic projections
import cartopy.feature as cfeature  # For adding features to the map
# from glob     import glob
# import pdb


#default SDSs if not specifically input by user
SDS = dict (
      GEODATA = ('latitude','longitude','latitude_bounds','longitude_bounds','delta_time'),
      SCIDATA = ('UVAerosolIndex','FinalAlgorithmFlags','FinalAerosolLayerHeight',\
                 'SurfaceType','Reflectivity')
            )

ALIAS = dict ( latitude            ='lat',
               longitude           ='lon',
               latitude_bounds     ='lat4',
               longitude_bounds    ='lon4',
               delta_time          ='dtime', # seconds since start day of this orbit 
#               time                ='daystart', # seconds since 2010-01-01
               UVAerosolIndex      ='uvai',
               FinalAerosolLayerHeight='ztro',
               FinalAlgorithmFlags ='algQA',
               SurfaceType         ='surftype',
               Reflectivity        ='ref')

# KEYS=[1,2] # SDS to be selected in [1,2]='GEODATA','SCIDATA', note this is 1-order
# SELGROUP=['GEODATA','SCIDATA'] # SDS to be selected
# SELGROUP=['GEODATA'] # SDS to be selected

CHANNELS=(354.0 ,388.0 ,500.0)
DATE_START = datetime(2010,1,1,0,0,0)
#...........................................................................

class TROPOAER_L2(object):
    """
    This class loads TROPOMI data, for variable input types: single files, several directories.
    Based on existing codes in GEOSPYobs aura.py and mxd04.py
    """

    def __init__ (self,Path,SDS=SDS,Verbose=0,GEO=False,alias=None):
        """
        Put General Ancilliary info here

        """

        # Initially are lists of numpy arrays for each granule
        # ----------------------------------------------------
        self.verb = Verbose
        self.SDS = SDS
        # if input GEO is TRUE then only get Geolocation group
        if GEO:
            self.SELGROUP=['GEODATA']
            if self.verb: print('User selected GEODATA only')
        else:
            self.SELGROUP=['GEODATA','SCIDATA']


        # Variable names
        # --------------
        self.Names = []
        for group in list(self.SDS.keys()):
            for name in self.SDS[group]:
                self.Names.append(name)
        self.Names += ['nymd','nhms'] # inherited from old code. Might be useful. 
        
       # Add/Substitute some aliases if given
       # ------------------------------------
        self.ALIAS = ALIAS.copy()
        if alias is not None:
           for a in alias: self.ALIAS[a] = alias[a]
        # print('self.ALIAS 2 ', self.ALIAS)
        
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
        self._readList(Path) ##### note this read the actual file, see example in mxd04.py


    # Make each attribute a single numpy array
    # ----------------------------------------
        for igroup in self.SELGROUP:
            key_content_names=tuple(self.SDS.get(igroup))
            for sds_name in key_content_names:
                try :
                     self.__dict__[sds_name] = np.concatenate(self.__dict__[sds_name])
                except:
                    print("tropomi_le_reader: Failed concatenating,name=  "+sds_name)
                    print("tropomi_le_reader: Failed concatenating,group=  "+igroup)

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
       
       # remove temporary attribute that does not need to be output  
        delattr(self,'SELGROUP')

    def _readList(self,List):
        """
        Recursively, look for files in list; list items can
        be files or directories.
        """
        for item in List:
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
    # Open the hdf file and loop over the datasets,
    # extracting GEOLOCATION and Data fields
    # ----------------------------------------------
        if self.verb: print("\n ... working with file <%s>"%filename)
        f = h5py.File(filename,mode='r')
        for igroup in self.SELGROUP:
            key_content_names=self.SDS[igroup]
            if self.verb: 
                print(f"The name of the selected group is {igroup} and it contains the following names: {key_content_names}")
            ### now open corresponding group
            g=f.get(igroup)
            ### now loop over the SDS names stores in key_content_name
        # Read select variables (reshape to allow concatenation later)
        # ------------------------------------------------------------
            for ig in key_content_names:
<<<<<<< HEAD
=======
                #print('hello2 :',ig,igroup)
>>>>>>> feature/sgassoumd/mpl_reader
                if ig == 'delta_time': #  "MILLIseconds  since beginning of day for orbit"
                   delta_time=g.get(ig)
                   delta_time=delta_time.astype('float')
                   self.__dict__[ig].append(delta_time[:]) # add reference time
                   Time      = g.get('time')
                   time_value = Time[()]  # or Time.value , = int32
                   self.time.append(time_value)
<<<<<<< HEAD
#                   # formatted_date = (datetime(2010, 1, 1) + timedelta(seconds=int(time_value))).strftime('%Y-%m-%dT%H:%M:%S')
                   # breakpoint()

=======

                   # breakpoint()
>>>>>>> feature/sgassoumd/mpl_reader
#                   nobs = len(delta_time)
#                   nymd  = np.ones(nobs).astype('int')
#                   nhms  = np.ones(nobs).astype('int')
                   # print('')
#                   for i in range(nobs):
#                     t_secs = delta_time[i]
#                     n_secs = timedelta(seconds=t_secs)
#                     t = DATE_START + n_secs
#                     yy, mm, dd = (t.year,t.month,t.day)
#                     h, m, s = (t.hour,t.minute,t.second)
#                     nymd[i] = 10000 * yy + 100 * mm + dd
#                     nhms[i] = 10000 * h  + 100 * m  + s
#                     self.time.append(t)

#                   self.delta_time.append(delta_time[:]+Time) # time as on file
                   # self.__dict__[ig].append(delta_time[:]+time_value) # add reference time
#                   self.nymd.append(nymd)
#                   self.nhms.append(nhms)
                   # breakpoint()
#                elif ig == 'time':
                   #breakpoint()

<<<<<<< HEAD
=======
#                   # formatted_date = (datetime(2010, 1, 1) + timedelta(seconds=int(time_value))).strftime('%Y-%m-%dT%H:%M:%S')
>>>>>>> feature/sgassoumd/mpl_reader
#                   # if self.verb: print('Time = ', time_value)
#                   if self.verb: print('Reference Day : ', (DATE_START + timedelta(seconds=int(time_value))).strftime('%Y-%m-%dT%H:%M:%S'))
#                   self.__dict__[ig].append(time_value)

                else:
                     data = g.get(ig)
                     self.__dict__[ig].append(data)
                # if self.verb: print('Read .... ',self.__dict__[ig])
            # print('End of _readOrbit \n')

<<<<<<< HEAD
=======


>>>>>>> feature/sgassoumd/mpl_reader
#............................................................................
###### -------- test area
if __name__ == "__main__":

    current_working_directory = os.getcwd()
    print('\nThis code is running from directory :', current_working_directory)
    # Get user inputs
    now_os   = platform.system().lower()
    now_computer = platform.node()
    print(f'This code is executed in computer {now_computer} \n')
    # Your existing get_path function remains the same
    if now_os=='win32':
        base_path    = 'C:/Users/sgasso/Downloads/'
        pth_fig_out='C:/Users/sgasso/OneDrive - NASA/ToShare/2025/GEOS/Pyfigures/'
    elif now_os == 'darwin':
        base_path= '/Volumes/ExtData1/SatData/Tropomi/Level2/2023/359/'
        pth_fig_out='/Users/sgasso/Library/CloudStorage/OneDrive-NASA/ToShare/2025/AOCH/PyFigures/'
    elif now_os == 'linux' and "calculon" in now_computer:
        base_path = '/nobackup/TROPOMAER/2023/359/'
        pth_fig_out = ''
    elif now_os == 'linux' and "discover" in now_computer:
        base_path = '/nobackup/CALIPSO/Level1.5/Santiago/'
    else:
        print('Current operating system no recognized.')
        print('Cannot set path to MPL  files. Terminate Run')
        sys.exit()

 
    filename_tro='TROPOMI-Sentinel-5P_L2-TROPOMAER_2023m1225t132102-o32123_v01-2023m1231t105113.nc'
    Files = base_path+filename_tro

    # Extract only GEO data
    TROP = TROPOAER_L2(Files,GEO=True,Verbose=0)
    ## extract defaule list of variables , see list in tropomaer.py
    # TROP2 = TROPOAER_L2(Files,Verbose=1)


####  Miscellanous  Check date
    # Both conversions in one-liners
    date1    = (datetime(2010, 1, 1) + timedelta(seconds=int(TROP.time[0])))
    date1_str= date1.strftime('%Y-%m-%dT%H:%M:%S')
    date2     = date1 + timedelta(milliseconds=int(TROP.dtime[0]))
    date2_str = date2.strftime('%Y-%m-%dT%H:%M:%S')
    # Print results
    print(f"Reference day: {date1}")
    print(f"Time first scanline : {date2}")