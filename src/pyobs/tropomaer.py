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

class TROPO(object):
    """
    This class loads TROPOMI data, for variable input types: single files, several directories.
    Based on existing codes in GEOSPYobs aura.py and mxd04.py
    """

    def __init__ (self,Path,SDS,Verbose=0,only_good=True,alias=None):
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

       ### clean up output
        delattr(self,'SDS'); delattr(self,'Names')
        delattr(self,'ALIAS');delattr(self,'verb')


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
            # if isinstance(self.SDS.get(igroup), tuple):
            #     key_content_names=self.SDS.get(igroup)
            # elif isinstance(self.SDS.get(igroup), str):
            #     key_content_names=self.SDS.get(igroup)
            # else :
            key_content_names=tuple(self.SDS.get(igroup))

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

def load_netcdf_arrays(file_path,array_names):
    """
    Loads specified arrays from a NetCDF file into memory using netCDF4.
    Args:
        file_path (str): Path to the NetCDF file.
        array_names (list of str) : names of variables to read
    Returns:
        dict: A dictionary containing the loaded arrays.
    """
    try:
        # Open the NetCDF file
        dataset = nc.Dataset(file_path, "r")
        # Initialize an empty dictionary to store the arrays
        loaded_arrays = {}
        # Load each specified array
        for array_name in array_names:
            if array_name in dataset.variables:
                loaded_arrays[array_name] = dataset.variables[array_name][:]  # Load data
            else:
                print(f"Array '{array_name}' not found in the NetCDF file.")
        # Close the dataset
        dataset.close()
        return loaded_arrays
    except Exception as e:
        print(f"Error loading NetCDF arrays: {e}")
        return None

def haversine(lat1, lon1, lat2, lon2):
        """
        Calculate the great circle distance between two points
        on the earth (specified in decimal degrees)

        lat1, lon1, lat2, lon2 are arrays
        """
        # convert decimal degrees to radians
        lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

        # haversine formula
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
        c = 2 * np.arcsin(np.sqrt(a))
        r = 6371 # Radius of earth in kilometers. Use 3956 for miles
        return c * r


def find_first_last_zero(A):
    left, right = 0, len(A) - 1
    first = -1
    last = -1

    # Find the first occurrence of 0
    while left <= right:
        mid = (left + right) // 2
        if A[mid] == 0:
            first = mid
            right = mid - 1
        elif A[mid] < 0:
            left = mid + 1
        else:
            right = mid - 1

    # Find the last occurrence of 0
    left, right = 0, len(A) - 1
    while left <= right:
        mid = (left + right) // 2
        if A[mid] == 0:
            last = mid
            left = mid + 1
        elif A[mid] < 0:
            left = mid + 1
        else:
            right = mid - 1

    return first, last


def f_find_cal_in_tro_pix(scan0,saveit_cal0):
    ######## Now find index of first and last CALIOP
    ######## profile inside each TROPO pixel
    # Start timer
    start_time = time.perf_counter()
    d0=100
    Ntrop=scan0[1]-scan0[0]+1
    scan_lines = np.arange(scan0[0],scan0[1]+1)
    saveit_trop0=-999*np.ones((Ntrop, 2), dtype=int)
    last_line= Ntrop - d0
    for itro,line in enumerate(scan_lines): #itrop=0,1,2,...., line=scan[0],scan[0]+1,scan[0]+2,....
        ### take advantage indexes are continous and appear once as group
        ### preselect a group in Caliop array
        if itro == 0:
            # tmp0= line - d0 # get this value
            tmp1=saveit_cal0[0:d0,0]
            shift = 0
        elif line < last_line :
            ### get the index of last element found in previous search
            ### and get a block of data starting from ther
            shift=saveit_trop0[itro-1,1]+1
            tmp1= saveit_cal0[shift:shift+d0,0]
            # tmp1= cal_ind_fine[itro:itro+d0]
        else :
            shift=saveit_trop0[itro-1,1]+1
            tmp1= saveit_cal0[shift:-1,0]
            # tmp1= cal_ind_fine[start:-1]
        # tmp1 is an array of ~d0 lenght with indexes
        tmp2=tmp1 - line # where tmp2=0 save those indexes.
        first,last = find_first_last_zero(tmp2)
        saveit_trop0[itro,:] = cal_ind_fine[first+shift],cal_ind_fine[last+shift]
        if saveit_trop0[itro,0] < 0 : sys.exit('No found. Stop. Line 547')

    # End timer
    end_time = time.perf_counter()
    # Calculate elapsed time
    elapsed_time = end_time - start_time
    print("  Secs. Elapsed time for matching CALIOP ranges to respectve TROPO pixel: %10.4f" %elapsed_time)
    return saveit_trop0

def f_save_cal_trop(file1,file2,var1,var2):
        #####  now save indexes in output file
    # Create DIMENSIONS to save
    dimcal_1=var1.shape[0]
    dimcal_2=var1.shape[1]
    dimtro_1=var2.shape[0]
    dimtro_2=var2.shape[1]


    # Ancillary information

    DESCRIPTION = """
    (Edited Jul-04-2024)

    File with indexes that map TROPO and CALIOP pixels into each other orbits
    TROPO_filename    = filename of original CALIPSO orbit used
    CALIOP_filename   = filename of original TROPOMAER orbit used
    cali_mapped2_trop = [Ncali x 2] where [i,1&2] are line and row numbers of pixel in
        original TROPO file that contains the  CALIOP profile i.
    trop_mapped2_cali = [Ntrop x 2] where [i,1&2] are the start and end index of the CALIOP
        profiles that are inside the TROPO pixel i

    Ncali are the number of pixels (~50K) that were read from file output.nc
        (from test_read_lat_calipso_v7.py)
    Ntrop are number of scan lines in TROPO file (about 3157)

    indexes in trop_mapped2_cali refer to pixels in output.nc, not to original caliop file
    indexes in cali_mapped2_trop do map to original TROPO file.


    """

    ###  make output file : Date : TropoORBIT and time and CALIOP time
    filenameout="collocation_"+file1[33:42]+'_TRO'+file1[43:49]+'_CAL'+file2[37:45]

    # Create a new NetCDF file
    with nc.Dataset(filenameout+'.nc', 'w', format='NETCDF4') as ncfileout:

        # Create dimensions
        ncfileout.createDimension('caldim_1', dimcal_1)
        ncfileout.createDimension('caldim_2', dimcal_2)
        ncfileout.createDimension('trodim_1', dimtro_1)
        ncfileout.createDimension('trodim_2', dimtro_2)

        # Create variables
        var2d_int1 = ncfileout.createVariable('cali_mapped2_trop', np.int32, ('caldim_1', 'caldim_2'))
        var2d_int2 = ncfileout.createVariable('trop_mapped2_cali', np.int32, ('trodim_1', 'trodim_2'))

        # Assign data to variables
        var2d_int1[:] = var1
        var2d_int2[:] = var2

        # Add ancfileoutillary information as global attributes
        ncfileout.TROPO_filename = filename_tro
        ncfileout.CALIOP_filename = filename_cali
        ncfileout.README = DESCRIPTION

        # Add creation date as an attribute
        ncfileout.creation_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    print('\n NetCDF file ' + filenameout + ' has been created successfully.')
    print('\n File Contents')
    print(DESCRIPTION)
    return
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
   # filename='TROPOMI-Sentinel-5P_L2-TROPOMAER_2019m0814t095613-o09506_v01-2021m1015t022352.nc'
   #    syn_time = datetime(2008,6,30,0,0,0)
   syn_time = datetime(2019,8,14,13,19,12)
   Files = pthin+filename_tro
   # print(Files)
   # Files = granules('/discover/nobackup/dao_ops/intermediate/flk/modis','MOD04',syn_time,coll='006')
   # select group of data from dictionary defined begining of code
   keys=[1,2,3] # select GEODATA amd SCIDATA, using index  start in 1


   TROP = TROPO(Files,SDS,Verbose=1,only_good=True)
   # print(dir(TROP))
   tro_lat =TROP.lat
   tro_lat4=TROP.lat4
   tro_lon =TROP.lon
   tro_lon4=TROP.lon4
   # sys.exit()
   #### Get Caliop
   print("Get CALIOP data")
   netcdf_file_path = "output_Calc.nc"
   # Specify the array names to load
   var_names = ['lat', 'lon','zkm','timeUTC', 'tab532','filename']
   CALI = load_netcdf_arrays(netcdf_file_path,var_names)
   ## get original filename
   print('\n  Loaded data from original file:       \n',CALI['filename'].data.tobytes().decode('utf-8'))
   filename_cali=CALI['filename'].data.tobytes().decode('utf-8')
   if CALI:
    print("\n    CALIPSO arrays loaded in structure CALI:")
    for array_name, array_data in CALI.items():
       print(f"    {array_name}: Shape {array_data.shape}")


   ### Done loading . Now both Tropo and Caliop are in memory.


   ### Coarse Search
   ### Now get 6 points representative from CALIOP
   print('\nCoarse Search ..... ')
   cal_lat_data=CALI['lat'].data
   cal_lon_data=CALI['lon'].data[:,0]
   Ncal=len(CALI['lat'])
   print('  N profiles in CALIOP curtain: ', Ncal)
   # cal_ind_coa=np.append(np.arange(0, Ncal , int(Ncal/8)),Ncal-1)
   cal_ind_coa=np.arange(0, Ncal , int(Ncal/8))
   ### now get these cal lat/lon and use geometric distance (not actual distance)
   ### to find the range of lines and cols., actuall we need to the cols. We already know
   ### the lines span the whole range.
   cal_lat_range=[cal_lat_data.min(),cal_lat_data.max()]
   print("  CALIOP min and max latitudes: %8.4f,%8.4f" %(cal_lat_range[0],cal_lat_range[1]))
   trop_lat_range=[tro_lat.min(),tro_lat.max()]
   print("  TROPOMI min and max latitudes: %8.4f,%8.4f" %(trop_lat_range[0],trop_lat_range[1]))
   if cal_lat_range[1] > trop_lat_range[1] :
    print("        MAX CALIOP latitude farther northe than MAX TROPO latitude")
    sys.exit()

   ### now Select
   lat_tmp=cal_lat_data[cal_ind_coa]
   ### compute distance between TROPOMI and Selected CALIOP points.
   ### Then save respective lines and cols from TRopomi
   min_save=np.empty((len(cal_ind_coa), 2), dtype=int)
   for i,indx in enumerate(cal_ind_coa):
      point=[cal_lat_data[indx],cal_lon_data[indx]]
     #### Compute min lat/lon distance
      tmp=np.square(tro_lat-point[0])+np.square(tro_lon-point[1])
      # Find the minimum value
      min_value = np.min(tmp)
      min_row, min_col = np.unravel_index(np.argmin(tmp), tmp.shape)
      min_save[i,:]=[min_row, min_col]
      # min_row, min_col = min_indices
      # print(f"    Minimum value: {min_value}")
      # print(f"    Row of first minimum: {min_row}")
      # print(f"    Column of first minimum: {min_col}")

   scan = [min_save[:,0].min(),min_save[:,0].max()]
   cols = [min_save[:,1].min(),min_save[:,1].max()]
   if scan[1]> 5000 :sys.exit("Number of scan lines found do not make sense. Stop")
   print('  TROPO Scan lines min/max : ', scan)
   print('  TROPO Column     min/max : ', cols)

   print('\nFine Search....')
   print('  Scan over each Caliop profile and match it to each TROPO pixel .. ')
    ####### Fine Search
    ### switch to fine search for each Caliop pixel
    #### First subset Tropomi lat/lon
    #### Check for outliers
   if cols[1] >450 : cols[1]=450
    #### Subset
   tro_sub_lat=tro_lat[scan[0]:scan[1]+1,cols[0]:cols[1]+1]
   tro_sub_lon=tro_lon[scan[0]:scan[1]+1,cols[0]:cols[1]+1]
   cal_ind_fine=np.arange(0, Ncal )
    ####  Loop over each CALIOP point , find shortest distance
    #### and save each tropomi pixel to which this CALIOP pixel is assigned
    ####
   start_time = time.perf_counter()
    ## Initialize a few things before loop. This array will save TROPOMI granule
    ## scan line that contains the cal_ind_fine[ical] profile
   saveit_cal=-999*np.ones((Ncal, 2), dtype=int)
   tro_shift_lines=20
    ### optimize the loop by making the distance computation with a smaller TROPO
    ### array (from 3500 lines to 20) and use the fact that the index from the previous iteration
    ### can be used in current iteration to subset the TROPO lat/lon grid
   start_line=0
   Nlines_sub=scan[1]-scan[0]+1
   for ical in cal_ind_fine:
      cal0=[cal_lat_data[ical],cal_lon_data[ical]]
      if start_line < Nlines_sub+1:
         #### Subset TROPOMI array to a smaller array
         #### Check if first , last and too close to the end by tro_shift_lines
         if start_line==0:
            trlat_small=tro_sub_lat;
            trlon_small=tro_sub_lon
         elif start_line < Nlines_sub - tro_shift_lines:
            trlat_small=tro_sub_lat[start_line:start_line+tro_shift_lines,:]
            trlon_small=tro_sub_lon[start_line:start_line+tro_shift_lines,:]
         else:
            trlat_small=tro_sub_lat[start_line:-1,:]
            trlon_small=tro_sub_lon[start_line:-1,:]
      #### Compute actual distance
         dist=haversine(trlat_small,trlon_small,cal0[0],cal0[1])
         tmp=np.unravel_index(np.argmin(dist), dist.shape)
         if ical==0:
           saveit_cal[ical,:]=[tmp[0]+scan[0],tmp[1]+cols[0]]
         else:
            saveit_cal[ical,:]=[tmp[0]+saveit_cal[ical-1,0],tmp[1]+cols[0]]
            start_line=saveit_cal[ical,0]-scan[0]
      else:
            print("Reached last line of TROPOMI daylight data. Ending matching process ",ical)
      if (ical % 10000) == 0 : print('   N = %8i of %8i'  %(ical,Ncal))


    # End timer
   end_time = time.perf_counter()
    # Calculate elapsed time
   elapsed_time = end_time - start_time
   print("  Secs. Elapsed time for mapping each Caliop pixe with TROPO: %10.4f" %elapsed_time)

    #### now loop over TROPO pixels and find CALIOP pixel inside each one
   saveit_trop=f_find_cal_in_tro_pix(scan,saveit_cal)
    #### save indexes in output file
   f_save_cal_trop(filename_tro,filename_cali,saveit_cal,saveit_trop)
    #### Now save TROPO data along the CALIOP track
    #### Save lline,row,lat,long,SurfaceType,Reflectance,uvai,ztro,ref,algQA

    ## First get unique values of TROPO lines and rows
    # Convert to a structured array to handle uniqueness correctly
   dtype = np.dtype([('f0', saveit_cal.dtype), ('f1',saveit_cal.dtype)])
   structured_VAR = saveit_cal.view(dtype)
    # Use np.unique to find unique rows
   unique_structured_VAR = np.unique(structured_VAR)
    # Convert back to a regular 2D array
   trop_unique_indx = unique_structured_VAR.view(saveit_cal.dtype).reshape(-1, 2)

    # Now loop over each pair, gather TROP data and save it into output
   Nout=len(trop_unique_indx)
   out_array1= -999*np.ones((Nout, 3), dtype=int) # line,row,algQA,ymd,hms
   out_array2= -999.*np.ones((Nout, 8), dtype=float)#lat,long,SurfaceType,Reflectance(1:2),uvai,ztro

   #### Now gather all data according to unique lines,rows in trop_unique_indx
   for i,v in enumerate(tuple(trop_unique_indx)):
      out_array1[i,:]=[v[0],v[1],TROP.algQA[v[0],v[1]]]
      out_array2[i,:]=[TROP.dtime[0],\
             TROP.lat[v[0],v[1]],\
             TROP.lon[v[0],v[1]],\
             TROP.surftype[v[0],v[1]],\
             TROP.uvai[v[0],v[1]],\
             TROP.ztro[v[0],v[1]],\
             TROP.ref[v[0],v[1],0],\
             TROP.ref[v[0],v[1],1]\
            ]

    #### Now saved


    #####  now save indexes in output file
   # Create DIMENSIONS to save
   dim_Lines=Nout
   dim_Cols1=out_array1.shape[1]
   dim_Cols2=out_array2.shape[1]


   # Ancillary information

   DESCRIPTION = """
   (Edited Jul-05-2024)

   File with TROPOMI L2 OMAERUV retrievals in pixels that were profiled by CALIPSO
   There are 2 arrays.
   array1=[tro_line,tro_row,algQA,delta_time]
   array2=[lat,long,SurfaceType,Reflectance(1:2),uvai,ztro]

   Including of original CALIOP and TROP files .
   NOTE :
    delta_time=time measurement of each scan line in milliseconds
    since 2019-08-14 00:00:00
    Reflectance at the two retrieval bands are reported

   """

   ###  make output file : Date : TropoORBIT and time and CALIOP time
   filenameout= 'tropoin_CALtrack__' + filename_tro[33:42] + '_TRO'+filename_tro[43:49]+'_CAL'+filename_cali[37:45]

   # Create a new NetCDF file
   with nc.Dataset(filenameout+'.nc', 'w', format='NETCDF4') as ncfileout:

       # Create dimensions
      ncfileout.createDimension('Nlines', dim_Lines)
      ncfileout.createDimension('dim_Cols1', dim_Cols1)
      ncfileout.createDimension('dim_Cols2', dim_Cols2)

       # Create variables
      var2d_int1 = ncfileout.createVariable('trop_int', int, ('Nlines', 'dim_Cols1'))
      var2d_int2 = ncfileout.createVariable('trop_float',float, ('Nlines', 'dim_Cols2'))

       # Assign data to variables
      var2d_int1[:] = out_array1
      var2d_int2[:] = out_array2

       # Add ancfileoutillary information as global attributes
      ncfileout.TROPO_filename = filename_tro
      ncfileout.CALIOP_filename = filename_cali
      ncfileout.README = DESCRIPTION

       # Add creation date as an attribute
      ncfileout.creation_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

   print('\n NetCDF file ' + filenameout + ' has been created successfully.')
   print('\n File Contents')
   print(DESCRIPTION)
sys.exit()
#####
# Plot
   # Parameters for UVAI - Full orbit
   vmin = 1;vmax = 3
   num_colors = 10
   bad_value = TROP.uvai.min()  #

   ## select data to plot:
   lat_omi=TROP.lat;lon_omi=TROP.lon;uvai=TROP.uvai

   # Handle bad values
   uvai[TROP.uvai == bad_value] = np.nan

   # Create discrete colormap
   base = plt.cm.get_cmap('jet')
   color_list = base(np.linspace(0, 1, num_colors))
   cmap = ListedColormap(color_list, name='jet' + str(num_colors))
   cmap.set_under('lightgray')  # Color for values below vmin
   cmap.set_over('brown')  # Color for values above vmax
   cmap.set_bad('white')  # Color for values bad values

   # Define boundaries for discrete levels
   boundaries = np.linspace(vmin, vmax, num_colors + 1)
   norm = BoundaryNorm(boundaries, cmap.N)
   # norm = BoundaryNorm(np.linspace(vmin, vmax, num_colors + 1), cmap.N)



   # Plot
   plt.figure(figsize=(14, 6))
   ax = plt.axes(projection=ccrs.PlateCarree())
   ax.add_feature(cfeature.COASTLINE);ax.add_feature(cfeature.BORDERS)

   # Using scatter for plotting points
   sc = plt.scatter(lon_omi, lat_omi, c=uvai, cmap=cmap, norm=norm,
                    s=5, edgecolor='none', transform=ccrs.PlateCarree())
   # Create colorbar with adjusted width and position
   cbar = plt.colorbar(sc, orientation='horizontal', extend='both', shrink=0.7, pad=0.05, ticks=boundaries)
   cbar.ax.minorticks_off()  # Remove minor ticks
   cbar.ax.set_xlabel('UV Aerosol Index')
   # cbar.ax.set_xticklabels([f'{val:.1f}' for val in boundaries])
   # Set colorbar ticks at boundaries
   cbar.set_ticks(boundaries)

   # Add gridlines and labels
   gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray')
   gl.top_labels = False;gl.right_labels = False


   ### now add line for CALIOP Track
   plt.plot(CALI['lon'][::100],CALI['lat'][::100],color='black')
   # plt.xlim(-180, 180);plt.ylim(-60, 60)
   plt.xlim(-160, 160);plt.ylim(-60, 60)
   plt.title('TROPOMI Aerosol Index \n' + filename_tro[23:-24] ,fontsize=16)


###### NOW display Tropomi Data alon the CALIOP track
   #### Subset OMI data that includes CALIOP track
   #### get line and column numbers

   # Plot
   plt.figure(figsize=(6, 12))
   ax = plt.axes(projection=ccrs.PlateCarree())
   ax.add_feature(cfeature.COASTLINE);ax.add_feature(cfeature.BORDERS)

   u=unique_structured_VAR;# type u in the command line to see the fields f0 and f1 in this structure
   xlon=lon_omi[u['f0'],u['f1']];ylat=lat_omi[u['f0'],u['f1']]
   cuvai=uvai[u['f0'],u['f1']]
   sc = plt.scatter(xlon[::2],ylat[::2],c=cuvai[::2],cmap=cmap, norm=norm,marker='s',s=4,\
                    transform=ccrs.PlateCarree())
   # Create colorbar with adjusted width and position
   cbar = plt.colorbar(sc, orientation='vertical', extend='both', shrink=0.7, pad=0.05, ticks=boundaries)
   cbar.ax.minorticks_off()  # Remove minor ticks
   cbar.ax.set_ylabel('UV Aerosol Index')
   # cbar.ax.set_xticklabels([f'{val:.1f}' for val in boundaries])
   # Set colorbar ticks at boundaries
   cbar.set_ticks(boundaries)
   # plt.xlim(np.ceil(xlon.min()), np.floor(xlon.max()));
   plt.xlim(-40, 20);
   # Add gridlines and labels
   gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray')
   gl.top_labels = False;gl.right_labels = False

   plt.ylim(-60, 60)
   plt.title('TROPOMI Aerosol Index \n' + filename_tro[23:-24] ,fontsize=12)