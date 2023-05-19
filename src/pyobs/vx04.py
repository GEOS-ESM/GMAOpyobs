"""
Reads Level 2 Deep Blue and Dark Target retrievals applied to VIIRS
granules for a single day and returns a
single object with the relevant data.

Builds on the heritage of mxd04.py

This software is hereby placed in the public domain.
patricia.castellanos@nasa.gov
"""

import os
import sys
import numpy as np
from datetime import date, datetime, timedelta
from glob     import glob
from dateutil.parser import isoparse
from pyobs.npz import NPZ
from binObs_   import binobs2d, binobs3d, binobscnt3d

try:
    from pyods import ODS # must be inported before pyhdf
except:
    pass

from netCDF4 import Dataset
from bits import BITS

#---  

DATE_START = datetime(1993,1,1,0,0,0)

SDS = dict (
     DT_META = ('longitude', 'latitude',
              'sensor_zenith_angle', 'sensor_azimuth_angle',
              'solar_zenith_angle', 'solar_azimuth_angle',
              'Scattering_Angle', 'Glint_Angle'),
     DT_LAND = ( 'Mean_Reflectance_Land',
               'Corrected_Optical_Depth_Land',
               'Surface_Reflectance_Land',
               'Aerosol_Cloud_Fraction_Land',
               'Land_Ocean_Quality_Flag',
               'Number_Pixels_Used_Land',
               'STD_Reflectance_Land'),
     DT_OCEAN = ( 'Effective_Optical_Depth_Average_Ocean',
               'Optical_Depth_Small_Average_Ocean',
               'Effective_Radius_Ocean',
               'Asymmetry_Factor_Average_Ocean',
               'Angstrom_Exponent_1_Ocean',
               'Angstrom_Exponent_2_Ocean',
               'Aerosol_Cloud_Fraction_Ocean',
               'Mean_Reflectance_Ocean',
               'Land_Ocean_Quality_Flag' ),
     DB_META = ('Longitude', 'Latitude', 'Scan_Start_Time',
              'Viewing_Zenith_Angle', 'Relative_Azimuth_Angle',
              'Solar_Zenith_Angle',
              'Scattering_Angle' ),
     DB_LAND = ('Aerosol_Optical_Thickness_550_Land',
               'Spectral_Aerosol_Optical_Thickness_Land',
               'Spectral_Single_Scattering_Albedo_Land',
               'Spectral_Surface_Reflectance',
               'Spectral_TOA_Reflectance_Land',
               'Aerosol_Optical_Thickness_550_Land_Best_Estimate',
               'Aerosol_Optical_Thickness_550_STDV_Land',
               'Aerosol_Optical_Thickness_QA_Flag_Land',
               'Algorithm_Flag_Land',
               'Aerosol_Type_Land',
               'Number_Of_Pixels_Used_Land',
               'Number_Valid_Pixels',
               'Unsuitable_Pixel_Fraction_Land_Ocean',
               'TOA_NDVI',
               'Total_Column_Ozone',
               'Cell_Average_Elevation_Land',
               'Precipitable_Water'),
     DB_OCEAN = ('Aerosol_Optical_Thickness_550_Ocean', 
               'Spectral_Aerosol_Optical_Thickness_Ocean',
               'Spectral_TOA_Reflectance_Ocean',
               'Aerosol_Optical_Thickness_550_Ocean_Best_Estimate',
               'Aerosol_Optical_Thickness_550_STDV_Ocean',
               'Aerosol_Optical_Thickness_QA_Flag_Ocean',
               'Algorithm_Flag_Ocean',
               'Aerosol_Type_Ocean',
               'Number_Of_Pixels_Used_Ocean',
               'Number_Valid_Pixels',
               'Unsuitable_Pixel_Fraction_Land_Ocean',
               'Wind_Speed',
               'Wind_Direction',
               'Cell_Average_Elevation_Ocean',
               'Total_Column_Ozone',
               'Precipitable_Water')
        )
# NOTE: DEEP BLUE does not have cloud information in their files.



# AOD Channels
CHANNELS = dict (
                   DT_LAND = ( 480., 550., 670., 2250.),
                   DT_OCEAN = ( 480., 550., 670., 860., 1240., 1600., 2250. ),
                   DB_LAND = ( 412, 488, 550, 670 ), # file has 550 separate
                   DB_OCEAN = (488.,  550.,  670.,  865., 1240., 1640., 2250.),
                   DT_SREF = ( 480., 670., 2250. ),
                   DB_SREF = (412., 488., 670. ),
                )

ALIAS = dict (  Longitude = 'lon',
                Latitude = 'lat',
                longitude = 'lon',
                latitude  = 'lat',
                Viewing_Zenith_Angle = 'SensorZenith',
                sensor_zenith_angle  = 'SensorZenith',
                sensor_azimuth_angle = 'SensorAzimuth',
                Relative_Azimuth_Angle = 'RelativeAzimuth',
                Solar_Zenith_Angle = 'SolarZenith',
                solar_zenith_angle = 'SolarZenith',
                solar_azimuth_angle = 'SolarAzimuth',
                Scattering_Angle = 'ScatteringAngle',
                Glint_Angle = 'GlintAngle',
                Mean_Reflectance_Land = 'reflectance',
                Surface_Reflectance_Land = 'sfc_reflectance',
                Corrected_Optical_Depth_Land = 'aod',
                Aerosol_Cloud_Fraction_Land = 'cloud',
                Effective_Optical_Depth_Average_Ocean = 'aod',
                Optical_Depth_Small_Average_Ocean = 'aod_fine',
                Aerosol_Cloud_Fraction_Ocean = 'cloud',
                Mean_Reflectance_Ocean = 'reflectance',
                Spectral_Aerosol_Optical_Thickness_Land = 'aod3ch',
                Aerosol_Optical_Thickness_550_Land = 'aod550',
                Spectral_Surface_Reflectance = 'sfc_reflectance',
                Spectral_TOA_Reflectance_Land = 'reflectance',
                TOA_NDVI = 'NDVI',
                Spectral_Single_Scattering_Albedo_Land = 'ssa',
                Algorithm_Flag_Land = 'algflag',
                Aerosol_Type_Land = 'atype',
                Angstrom_Exponent_Land = 'angstrom',
                Spectral_Aerosol_Optical_Thickness_Ocean = 'aod',
                Aerosol_Optical_Thickness_550_Ocean = 'aod550',
                Spectral_TOA_Reflectance_Ocean = 'reflectance',
                Algorithm_Flag_Ocean = 'algflag',
                Aerosol_Type_Ocean = 'atype',
                Angstrom_Exponent_Ocean = 'angstrom',                
                Number_Of_Pixels_Used_Ocean = 'npixels_used',
                Number_Of_Pixels_Used_Land = 'npixels_used',
                Number_Valid_Pixels = 'npixels_valid',
                Aerosol_Optical_Thickness_QA_Flag_Land = 'qa_flag',
                Aerosol_Optical_Thickness_QA_Flag_Ocean = 'qa_flag',
                Land_Ocean_Quality_Flag = 'qa_flag',
                Scan_Start_Time = 'Time',
                Cell_Average_Elevation_Land = 'pixel_elevation',
             )

BAD, MARGINAL, GOOD, BEST = ( 0, 1, 2, 3 ) # DT QA marks
# for DB 0 = no retrieval, 1 = poor, 2 = moderate, 3 = good

translate_sat = {'Suomi-NPP': 'SNPP'}


KX = dict ( SNPP_DT_OCEAN = 401,
            SNPP_DT_LAND  = 402,
            SNPP_DB_OCEAN  = 403,
            SNPP_DB_LAND  = 404, 
          )

KT = dict ( AOD = 45, )

IDENT = dict ( SNPP_DT_OCEAN = 'vsnppdto',
               SNPP_DT_LAND  = 'vsnppdtl',
               SNPP_DB_OCEAN  = 'vsnppdbo',
               SNPP_DB_LAND  = 'vsnppdbl',
          )

MISSING = 999.999

#...........................................................................

class Vx04_L2(object):
    """
    This class implements the VIIRS Level 2 AEROSOL products.
    VIIRS currently flies on the Suomi NPP (SNPP) and NOAA-20 (FKA JPSS-1) satellites
    VIIRS will be flown on JPSS-2,3, and 4
    The Level2 AEROSOL products are implementations of the MODIS Deep Blue and Dark Target algoritms,
    referred to as MOD04 (TERRA satellite) and MYD04 (AQUA satellite).
    Therefore, we will refer to these products at VSNPP04 and VN2004, respectively
    """

    def __init__ (self,Path,algo,syn_time=None,nsyn=8,Verb=0,
                  only_good=True,SDS=SDS,alias=None,anet_wav=False):
       """
       Reads individual granules or a full day of Level 2 Vx04 files
       present on a given *Path* and returns a single object with
       all data concatenated for a given algorithm. On input, 

       Required parameters:
         Path -- can be a single file, a single directory, of a list
                 of files and directories.  Directories are
                 transversed recursively. If a non Vx04 Level 2
                 file is encountered, it is simply ignored.
         algo -- Algorithm: DT_LAND, DT_OCEAN, DB_LAND or DB_OCEAN

       Optional parameters:
         syn_type  --- synoptic time
         nsyn      --- number of synoptic times per day
         only_good --- keep only *good* observations
         Verb      -- Verbose level:
                 0 - really quiet (default)
                 1 - Warns if invalid file is found
                 2 - Prints out non-zero number of aerosols in each file.
         SDS      --- Variables to be read from L2 Aerosol files.  Must 
                      be a dictionary with keys '{Algo}_META' and '{Algo}_{Surface}'
         ALIAS    --- dictionary of alises for SDSs
         anet_wav --- angstrom interpolate retrieved AOD to common AERONET wavelengths
                      e.g. 480 --> 490, 1600 --> 1640, 860 --> 870


       """

       if algo not in ('DT_LAND', 'DT_OCEAN', 'DB_LAND', 'DB_OCEAN'):
           raise ValueError, "invalid algorithm "+algo+" --- must be DT_LAND, DT_OCEAN, DB_LAND, DB_OCEAN"

#      Initially are lists of numpy arrays for each granule
#      ------------------------------------------------
       self.verb = Verb
       self.sat  = None # Satellite name
       self.col  = None # collection, e.g., 011
       self.algo = algo
       Algo, Surface = algo.split('_')
       self.SDS = SDS['{}_META'.format(Algo)] + SDS[algo]
       self.SDS_META = SDS['{}_META'.format(Algo)]

       # Add/Substitute some aliases if given
       # ------------------------------------
       self.ALIAS = ALIAS.copy()
       if alias is not None:
           for a in alias: self.ALIAS[a] = alias[a]  
       

       # Create empty lists for SDS to be read from file
       # -----------------------------------------------
       for name in self.SDS:
           self.__dict__[name] = []

       # Read each granule, appending them to the list
       # ---------------------------------------------
       if type(Path) is list:
           if len(Path) == 0:
               self.nobs = 0
               print "WARNING: Empty Vx04_L2 object created"
               return
       else:
           Path = [Path, ]
       self._readList(Path)

       #Protect against empty VX04 files
       # --------------------------------
       if len(self.Scattering_Angle) == 0:
           self.nobs = 0
           print "WARNING: Empty MxD04_L2 object created"
           return           

       # Make each attribute a single numpy array
       # ----------------------------------------
       if 'DT' in self.algo:
           self.SDS += ('Scan_Start_Time',)
       for sds in self.SDS:
           try:
               self.__dict__[sds] = np.ma.concatenate(self.__dict__[sds])
           except:
               print "Failed concatenating "+sds

       
       # Determine index of "good" observations
       # --------------------------------------
       if self.algo == 'DT_LAND':
           self.iGood = (self.Land_Ocean_Quality_Flag == BEST) & (~self.Corrected_Optical_Depth_Land.mask[:,1])
       elif self.algo == 'DT_OCEAN':
           self.iGood = (self.Land_Ocean_Quality_Flag > BAD) & (~self.Effective_Optical_Depth_Average_Ocean.mask[:,1])
       elif self.algo == 'DB_LAND':
           self.iGood = self.Aerosol_Optical_Thickness_QA_Flag_Land > BAD # for now
       elif self.algo == 'DB_OCEAN':
           self.iGood = self.Aerosol_Optical_Thickness_QA_Flag_Ocean == BEST
       else:
           raise ValueError, 'invalid algorithm (very strange)'


       # Keep only "good" observations
       # -----------------------------
       if only_good:
           m = self.iGood
           for sds in self.SDS:
               rank = len(self.__dict__[sds].shape)
               if rank == 1:
                   self.__dict__[sds] = self.__dict__[sds][m]
               elif rank == 2:
                   self.__dict__[sds] = self.__dict__[sds][m,:]
               else:
                   raise IndexError, 'invalid rank=%d'%rank
           self.iGood = self.iGood[m]

       
       # Make aliases for compatibility with older code 
       # ----------------------------------------------
       Alias = self.ALIAS.keys()
       for sds in self.SDS:
           if sds in Alias:
               self.__dict__[self.ALIAS[sds]] = self.__dict__[sds] 

       # Calculate Glint angle for Deep Blue
       # ----------------------
       if 'DB' in self.algo:
           sza = np.radians(self.SolarZenith)
           vza = np.radians(self.SensorZenith)
           raa = np.radians(self.RelativeAzimuth)
           cglint = np.cos(sza)*np.cos(vza) + np.sin(sza)*np.sin(vza)*np.cos(raa)
           self.GlintAngle = np.degrees(np.arccos(cglint))
           self.SDS += ('GlintAngle',)
       elif 'DT' in self.algo:
           # this kind of seems to match the DB RAA
           # the DT sensor and solar azimuth angles have a bug
           # this needs to be rechecked with updated version of DT dataset
           # DB definition for RAA is 'Gordon Convention', RAA=0 is back scattering
           saa = self.SolarAzimuth
           vaa = self.SensorAzimuth
           i = saa < 0
           saa[i] = 360. + saa[i]
           i = vaa < 0
           vaa[i] = 360. + vaa[i]
           raa = saa - vaa + 180.
           i = raa > 180.
           raa[i] = 360. - raa[i]

           self.RelativeAzimuth = raa
           self.SDS += ('RelativeAzimuth',)

       # Create corresponding python time
       # --------------------------------
       if 'DB' in self.algo:
           self.Time = np.array([DATE_START+timedelta(seconds=s) for s in self.Scan_Start_Time])
       else:
           self.Time = np.array(self.Time)   # masked datetime arrays aren't friendly

       # ODS friendly attributes
       # -----------------------
       self.nobs = self.Scattering_Angle.shape[0]
       self.kx = KX[self.sat+'_'+self.algo]
       self.ident = IDENT[self.sat+'_'+self.algo]
       self.channels = np.array(CHANNELS[self.algo])
       if Surface == 'LAND':
           self.sChannels = CHANNELS["{}_SREF".format(Algo)]   # LAND surface reflectivity (not the same as algo)           

       if 'DB' in self.algo:
           self.rChannels = self.Reflectance_Bands  # [ 412.,  488.,  550.,  670.,  865., 1240., 1640., 2250.]
       elif self.algo == 'DT_LAND':
           self.rChannels = np.array([480.,670.,2250.])
       elif self.algo == 'DT_OCEAN':
           self.rChannels = np.array([480.,550.,670.,860.,1240.,1600.,2250.])

       # Concatenate AOD channels for Deep Blue
       # --------------------------------------
       if self.algo == 'DB_LAND':
           try:
               self.aod = np.ones((self.nobs,4))
               self.aod[:,0] = self.aod3ch[:,0]
               self.aod[:,1] = self.aod3ch[:,1]
               self.aod[:,2] = self.aod550[:]
               self.aod[:,3] = self.aod3ch[:,2]
           except:
               pass # don't fuss, aod3ch may not have been read


       # Angstrom interpolate to AERONET wavelengths for ease of comparison
       # and to faciliate NNR ODS files
       # DT: 480 --> 470, 860 --> 870, 1600 --> 1640
       # DB: 488 --> 470, 865 --> 870
       # ------------------------------------------------------------------
       self.anet_wav = anet_wav
       if anet_wav:
           if 'DT' in self.algo:
               tch = 470

               ch = 550
               ich = np.argmin(np.abs(self.channels - ch))
               aod550 = self.aod[:,ich]
               ch = 480
               ich = np.argmin(np.abs(self.channels - ch))
               aod480 = self.aod[:,ich]

               aod_ = self.aodInterpAngs(tch,aod480,aod550,480,550)
               self.channels[ich] = tch
               self.aod[:,ich] = aod_

           elif self.algo == 'DT_OCEAN':
               tch = 870

               ch = 670
               ich = np.argmin(np.abs(self.channels - ch))
               aod670 = self.aod[:,ich]
               ch = 860
               ich = np.argmin(np.abs(self.channels - ch))
               aod860 = self.aod[:,ich]               

               aod_ = self.aodInterpAngs(tch,aod670,aod860,670,860)
               self.channels[ich] = tch
               self.aod[:,ich] = aod_

               tch = 1640

               ch = 2250
               ich = np.argmin(np.abs(self.channels - ch))
               aod2250 = self.aod[:,ich]
               ch = 1600
               ich = np.argmin(np.abs(self.channels - ch))
               aod1600 = self.aod[:,ich]               

               aod_ = self.aodInterpAngs(tch,aod1600,aod2250,1600,2250)
               self.channels[ich] = tch
               self.aod[:,ich] = aod_               

           elif 'DB' in self.algo:
               tch = 470
               
               ch = 550
               ich = np.argmin(np.abs(self.channels - ch))
               aod550 = self.aod[:,ich]
               ch = 488
               ich = np.argmin(np.abs(self.channels - ch))
               aod488 = self.aod[:,ich]               

               aod_ = self.aodInterpAngs(tch,aod488,aod550,488,550)
               self.channels[ich] = tch
               self.aod[:,ich] = aod_               

           elif self.algo == 'DB_OCEAN':
               tch = 870

               ch = 670
               ich = np.argmin(np.abs(self.channels - ch))
               aod670 = self.aod[:,ich]
               ch = 865
               ich = np.argmin(np.abs(self.channels - ch))
               aod865 = self.aod[:,ich]

               aod_ = self.aodInterpAngs(tch,aod670,aod865,670,865)
               self.channels[ich] = tch
               self.aod[:,ich] = aod_
               


       # ODS time variables
       # -------------------
       if syn_time == None:
           self.syn_time = None
           self.time = None
           self.nymd = None
           self.nhms = None
           self.nsyn = None
       else:
           Dt = [ t-syn_time for t in self.Time ]
           self.time = np.array([ (86400.*dt.days+dt.seconds)/60. for dt in Dt ])   # in minutes
           self.syn_time = syn_time
           self.nymd = 10000 * syn_time.year + 100*syn_time.month  + syn_time.day
           self.nhms = 10000 * syn_time.hour + 100*syn_time.minute + syn_time.second
           self.nsyn = nsyn # number of synoptic times per day


       # Create a pseudo cloud fraction for Deep Blue
       if Algo == 'DB':
           self.cloud = 1. - self.npixels_used.astype(float)/self.npixels_valid.astype(float)


#---
    def aodInterpAngs(self,lambda_,tau1,tau2,lambda1,lambda2):
        """
           Angstrom-interpolated AOD.
           lambda_ = target wavelength for AOD output
        """
        tau = MISSING * np.ones(len(tau1))
        I = (tau1>0) & (tau2>0)
        if any(I):
            angstrom = -np.log(tau1[I]/tau2[I])/np.log(float(lambda1)/float(lambda2))
            tau[I] = tau2[I] * (lambda2/lambda_)**angstrom
        return tau

#---
    def _readList(self,List):
        """
        Recursively, look for files in list; list items can
        be files or directories.
        """
        for item in List:
            if os.path.isdir(item):      self._readDir(item)
            elif os.path.isfile(item):   
                if 'DB' in self.algo:
                    self._readGranuleDB(item)
                else:
                    self._readGranuleDT(item)
            else:
                print "%s is not a valid file or directory, ignoring it"%item
#---
    def _readDir(self,dir):
        """Recursively, look for files in directory."""
        for item in os.listdir(dir):
            path = dir + os.sep + item
            if os.path.isdir(path):      self._readDir(path)
            elif os.path.isfile(path):   self._readGranule(path)
            else:
                print "%s is not a valid file or directory, ignoring it"%item

#---
    def _readGranuleDB(self,filename):
        """Reads one Vx04 Deep Blue granule with Level 2 aerosol data."""

        # Don't fuss if the file cannot be opened
        # ---------------------------------------
        try:
            if self.verb:
                print "[] Working on "+filename
            nc = Dataset(filename)
        except:
            if self.verb > 2:
                print "- %s: not recognized as an netCDF file"%filename
            return 

        # Read select variables (reshape to allow concatenation later)
        # ------------------------------------------------------------
        for sds in self.SDS:
            v = nc.variables[sds][:]
            a = nc.variables[sds].ncattrs()
            if 'scale_factor' in a:
                scale = nc.variables[sds].getncattr('scale_factor')
                v = scale*v
            if 'add_offset' in a:
                add = nc.variables[sds].getncattr('add_offset')
                v = v + add

            if len(v.shape) == 3:
                if "TOA_Reflectance" in sds:
                    i, j, k = v.shape
                    v = v.reshape((i*j,k))
                else:
                    i, j, k = v.shape
                    v = v.reshape((i,j*k)).T
            elif len(v.shape) == 2:
                v = v.ravel()
            else:
                raise IndexError, "invalid shape for SDS <%s>"%sds
            self.__dict__[sds].append(v) 


#       Satellite name
#       --------------
        if self.sat is None:
            self.sat = translate_sat[nc.platform]
            
#       Collection
#       ----------
        if self.col is None:
            self.col = nc.product_name.split('.')[-3]

#       Reflectance Bands
#       ------------------
        self.Reflectance_Bands = nc.variables['Reflectance_Bands'][:]

#---
    def _readGranuleDT(self,filename):
        """Reads one Vx04 Dark Target granule with Level 2 aerosol data."""

        # Don't fuss if the file cannot be opened
        # ---------------------------------------
        try:
            if self.verb:
                print "[] Working on "+filename
            nc = Dataset(filename)
            data  = nc.groups['geophysical_data']
            loc   = nc.groups['geolocation_data']
        except:
            if self.verb > 2:
                print "- %s: not recognized as an netCDF file"%filename
            return

        # Read select variables (reshape to allow concatenation later)
        # ------------------------------------------------------------
        for sds in self.SDS:
            if sds in self.SDS_META:
                vobj = loc.variables[sds]
            else:
                vobj = data.variables[sds]

            a = vobj.ncattrs()
            v = vobj[:]

            # commenting this out for now
            # all the scale factors in the DT files seem to be wrong
            # maybe it will be fixed in future versions, so keeping
            # code for now
#            if 'scale_factor' in a:
#                scale = vobj.getncattr('scale_factor')
#                v = scale*v
#            if 'add_offset' in a:
#                add = vobj.getncattr('add_offset')
#                v = v + add

            if len(v.shape) == 3:
                i, j, k = v.shape
                v = v.reshape((i*j,k))
            elif len(v.shape) == 2:
                v = v.ravel()
            else:
                raise IndexError, "invalid shape for SDS <%s>"%sds
            self.__dict__[sds].append(v)


#       Satellite name
#       --------------
        if self.sat is None:
            self.sat = translate_sat[nc.platform]

#       Collection
#       ----------
        if self.col is None:
            self.col = nc.product_name.split('.')[-3]

#       Time
        if 'Scan_Start_Time' not in self.__dict__:
            self.Scan_Start_Time = [np.repeat(isoparse(nc.time_coverage_start[:-1]),i*j)]
        else:
            self.Scan_Start_Time.append(np.repeat(isoparse(nc.time_coverage_start[:-1]),i*j))

#---

    def reduce(self,I):
        """
        Reduce observations according to index I. 
        """
        Nicknames = self.ALIAS.values()
        for name in self.__dict__:
            if name in Nicknames:
                continue # alias do not get reduced
            q = self.__dict__[name]
            if type(q) is type(self.lon):
                if len(q) == self.nobs:
                    # print "{} Reducing "+name
                    self.__dict__[name] = q[I]

        Alias = self.ALIAS.keys()
        for sds in self.SDS:
            if sds in Alias:
                self.__dict__[self.ALIAS[sds]] = self.__dict__[sds] # redefine aliases

        self.nobs = len(self.lon)
        # Create corresponding python time
        # --------------------------------
        if 'DB' in self.algo:
            self.Time = np.array([DATE_START+timedelta(seconds=s) for s in self.Scan_Start_Time])
        else:
            self.Time = np.array(self.Time)   # masked datetime arrays aren't friendly        

#---
    def write(self,filename=None,dir='.',expid=None,Verb=1):
        """
        Writes the un-gridded data to a numpy npz file. 
        """

        # Stop here is no good obs available
        # ----------------------------------
        if self.nobs == 0:
            return # no data to work with
        if any(self.iGood) == False:
            return # no good data to work with

        if expid == None:
            expid = self.algo

        if filename is None:
            filename = '%s/%s.viirs.%d_%02dz.npz'%(dir,expid,self.nymd,self.nhms/10000)

        version = 1 # File format version
        meta = [self.nymd,self.nhms,self.nobs,self.nch,self.kx,version]
        savez(filename,
                            meta = meta,
                             lon = self.lon,
                             lat = self.lat,
                              ks = self.ks,
                        channels = self.channels,
                         qa_flag = self.qa_flag,
                     SolarZenith = self.SolarZenith,
                    SensorZenith = self.SensorZenith,
                  RelativeAzimuth = self.RelativeAzimuth,
                 ScatteringAngle = self.ScatteringAngle,
                           cloud = self.cloud,
                             aod = self.aod,
                     reflectance = self.reflectance)

        if Verb >=1:
            print "[w] Wrote file "+filename
#---

    def writeODS(self,filename=None,dir='.',expid=None,channels=None,
                 revised=False,nsyn=8,Verb=1):
        """
        Writes the un-gridded data to an ODS file. If *revised*
        is True, the revised *aod_* parameter is written to file.
        """
        
        if self.syn_time == None:
            raise ValuError, "synoptic time missing, cannot write ODS"
            
        # Stop here if no good obs available
        # ----------------------------------
        if self.nobs == 0:
            return # no data to work with
        if any(self.iGood) == False:
            return # no good data to work with

        if expid == None:
            expid = self.algo

        if filename is None:
            filename = '%s/%s.obs.%d_%02dz.ods'%(dir,expid,self.nymd,self.nhms/10000)

        if channels is None:
            channels = self.channels

        # Create and populated ODS object
        # -------------------------------
        ns = self.nobs
        nobs = len(channels) * ns
        ods = ODS(nobs=nobs, kx=self.kx, kt=KT['AOD'])
        i = 0
        ks = np.arange(ns) + 1
        for ch in channels:
            I = range(i,i+ns)
            j = list(channels).index(ch) # index of channel
            ods.ks[I]  = ks
            ods.lat[I] = self.lat[:]
            ods.lon[I] = self.lon[:]
            ods.time[I] = self.time[:].astype('int')
            ods.lev[I] = ch
            ods.qch[I] = self.qa_flag[:].astype('int')
            if revised:
                ods.obs[I]  = self.aod_[:,j].astype('float32')
                ods.xvec[I] = self.aod[:,j].astype('float32')
            else:
                ods.obs[I] = self.aod[:,j].astype('float32')
            i += ns

        # Handle corrupted coordinates
        # ----------------------------
        iBad = (ods.lon<-180) | (ods.lon>180.) | \
               (ods.lat<-90)  | (ods.lat>90.)  | \
               (abs(ods.time)>1440./self.nsyn)
        ods.lon[iBad] = 0.
        ods.lat[iBad] = -90.
        ods.time[iBad] = 0.
        ods.qcx[iBad] = 2

        # Exclusion flag
        # --------------
        iGood = (ods.qch>0) & (ods.obs<10.) & (ods.qcx==0)
        ods.qcx[:] = 1     # All bad...
        ods.qcx[iGood] = 0 # ... unless good

        ods_ = ods.select(qcx=0)
        if Verb >=1:
            print "[w] Writing file <"+filename+"> with %d observations"%ods_.nobs

        ods_.write(filename,self.nymd,self.nhms,nsyn=8,ftype='pre_anal')
        
#---
    def writeg(self,filename=None,dir='.',expid=None,refine=8,res=None,
               channels=None,Verb=1):
       """
        Writes gridded MODIS measurements to file.

         refine  -- refinement level for a base 4x5 GEOS-5 grid
                       refine=1  produces a   4  x  5    grid
                       refine=2  produces a   2  x2.50   grid
                       refine=4  produces a   1  x1,25   grid
                       refine=8  produces a  0.50x0.625  grid
                       refine=16 produces a  0.25x0.3125 grid
        Alternatively, one can specify the grid resolution with a
        single letter:

         res     -- single letter denoting GEOS-5 resolution,
                       res='a'  produces a   4  x  5    grid
                       res='b'  produces a   2  x2.50   grid
                       res='c'  produces a   1  x1,25   grid
                       res='d'  produces a  0.50x0.625  grid
                       res='e'  produces a  0.25x0.3125 grid

                   NOTE: *res*, if specified, supersedes *refine*.

         Verb -- Verbose level:
                 0 - really quiet (default)
                 1 - Warns if invalid file is found
                 2 - Prints out non-zero number of fires in each file.


       """
       from gfio import GFIO
       
       # Stop here is no good obs available
       # ----------------------------------
       if self.nobs == 0:
           return # no data to work with
       if any(self.iGood) == False:
           return # no good data to work with

       if expid == None:
           expid = self.algo

#      Output grid resolution
#      ----------------------
       if res is not None:
           if res=='a': refine = 1 
           if res=='b': refine = 2
           if res=='c': refine = 4
           if res=='d': refine = 8
           if res=='e': refine = 16

#      Lat lon grid
#      ------------
       dx = 5. / refine
       dy = 4. / refine
       im = int(360. / dx)
       jm = int(180. / dy + 1)

       glon = np.linspace(-180.,180.,im,endpoint=False)
       glat = np.linspace(-90.,90.,jm)

       if channels is None:
           channels = self.channels
 
       levs = np.array(channels)

       nch = len(channels)
       nymd = self.nymd
       nhms = self.nhms

       vtitle = [ 'Aerosol Optical Depth',
                  'Aerosol Optical Depth (Revised)',
                  'Aerosol Optical Depth (Fine Mode)',
                  'Aerosol Optical Depth Obs Count',
                  'Aerosol Optical Depth (Revised) Obs Count',
                  'Cloud Fraction' ]

       vname  = ['tau', 'tau_', 'tau_fine', 'count_tau', 'count_tau_', 'cloud' ]
       vunits = [ '1',    '1',     '1',       '1',            '1',       '1',  ]
       kmvar  = [ nch,    nch,     nch,       nch,            nch,        0    ]

       title = 'Gridded MODIS Aerosol Retrievals'
       source = 'NASA/GSFC/GMAO GEOS-5 Aerosol Group'
       contact = 'arlindo.dasilva@nasa.gov'

       if filename is None:
           filename = '%s/%s.sfc.%d_%02dz.nc4'%(dir,expid,self.nymd,self.nhms/10000)

       # Create the file
       # ---------------
       f = GFIO()
       f.create(filename, vname, nymd, nhms,
                lon=glon, lat=glat, levs=levs, levunits='nm',
                vtitle=vtitle, vunits=vunits,kmvar=kmvar,amiss=MISSING,
                title=title, source=source, contact=contact)

       # Subset AOD at specified channels
       # --------------------------------
       I = []
       for ch in channels:
           i = list(self.channels).index(ch)
           I = I + [i,]
       aod = self.aod[:,I]

       # Fine mode
       # ---------
       try:
           aod_fine = self.aod_fine[:,I]
       except:
           aod_fine = MISSING * np.ones(aod.shape) # will compress like a charm
       
       # The Revised AOD may not exist
       # -------------------------------
       try:
           aod_ = self.aod_[:,I]
       except:
           aod_ = MISSING * np.ones(aod.shape) # will compress like a charm

       # Grid variable and write to file
       # -------------------------------
       f.write('tau', nymd, nhms, 
               binobs3d(self.lon,self.lat,aod,im,jm,MISSING) )
       f.write('tau_', nymd, nhms, 
               binobs3d(self.lon,self.lat,aod_,im,jm,MISSING) )
       f.write('tau_fine', nymd, nhms, 
               binobs3d(self.lon,self.lat,aod_fine,im,jm,MISSING) )
       f.write('count_tau', nymd, nhms,
               binobscnt3d(self.lon,self.lat,aod,im,jm,MISSING) )
       f.write('count_tau_', nymd, nhms,
               binobscnt3d(self.lon,self.lat,aod_,im,jm,MISSING) )
       f.write('cloud', nymd, nhms, 
               binobs2d(self.lon,self.lat,self.cloud,im,jm,MISSING) )
           
#       try:
#           f.close()
#       except:
#           pass

       if Verb >=1:
           print "[w] Wrote file "+filename

#---
    def addVar(self,ga,expr='mag(u10m,v10m)',vname='wind',clmYear=None,tight=True):
        """
        Given a grads object *ga* having the correct MERRA file as default,
        interpolates *var* to obs location and saves it as an attribute
        named *vname*.

        If *tight* is True, domain will be restricted conserve memory. This feature
        has proven somewhat unstable for reasons yet TBD.
        """

        U = MISSING * np.ones(self.nobs)
        if vname == None:
            vname = expr

        # nearest time
        # ------------
        t = _gatime(self.nymd,self.nhms)
        if clmYear != None:
            t = t[:-4] + str(clmYear) # replace year
        ga('set time '+t,Quiet=True)

        # To conserve memory, restrict domain with 1 gridpoint halo
        # ---------------------------------------------------------
        if tight:
            fh = ga.query("file")
            x1, x2  = self.lon.min(),self.lon.max()
            y1, y2  = self.lat.min(),self.lat.max()
            ga('set lon %f %f'%(x1,x2),Quiet=True)
            ga('set lat %f %f'%(y1,y2),Quiet=True)
            qh = ga.query("dims")
            x1, x2 = (qh.xi[0]-1,qh.xi[1]+1)
            y1, y2 = (max(1,qh.yi[0]-1),min(fh.ny,qh.yi[1]+1)) # in [1,ny]
            ga('set x %d %d'%(x1,x2),Quiet=True) # make sure x range is int
            ga('set y %d %d'%(y1,y2),Quiet=True) # make sure y range is int
            expr_ = ga.exp(expr)
        else:
            expr_ = ga.expr(expr)
        u, levs = ga.interp(expr_, self.lon, self.lat )
        U = u.data
        if len(np.shape(U)) == 0:
             U = U * np.ones(1) # so that we can slice it later

        self.__dict__[vname] = U

#---
    def getCoxMunk(self,filename='/nobackup/NNR/Misc/coxmunk_lut.npz',channel=550):
        """
        Returns ocean albedo.
        """
        
        # Get precomputed albedo LUT
        # --------------------------
        lut = NPZ(filename)
        
        # Trimmed wind speed
        # ------------------
        w10m = self.wind.copy()
        w10m[w10m<0] = 0
        w10m[w10m>50.] = 50.

        j = list(lut.channels).index(channel)

        # Interpolate albedo
        # ------------------
        albedo = np.zeros(len(w10m))
        albedo[:] = np.interp(w10m,lut.speed,lut.albedo[:,j])

        self.albedo = albedo

#............................................................................

def granules ( path, algo, sat, syn_time, coll='011', nsyn=8, verbose=False ):
    """
    Returns a list of Vx04 granules for a given algorithm at given synoptic time.
    On input,

    path      ---  mounting point for the MxD04 Level 2 files
    algo      ---  either DT_LAND, DT_OCEAN, DB_LAND or DB_OCEAN
    sat       ---  SNPP
    syn_time  ---  synoptic time (timedate format)

    coll      ---  collection: 011 (optional)
    nsyn      ---  number of synoptic times per day (optional)

    """

    # Get product name
    # -----------------
    Algo = algo.split('_')[0]
    prod = 'AER{}'.format(Algo)

    # Determine synoptic time range
    # -----------------------------
    dt = timedelta(seconds = 12. * 60. * 60. / nsyn)
    t1, t2 = (syn_time-dt,syn_time+dt)

    # Find VIIRS granules in synoptic time range
    # ------------------------------------------
    dt = timedelta(minutes=6)
    t = datetime(t1.year,t1.month,t1.day,t1.hour,0,0)
    Granules = []
    while t < t2:
        if t >= t1:
            doy = t.timetuple()[7]
            basen = "%s/%s/%s/%s/Level2/%04d/%03d/AER%s_L2_VIIRS_%s.A%04d%03d.%02d%02d.%s.*.nc"\
                     %(path,prod,sat,coll,t.year,doy,Algo,sat,t.year,doy,t.hour,t.minute,coll)
            try:
                filen = glob(basen)[0]
                Granules += [filen,]
                if verbose:
                    print " [x] Found "+filen
            except:
                pass
        t += dt

    if len(Granules) == 0:
        print "WARNING: no %s collection %s granules found for"%(algo,coll), syn_time

    return Granules

#--

def print_stats(name,x=None):
    "Prints simple stats"
    from pylab import prctile
    if type(name) is not str:
        x = name
        name = 'mean,stdv,rms,min,25%,median,75%,max: '
    if name == '__header__':
        print ''
        n = (80 - len(x))/2
        print n * ' ' + x
        print n * ' ' + len(x) * '-'
        print ''
        print '   Name       mean      stdv      rms      min     25%    median     75%      max'
        print ' ---------  -------  -------  -------  -------  -------  -------  -------  -------'
    elif name == '__sep__':
        print ' ---------  -------  -------  -------  -------  -------  -------  -------  -------'
    elif name == '__footer__':
        print ' ---------  -------  -------  -------  -------  -------  -------  -------  -------'
        print ''
    else:
        ave = x.mean()
        std = x.std()
        rms = np.sqrt(ave*ave+std*std)
        prc = prctile(x)
        print '%10s  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  '%\
            (name,ave,std,rms,prc[0],prc[1],prc[2],prc[3],prc[4])

#--

def _gatime(nymd,nhms):
        Months = ('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')
        cymd = "%08d"%int(nymd)
        chms = "%06d"%int(nhms)
        t = chms[0:2]+":"+chms[2:4]+"Z"+\
            cymd[6:8]+Months[int(cymd[4:6])-1]+cymd[0:4]
        return t
    
#............................................................................

if __name__ == "__main__":

    syn_time = datetime(2013,10,26,10,0,0)
    Files = granules('/nobackup/VIIRS/','DB','SNPP',syn_time,coll='011')

    db_ocean = Vx04_L2(Files,'DB','OCEAN',syn_time,Verb=1,only_good=True)
    
