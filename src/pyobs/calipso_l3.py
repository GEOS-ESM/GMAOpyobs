#!/usr/bin/env python
"""
   Implements Python interface to CALIPSO Level-3 data.
   Data is obtained from: https://asdc.larc.nasa.gov/project/CALIPSO/CAL_LID_L3_Tropospheric_APro_AllSky-Standard-V4-20_V4-20
   Presented in monthly means at 2 degree latitude (85 nodes) x
                                 5 degree longitude (72 nodes) x
                                 60 m vertical (208 nodes)
   Code is adapted from calipso_l2.py
   
"""

import os
from glob     import glob
import numpy as np
from   datetime import date, datetime, timedelta
import sys
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

MISSING = -9999.0

months = ['January','February','March','April','May','June',
          'July','August','September','October','November','December']

ALIAS = dict (
        
                                                  Latitude_Midpoint = 'lat' ,
                                                 Longitude_Midpoint = 'lon' ,
                                                  Altitude_Midpoint = 'alt',
                                                      Pressure_Mean = 'pres',
                                           Surface_Elevation_Median = 'elev',
                                                   Temperature_Mean = 't',
                                    Extinction_Coefficient_532_Mean = 'ext532',
                                                             AOD532 = 'aod532' )


SDS = list(ALIAS.keys())

#.........................................................................................
def getfiles(yy,mm,daynight="N",PATH="/gpfsm/dnb04/projects/p22/aerosol/data/CALIPSO/Level3/LID_L3_Tropospheric_APro_CloudFree-Standard-V5-00/"):
   filen = []
   for m in mm:
       filen_ = f"{PATH}/{yy:04d}/CAL_LID_L3_Tropospheric_APro_CloudFree-Standard-V5-00.{yy:04d}-{m:02d}{daynight}.hdf"
       if not os.path.exists(filen_):
           print(f"File not found: {filen_}. Skipping...")
           continue
       filen.append(filen_)
   return filen
   

class CALIPSO_L3(object):
   """
    Base class for generic CALIPSO object.
   """

   def __init__ (self,Path,Verbose=0,only_good=True):
     """
       Creates an CALIPSO object defining the attributes corresponding
       to the SDS's on input.

     """

     #  Initially are lists of numpy arrays for each granule
     # ----------------------------------------------------
     self.verb = Verbose
     self.SDS = SDS

     # Variable names
     # --------------
     self.Names = []

     for name in SDS:
             
             self.Names.append(name)
     self.Names += ['nymd','nhms']

     # Create empty lists for SDS to be read from file;
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
            print("WARNING: Empty CALIPSO object created")
            return
     else:
         Path = [Path, ]

     self._readList(Path)
     
     """
     # Make each attribute a single numpy array
     # ----------------------------------------
     for name in self.Names:
            try:
                self.__dict__[name] = np.ma.concatenate(self.__dict__[name])
                
            except:
                print("Failed concatenating "+name)

     self.time = np.array(self.time)

     # Determine index of "good" observations
     # --------------------------------------
     pass # to do
     # Keep only "good" observations
     # -----------------------------
     if only_good:
         pass
     """
 
     # Make aliases for compatibility with older code
     # ----------------------------------------------
#    Alias = ALIAS.keys()
     for name in self.Names:
         if self.verb:
             print('shape', name,np.array(self.__dict__[name]).shape)
         if name in SDS:
             self.__dict__[ALIAS[name]] = self.__dict__[name]
     
    
     
#---
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


#---
   def _readDir(self,dir):
        """Recursively, look for files in directory."""
        for item in os.listdir(dir):
            path = dir + os.sep + item
            if os.path.isdir(path):      self._readDir(path)
            elif os.path.isfile(path):   self._readOrbit(path)
            else:
                print("%s is not a valid file or directory, ignoring it"%item)

#---
   def _readOrbit(self,filename):
        """Reads one CALIPSO Level3 file."""

        # Reference time
        # --------------
        REF_DATE = datetime(1993,1,1,0,0,0)

        # Open the CALIPSO file and loop over the datasets,
        # extracting GEOLOCATION and Data fields
        # ----------------------------------------------

        if self.verb:
            print("[] working on <%s>"%filename)

        obs_ds = xr.open_dataset(filename,engine="netcdf4")

        for name in self.SDS:
            v = name
            
            if self.verb:
                print('v', v)
            if v == "AOD532":
                data = np.squeeze(obs_ds["Extinction_Coefficient_532_Mean"].values)
                data = np.nansum(data*0.06,axis=2)  # 60 m = 0.06 km vertical
            else:
                data = np.squeeze(obs_ds[v].values)
            self.__dict__[v].append(data)



          
#---
   def writeg(self,aer,syn_time,nsyn=8,filename=None,dir='.',expid='calipso_lev2',Verb=1):
       """
        Writes gridded CALIPSO measurements to file (same grid as GEOS-5 file).
        Verb -- Verbose level:
                 0 - really quiet (default)
                 1 - Warns if invalid file is found
                 2 - Prints out non-zero number of fires in each file.

       """
       from gfio     import GFIO
       from binObs_  import binobs3d, binobs3dp

       # Determine synoptic time range
       # -----------------------------
       dt = timedelta(seconds = 12. * 60. * 60. / nsyn)
       t1, t2 = (syn_time-dt,syn_time+dt)
      
#      Lat lon grid from GEOS-5 file
#      ------------
       im = aer.im
       jm = aer.jm
      
       glon = np.linspace(-180.,180.,im,endpoint=False)
       glat = np.linspace(-90.,90.,jm)
      
       nymd = 10000 * syn_time.year + 100 * syn_time.month  + syn_time.day
       nhms = 10000 * syn_time.hour + 100 * syn_time.minute + syn_time.second

       print('nymd=',nymd, 'nhms=',nhms) 
       km = aer.km             # vertical levels
       ptop = 1.               # ~ 1 Pa at the top

#      GEOS-5 edge pressure [Pa]
#      ------------------------- 
       pe = np.ones((im,jm,km+1)) 
       pe[:,:,0] = ptop
       for k in range(aer.km):
           pe[:,:,k+1] = pe[:,:,k] + aer.read('delp',nymd=nymd,nhms=nhms)[:,:,k]   
       pe = pe * 0.01         # in hPa

#      GEOS-5 mid-level pressure [Pa]
#      -----------------------------
       plev = np.ones((im,jm,km)) # mid-level pressure [Pa]
       plev[:,:,0] = ptop + aer.read('delp',nymd=nymd,nhms=nhms)[:,:,0]/2. 
       for k in range(aer.km-1):  
           plev[:,:,k+1] = plev[:,:,k] + aer.read('delp',nymd=nymd,nhms=nhms)[:,:,k]/2.\
           + aer.read('delp',nymd=nymd,nhms=nhms)[:,:,k+1]/2.      
       plev = plev * 0.01
       
       vtitle = [ 'tback',
                  'tback_err',
                  'extinction', 
                  'ext_err',
#                  'mol_aback',
                  'pressure' ]

       vname  = ['tback','tback_err', 'ext', 'ext_err', 'pressure' ]
       vunits = [ 'km-1 sr-1','km-1 sr-1',  'km-1',   'km-1', 'hPa' ]
       kmvar  = [km,      km,    km,    km,  km ]

       title = 'Gridded CALIPSO Level 2 version 3.01 data'
       source = 'NASA/GSFC/GMAO GEOS-5 Aerosol Group'
       contact = 'arlindo.dasilva@nasa.gov'

       if filename is None:
           filename = '%s/%s.obs_l3a.%d_%02dz.nc4'%(dir,expid,nymd,nhms/10000)
       
#      Create the file
#      ---------------
       f = GFIO()
       glevs=np.arange(km)
       f.create(filename, vname, nymd, nhms,
                lon=glon, lat=glat, levs=glevs, levunits='hPa',
                vtitle=vtitle, vunits=vunits,kmvar=kmvar,amiss=MISSING,
                title=title, source=source, contact=contact)

       # QA filtering
       # ------------
       I_bad = np.ones(self.tback.shape) # bad data
       I_bad = False
       
       # Time filter of data
       # -------------------
#       lon = self.lon.ravel()
#       lat = self.lat.ravel()

       lon = self.lon[:,1]  # choose the middle pulse
       lat = self.lat[:,1]


       tback = _timefilter(self.time,t1,t2,self.tback,I_bad)
       tback_err = _timefilter(self.time,t1,t2,self.tback_err,I_bad)
       ext = _timefilter(self.time,t1,t2,self.ext,I_bad)
       ext_err = _timefilter(self.time,t1,t2,self.ext_err,I_bad)
       pressure = _timefilter(self.time,t1,t2,self.plev,I_bad)
       
       gObs=binobs3dp(lon,lat,pressure,tback,pe,MISSING)
       
#      Grid variable and write to file
#      -------------------------------
       f.write('tback', nymd, nhms, binobs3dp(lon,lat,tback,pressure,pe,MISSING) )
       f.write('tback_err', nymd, nhms, binobs3dp(lon,lat,tback_err,pressure,pe,MISSING) )
       f.write('ext',    nymd, nhms, binobs3dp(lon,lat,ext,pressure,pe,MISSING) )
       f.write('ext_err', nymd, nhms, binobs3dp(lon,lat,ext_err,pressure,pe,MISSING) )  
       f.write('pressure', nymd, nhms, plev)

       if Verb >=1:
           print("[w] Wrote file "+filename)

#--
   def plotaod(self,yy,mm,varn="AOD532",extent=[-180,180,-90,90],box=None,prange=[0,1],
               var2=None, showplot=False):
       lon = self.lon[0]
       lat = self.lat[0]
       if (len(mm) > 1):
           fig, axs = plt.subplots(3,4,figsize=(24,14),layout='constrained',
                                   subplot_kw={'projection': ccrs.PlateCarree()})
           i = 0
           for ax in axs.flat:
              var = self.__dict__[varn][i]
              title = f"{yy:04d}-{mm[i]:02d}"
              im = self.plotone(lon,lat,var,prange=prange,extent=extent,box=box,title=title,axes=ax)
              i += 1
              if(i == len(mm)):
                 break
           cb = fig.colorbar(im, ax=axs, location='bottom',shrink=.6)
           cb.ax.tick_params(labelsize=16)
           fig.suptitle("AOD 532 nm", size=30)
           plt.savefig(f"CALIPSO_L3.{varn}.{yy:04d}annual.png")
       else:
           fig, axs = plt.subplots(1, 1, figsize=(8, 6),
                                    subplot_kw={'projection': ccrs.PlateCarree()})
           var = self.__dict__[varn][0]
           title = f"{yy:04d}-{mm[0]:02d}"
           im = self.plotone(lon,lat,var,prange=prange,extent=extent,box=box,title=title,axes=axs)
           cb = fig.colorbar(im, ax=axs, location='bottom',shrink=.6)
           cb.ax.tick_params(labelsize=16)
           fig.suptitle("AOD 532 nm", size=30)
           plt.savefig(f"CALIPSO_L3.{varn}.{yy:04d}{mm[0]:02d}.png")
           
       if showplot:
           plt.show()

#--
   def plotone(self,lon,lat,var,prange=[0,1],title=None,axes=None,
               extent=[-180,180,-90,90],box=None,
               var2=None):
       """
       Simple plot of single input map (could be model or satellite)
       """
#       fig, axes = plt.subplots(1, 1, figsize=(8, 6), subplot_kw={'projection': ccrs.PlateCarree()})
#       print(type(axes))

       # Common settings
       vmin_aod = prange[0]  # Min value for AOD
       vmax_aod = prange[1]  # Max value for AOD
       cmap_aod = 'Spectral_r'  # Colormap for AOD
       cmap_cnt = 'grey'

       im = axes.pcolormesh(lon,lat,var, cmap=cmap_aod, vmin=vmin_aod, vmax=vmax_aod,
                            transform=ccrs.PlateCarree())
#       axes.set_title(baselinen)
       axes.add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=2)
       axes.set_extent(extent)
       axes.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=1)
       axes.gridlines(draw_labels=True)
       axes.set_title(title,size=18)
#    if var2[0] != None:
#    im2 = axes.contour(lon, lat, var2, [40], colors='red',linewidths=1,transform=ccrs.PlateCarree())
#       im3 = axes.pcolormesh(lon, lat, var2, cmap=cmap_cnt, vmin=0.01,vmax=.1,transform=ccrs.PlateCarree(),alpha=.2)
       if box != None:
           from matplotlib.patches import Rectangle
           ll = (box[0],box[2])
           wd = box[1]-box[0]
           ht = box[3]-box[2]
           rect = Rectangle(ll, width=wd, height=ht, edgecolor='k',facecolor='none')
           axes.add_patch(rect)

       return im

#............................................................................

if __name__ == "__main__":

    yy = 2012
    mm = [1,2,3,4,5,6,7,8,9,10,11,12]
    Files = getfiles(yy,mm)
    cal = CALIPSO_L3(Files)
    cal.plotaod(yy,mm,extent=[-150,-60,0,55])
