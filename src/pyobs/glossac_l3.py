#!/usr/bin/env python
"""
   Implements Python interface to GloSSAC  Level-3 data.
   Data is obtained from: https://asdc.larc.nasa.gov/project/GloSSAC/GloSSAC_2.23
   Presented in monthly, zonal means
    time: 552 values in YYYYMM form - 197901 - 202412
    lat:   32 values                - from -77.5 to 77.5 in 5 degree resolution
    alt:   80 values                - from 0.5 - 40 km in 0.5 km resolution
    wavelengths_glossac: 4          - 386, 452, 525, 1020 nm
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
import matplotlib.dates as mdates
from matplotlib.dates import MonthLocator, WeekdayLocator, DateFormatter
import matplotlib.ticker as ticker

MISSING = -9999.0


#.........................................................................................

class GLOSSAC_L3(object):
   """
    Base class for generic GLOSSAC object.
   """

   def __init__ (self,Path=None,Verbose=0,only_good=True):
     """
       Creates an GLOSSAC object defining the attributes corresponding
       to the SDS's on input.

     """

     self.verb = Verbose

     # Path to data
     if Path == None:
        fname = "/home/pcolarco/geos_aerosols/pcolarco/GloSSAC_V2.23_NC4.nc"

     # Open file
     self.aer = xr.open_mfdataset(fname)

     
     # Variable names
     # --------------
     self.lat  = self.aer["lat"]
     self.time = self.aer["time"]
     self.alt  = self.aer["alt"]
     self.saod = self.aer["Glossac_Aerosol_Optical_Depth"]
     
          
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
   def plotsaod(self, showplot=False):
       lat  = self.lat.values
       saod = np.transpose(np.squeeze(self.saod.values[:,:,3]))
       yyyymm = self.time.values
       dates = []
       for x in yyyymm:
           dates.append(datetime.strptime(x.astype(str).tolist(),"%Y%m"))
       
       saodt = []
       for i in np.arange(0,len(yyyymm)):
           saodt.append(np.sum(np.squeeze(saod[:,i])*np.cos(lat*np.pi/180.))/np.sum(np.cos(lat*np.pi/180.)))

       fig, ax = plt.subplots(1, 1, figsize=(20, 5))
       clevs = np.arange(0,0.03,0.0002)
#       clevs = np.append(clevs,[0.03,0.05,0.1,0.2,0.5])
#       print(clevs)
       im = ax.contourf(dates,lat,saod,clevs,cmap="Spectral_r",extend="both")
#       cbar=plt.colorbar(im, ax=ax, location="bottom", pad=0.01, extend="max", shrink=.6)
#       cbar.ax.tick_params(labelsize=12)
#       plt.xticks(ticks=[2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023])
       plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
#       cbar.set_label(label='Stratospheric AOD 1020 nm',size=12)
#       plt.tight_layout()
       ax.set_xlim([date(2012, 1, 1), date(2026, 2, 1)])
       plt.ylim((-80,80))
       ax.set_ylabel("Latitude")

       ax2 = ax.twinx()
       ax2.set_ylim(0,0.02)
       ax2.plot(dates,saodt, color="black",linewidth=3)
       ax2.set_ylabel("Global Mean Stratospheric AOD 1020 nm")
       if showplot:
           plt.show()
       fig.savefig("glossac.png")

   def plotsaodt(self, showplot=False):
       lat  = self.lat.values
       saod = np.transpose(np.squeeze(self.saod.values[:,:,3]))
       yyyymm = self.time.values
       dates = []
       for x in yyyymm:
           dates.append(datetime.strptime(x.astype(str).tolist(),"%Y%m"))
       saodt = []
       for i in np.arange(0,len(yyyymm)):
           saodt.append(np.sum(np.squeeze(saod[:,i])*np.cos(lat*np.pi/180.))/np.sum(np.cos(lat*np.pi/180.)))
       
       fig, ax = plt.subplots(1, 1, figsize=(20, 5))
       im = ax.plot(dates,saodt)
       plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
       plt.tight_layout()
       ax.set_xlim([date(2012, 1, 1), date(2026, 2, 1)])
       ax.set_ylim([0,0.01])
       if showplot:
           plt.show()
       fig.savefig("glossact.png")

#............................................................................

if __name__ == "__main__":

    glo = GLOSSAC_L3()
    glo.plotsaod()
