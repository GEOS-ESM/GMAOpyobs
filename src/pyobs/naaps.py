#!/bin/env python
"""
   Implements Python interface to NRL NAAPS files
"""

import os
import sys

from types    import *
from pyhdf import SD
from glob     import glob
from   numpy    import ones, concatenate, array,linspace,arange, transpose
from   datetime import date, datetime, timedelta

from .config import strTemplate

MISSING = -9999.99

ALIAS = dict (latitude = 'lat' ,
              longitude = 'lon' ,
              elevation = 'zs' ,
                   time = 'Time')
        
ALIAS['532_attenuated_backscatter'] = 'taback'
ALIAS['532_attenuated_backscatter_error'] = 'taback_err'
ALIAS['532_attenuated_molecular_backscatter'] = 'mol_aback' 
                                                     
                     

SDS = list(ALIAS.keys())

#.........................................................................................

class NAAPS(object):
   """
    Base class for NAAPS object.
   """

   def __init__ (self,Path,keep=None,Verbose=0,only_good=True):
     """
       Creates an NAAPS object defining the attributes corresponding
       to the SDS's on input.
        The optional parameter *keep* is used to specify the number of scan
        lines (from the left of the swath) to keep. This is needed for
        coping with the row anomaly problem.

     """

     #  Initially are lists of numpy arrays for each granule
     # ----------------------------------------------------
     self.verb = Verbose
     self.keep = keep
     self.SDS = SDS

     # Variable names
     # --------------
     self.Names = []

     for name in SDS:
             
             self.Names.append(name)
     self.Names += ['nymd','nhms']

     # Create empty lists for SDS to be read from orbit file;
     #  each element of the list contains data for one orbit
     # ------------------------------------------------------

     for name in self.Names:
         self.__dict__[name] = []
     self.time_ = [] # to hold datetime objects
     

     # Read each orbit, appending them to the list
     # -------------------------------------------
     
     if type(Path) is ListType:
        if len(Path) == 0:
            self.nobs = 0
            print("WARNING: Empty NAAPS object created")
            return
     else:
         Path = [Path, ]

     self._readList(Path)
     

     # Make each attribute a single numpy array
     # ----------------------------------------
     for name in self.Names:
#            print 'name',name, 'donnees',self.__dict__[name]
            try:
                self.__dict__[name] = concatenate((self.__dict__[name]))
                
            except:
                print("Failed concatenating "+name)
     
    
 
     # Make aliases for compatibility with older code
     # ----------------------------------------------
#    Alias = ALIAS.keys()
     for name in self.Names:
         
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
        """Reads one CALIPSO orbit with Level 1.5 data."""

        # Reference time
        # --------------
        REF_DATE = datetime(1993,1,1,0,0,0)

        # Open the CALIPSO file and loop over the datasets,
        # extracting GEOLOCATION and Data fields
        # ----------------------------------------------

        if self.verb:
            print("[] working on <%s>"%filename)

        f = SD.SD(filename)

#       for group in self.SDS.keys():
        for name in self.SDS:         
          v = name
          print('v', v)

          if v == 'time':
              sd = f.select(v)
              Time  = sd.get()
              nobs  = len(Time)
              
              nymd  = ones(nobs).astype('int')
              nhms  = ones(nobs).astype('int')
              self.__dict__[v].append(Time)       # time as on file
              for i in range(nobs):
                yymmdd = Time[i]
                nymd0 = int(Time[i])               
                nd    = Time[i] - nymd0
                nd0   = nd * 24.0
                hh    = int(nd0)
                nd1   = nd0 - hh
                nd2   = nd1 * 60
                mm    = int(nd2)
                nd3   = nd2 - mm
                nd4   = nd3 * 60
                ss    = int(nd4)
                                
                nymd[i]  = 20000000 + nymd0
                nhms[i]  = ((hh * 100) + mm) * 100 + ss
                
                self.nymd.append(nymd)
                self.nhms.append(nhms)
                
                year = int(nymd[i]/10000)
                month = int((nymd[i] - 10000*year)/100)
                day = nymd[i] - (year*10000 + month * 100) 
                self.time_.append(datetime(year,month,day,hh,mm,ss))
                              
          else:

              sd = f.select(v)
              data  = sd.get()  # most of parameter : data = (nobs) or (nobs,km) except L2 feature type(nobs,km,4)
              data = transpose(data)
              print('data', data.shape)
              if self.keep != None:
                    self.__dict__[v].append(data[0:self.keep,:])
              else:
                    self.__dict__[v].append(data)
          
#---
   def writeg(self,g5,syn_time,nsyn=8,g5_h=None,g5_ab=None,filename=None,dir='.',expid='NAAPS',Verb=1):
       """
        Writes gridded CALIPSO measurements to file (same grid as GEOS-5 file).
        Verb -- Verbose level:
                 0 - really quiet (default)
                 1 - Warns if invalid file is found
                 2 - Prints out non-zero number of fires in each file.

       """

       from gfio     import GFIO
       from binObs_  import binobs3dh

       # Determine synoptic time range
       # -----------------------------
       dt = timedelta(seconds = 12. * 60. * 60. / nsyn)
       t1, t2 = (syn_time-dt,syn_time+dt)
      
#      Lat lon grid from GEOS-5 file
#      ------------
       im = 360
       jm = 181
       print('im,jm', im, jm)
       glon = linspace(-180.,180.,im,endpoint=False)
       glat = linspace(-90.,90.,jm)
       print('glon', glon, glat)       

       dLon = 360. / im
       dLat = 180. / ( jm - 1.)
       print('dlon', dLon, dLat)
       
       nymd = 10000 * syn_time.year + 100 * syn_time.month  + syn_time.day
       nhms = 10000 * syn_time.hour + 100 * syn_time.minute + syn_time.second

       print('nymd=',nymd, 'nhms=',nhms) 
       na_height = arange(0,8100,400) # height above sea level for NAAPS 100mfor night 400m forday
       
       print('na_height shape', na_height.shape, g5_h.shape)  

       g5_height = g5_h
       km = g5_height.shape[0] # because it is at the edge
       print('km', km, g5_height.shape, g5_height[:,0])    
       nobs = self.lon.shape

       vtitle = [ 'taback',
                  'taback_err',
                  'mol_aback',
                  'height' ]

       vname  = ['taback','taback_err',  'mol_aback']
       vunits = [ 'km-1 sr-1','km-1 sr-1', 'km-1 sr-1' ]
       kmvar  = [km,      km,    km]

       title = 'Gridded NAAPS attenuated backscatter coeff lev Geos5'
       source = 'NASA/GSFC/GMAO GEOS-5 Aerosol Group'
       contact = 'Virginie'

       if filename is None:
           filename = '%s/%s.day.calipso_l3a.%d_%02dz.nc4'%(dir,expid,nymd,nhms/10000)
       

       # QA filtering
       # ------------
       I_bad = ones(self.taback.shape) # bad data
       I_bad = False
       
       # Time filter of data
       # -------------------
       lon = self.lon
       lat = self.lat

       taback = _timefilter(self.time_,t1,t2,self.taback,I_bad)
       taback_err = _timefilter(self.time_,t1,t2,self.taback_err,I_bad)
       mol_aback = _timefilter(self.time_,t1,t2,self.mol_aback,I_bad)
#       height = _timefilter(self.time_,t1,t2,na_height,I_bad)

       print('taback', taback.shape)

#      Create the file
#      ---------------
       f = GFIO()
       glevs=arange(km)
       f.create(filename, vname, nymd, nhms,
                lon=glon, lat=glat, levs=glevs, levunits='m',
                vtitle=vtitle, vunits=vunits,kmvar=kmvar,amiss=MISSING,
                title=title, source=source, contact=contact)
         
#       gObs=binobs3dh(lon[13:14],lat[13:14],taback[13:14,:],na_height,g5_height[:,13:14],im,jm,MISSING)
       print('test', lon[10:11],lat[10:11],taback[10:11,:],na_height,g5_height[:,10:11])
       gObs=binobs3dh(lon[10:11],lat[10:11],taback[10:11,:],na_height,g5_height[:,10:11],im,jm,MISSING)
       
       print('gobs', gObs[357:358,101:102,:])

#      Grid variable and write to file
#      -------------------------------
       f.write('taback', nymd, nhms, binobs3dh(lon,lat,taback,na_height,g5_height,im,jm,MISSING) )
       f.write('taback_err', nymd, nhms, binobs3dh(lon,lat,taback_err,na_height,g5_height,im,jm,MISSING) )
       f.write('mol_aback', nymd, nhms, binobs3dh(lon,lat,mol_aback,na_height,g5_height,im,jm,MISSING) ) 
#       f.write('height', nymd, nhms, g5_height)

       if Verb >=1:
           print("[w] Wrote file "+filename)
           
#....................................................................

def _timefilter ( t, t1, t2, a, I_bad ):
    filler = MISSING * ones(a.shape[1:])
    b = a.copy()
    for i in range(len(t)):
        if (t[i]<t1) or (t[i]>=t2):
            b[i] = filler
    if len(b.shape) == 3:
        b[I_bad,:] = MISSING
    elif len(b.shape) == 2:
        b[I_bad] = MISSING
    else:
        raise IndexError("Invalid rank=%d for time filtering"%len(b.shape))
    return b
#---
def orbits (path, syn_time, nsyn=8, period='night', Verbose=0 ):
    """
    Returns a list of CALIPSO orbits for a given product at given synoptic time.
    On input,

    path      ---  mounting point for the CALIPSO Level 1.5 files
    syn_time  ---  synoptic time (timedate format)

    nsyn      ---  number of synoptic times per day (optional)

    """

    # Determine synoptic time range
    # -----------------------------
    dt = timedelta(seconds = 12. * 60. * 60. / nsyn)
    t1, t2 = (syn_time-dt,syn_time+dt)

    print("[*] ", t1,"|", t2)

    today = syn_time
    yesterday = today - timedelta(hours=24)

    Files = []
    for t in (yesterday,today):

        yy, mm, dd = (t.year,t.month,t.day)        
       
        dirn = "%s/%02d/%s"%(path,mm,period)        
        Files += glob("%s/naaps_caliop_assim_*.cdf"%(dirn))
#    print 'Files', dirn, Files     
    Orbits = []
    for f in Files:
        dirn, filen = os.path.split(f)
        
        tokens = filen.split('_')
        beg_yy = int(tokens[3][0:4])
        beg_mm = int(tokens[3][4:6])
        beg_dd = int(tokens[3][6:8])
        beg_h = int(tokens[3][8:10])
        beg_m = int(tokens[3][10:12])
        t_beg = datetime(beg_yy,beg_mm,beg_dd,beg_h,beg_m,0)
        t_end = t_beg + timedelta(minutes=90)
        # t_end = datetime(end_yy,end_mm,end_dd,end_h,end_m,0)
#        print 'year', beg_yy, 'month', beg_mm, 'day', beg_dd, 'hour', beg_h, 'min', beg_m
        if (t_beg>=t1 and t_beg<t2) or (t_end>=t1 and t_end<t2):
            print("[x] ", t_beg, '|', t_end)
            Orbits += [f,]
            
            if Verbose:
                print("[] ", f)

    return Orbits
#............................................................................

if __name__ == "__main__":

#    syn_time = datetime(2008,6,30,0,0,0)

    # Time interval snd time step
    # ---------------------------
    t_beg = datetime(2007,4,1,0)
    t_end = datetime(2007,4,1,21)
    dt = timedelta(seconds=3*60*60) # 3-hourly
    t = t_beg - dt
    while t < t_end:

        t += dt       
        syn_time = t
        Files = orbits('/nobackup/2/vbuchard/CALIPSO_L15/NAAPS/',syn_time,period='day',Verbose=1)
        print('files',Files)
#def hold():
    # NAAPS files
        naap  = NAAPS(Files,Verbose=1)
    # GEOS-5 file
        g_template = "/nobackup/2/vbuchard/CALIPSO_L15/GEOS-5/aback_63lay/Y%y4/M%m2/dR_MERRA-AA-r2_ext532nm_Nv_63layers.%y4%m2%d2_%h200z.nc4"
        g_fn = strTemplate(g_template,dtime=syn_time)
  
        lon=naap.lon
        lat=naap.lat
        g = GFIO(g_fn)
        g5_height = g.interp('h',lon,lat)
        g5_aback = g.read('taback')

        naap.writeg(g,syn_time,nsyn=8,g5_h=g5_height,g5_ab=g5_aback,filename=None,dir='/nobackup/2/vbuchard/CALIPSO_L15/',expid='NAAPS',Verb=1)



