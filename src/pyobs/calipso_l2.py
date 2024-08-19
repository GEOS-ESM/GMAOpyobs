#!/bin/env python
"""
   Implements Python interface to CALIPSO Level-2 data.
"""

import os
from pyhdf import SD, HDF
import pyhdf.VS
from glob     import glob
import numpy as np
from   datetime import date, datetime, timedelta
MISSING = -9999.0


ALIAS = dict (
        
                                                           Latitude = 'lat' ,
                                                          Longitude = 'lon' ,
                                                   Profile_UTC_Time = 'Time',
                                                           Pressure = 'plev',
                                       Surface_Elevation_Statistics = 'surfelev',
                                                        Temperature = 'T',
                                  Total_Backscatter_Coefficient_532 = 'tback',
                     Total_Backscatter_Coefficient_Uncertainty_532  = 'tback_err',
                                         Extinction_Coefficient_532 = 'ext',
                             Extinction_Coefficient_Uncertainty_532 = 'ext_err' )


SDS = list(ALIAS.keys())

#.........................................................................................

class CALIPSO_L2(object):
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
            print("WARNING: Empty CALIPSO object created")
            return
     else:
         Path = [Path, ]

     self._readList(Path)
     
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
 
     # Make aliases for compatibility with older code
     # ----------------------------------------------
#    Alias = ALIAS.keys()
     for name in self.Names:
         if self.verb:
             print('shape', name,np.array(self.__dict__[name]).shape)
         if name in SDS:
             self.__dict__[ALIAS[name]] = self.__dict__[name]
     
     # ODS friendly attributes
     # -----------------------
     pass # todo
     
     
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

        for name in self.SDS:
            v = name

#      ------------------------------------------------------------------
#            Temperory:

            if v == 'Profile_UTC_Time':
                sd = f.select(v)
                Time  = sd.get()
                nobs  = len(Time)
              
                nymd  = np.ones(nobs).astype('int')
                nhms  = np.ones(nobs).astype('int')
                self.__dict__[v].append(Time)       # time as on file
                for i in range(nobs):
                    yymmdd = Time[i,1]   # 3 times reported, choose the middle
                    nymd0 = int(Time[i,1])               
                    nd    = Time[i,1] - nymd0
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
                    self.time.append(datetime(year,month,day,hh,mm,ss))
                              
            else:
                if self.verb:
                    print('v', v)
                sd = f.select(v)
              
                data  = sd.get()  # most of parameter : data = (nobs) or (nobs,km) except L2 feature type(nobs,km,4)
                if v == 'Temperature':
                    data = data + 273.15

                # Read attributes.
                attrs       = sd.attributes(full=1)
                fva         = attrs["fillvalue"]
                fillvalue   = fva[0]
                ua          = attrs["units"]
                units       = ua[0]
                vra         = attrs["valid_range"]
                valid_range = vra[0].split('...')

                # Filter fill value and valid range. See Table 66 (p. 116) from [1]
                data = np.ma.array(data)
                data.mask = data == fillvalue
  
                # Apply the valid_range attribute.
                # sometimes the upper limit is not defined, so just try
                try:
                    invalid = (data < float(valid_range[0])) | (data > float(valid_range[1]))
                except:
                    invalid = (data < float(valid_range[0]))    
                data.mask[invalid] = True
                sd.endaccess()


                self.__dict__[v].append(data)


        # read altitude from metadata
        # only need to do this once:
        if not hasattr(self,'alt'):
            hdf = HDF.HDF(filename)
            vs  = hdf.vstart()
            meta = vs.attach("metadata")
            lidar_alt_field = meta.field('Lidar_Data_Altitudes')
            record_index = 0
            all_data = meta.read(meta._nrecs)[record_index]
            self.alt = np.array(all_data[lidar_alt_field._idx])

            meta.detach()
            vs.end()
            hdf.close()


        # clean up
        f.end()
          
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
           
#....................................................................

def _timefilter ( t, t1, t2, a, I_bad ):
    filler = MISSING * np.ones(a.shape[1:])
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
def orbits (path, syn_time, nsyn=8, Verbose=0 ):
    """
    Returns a list of CALIPSO orbits for a given product at given synoptic time.
    On input,

    path      ---  mounting point for the CALIPSO Level 2 files
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
       
        dirn = "%s/"%(path)        
        Files += glob("%s/CAL_LID_L2_05kmCPro-Prov-V3-01.*.hdf"%(dirn))
        
    Orbits = []
    for f in Files:
        dirn, filen = os.path.split(f)
        
        tokens = filen.split('-')
        beg_yy = int(tokens[3].split('.')[1])
        beg_mm = int(tokens[4])
        beg_dd = int(tokens[5].split('T')[0])
        beg_h = int(tokens[5].split('T')[1])
        beg_m = int(tokens[6])
        t_beg = datetime(beg_yy,beg_mm,beg_dd,beg_h,beg_m,0)
        t_end = t_beg + timedelta(minutes=90)
        # t_end = datetime(end_yy,end_mm,end_dd,end_h,end_m,0)
        
        if (t_beg>=t1 and t_beg<t2) or (t_end>=t1 and t_end<t2):
            print("[x] ", t_beg, '|', t_end)
            Orbits += [f,]
            
            if Verbose:
                print("[] ", f)

    return Orbits
#............................................................................

if __name__ == "__main__":

#    syn_time = datetime(2008,6,30,0,0,0)
    syn_time = datetime(2011,6,1,9,0,0)
#    Files = orbits('/nobackup/CALIPSO/Level1.5',syn_time,Verbose=1)
    
#    essai  = CALIPSO_L1p5(Files,Verbose=1)
    
#    essai.writeg(syn_time=syn_time,nsyn=8,filename=None,\
#    dir='/nobackup/2/vbuchard/',expid='calipso_lev2',refine=8,res=None,Verb=1)
