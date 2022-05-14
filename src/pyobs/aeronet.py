"""
   Classes for reading CSV files with AERONET Level 2 files.
"""

import os
import sys
from   types    import *
from   numpy    import loadtxt, ones, zeros, savez, pi, log, concatenate, \
                       arange, savez, shape, array, linspace
from   datetime import datetime, timedelta
from   dateutil.parser import parse as isoparse
from   glob     import glob

from .npz  import NPZ

MISSING = -999.

BAD, MARGINAL, GOOD, BEST = list(range(4))

VARS = ( 'AERONET_Site',
         'Longitude',
         'Latitude',
         'Elevation',
         'Date',
         'Time', 
         'AOT_1640' , 
         'AOT_1020', 
         'AOT_870', 
         'AOT_675', 
         'AOT_667',
         'AOT_555', 
         'AOT_551', 
         'AOT_532', 
         'AOT_531', 
         'AOT_500',
         'AOT_490'  ,
         'AOT_443' , 
         'AOT_440', 
         'AOT_412', 
         'AOT_380', 
         'AOT_340',
         'Water',
         'Solar_Zenith_Angle',
         )

 
ALIAS = dict (          Longitude = 'lon',
                        Latitude = 'lat',
                      Julian_Day = 'doy',
                    AERONET_Site = 'Location', 
              )

KX = 323 
KT = dict ( AOD = 45, )

#---- 
class AERONET_L2(object):
    """Base class for AERONET Level 2 data."""

    def __init__ (self,Path,Vars=VARS,Verbose=False):
        """
        Base class for generic AERONET dataset.
        """

        self.verb = Verbose
        self.columns = None
        self.nobs = 0
        self.Vars = Vars

        # Past is string or list
        # ----------------------
        if type(Path) is ListType:
            if len(Path) == 0:
                print("WARNING: Empty AERONET object created")
                return
        else:
            if Path[-4:] == '.npz': # Special handling of npz files
                npz = NPZ(Path)
                self.__dict__ = npz.__dict__
                self.nobs = self.lon.size
                return
            Path = [Path, ]
            
        # Create empty lists for SDS to be read from file
        # -----------------------------------------------
        for name in self.Vars:
            self.__dict__[name] = []

        # Read each granule, appending them to the list
        # ---------------------------------------------
        self._readList(Path)

        # Make each attribute a single numpy array
        # ----------------------------------------
        for var in self.Vars:
            try:
                self.__dict__[var] = concatenate(self.__dict__[var])
            except:
                print("Failed concatenating "+var)

        # Make aliases
        # ------------
        Alias = list(ALIAS.keys())
        for var in self.Vars:
            if var in Alias:
                self.__dict__[ALIAS[var]] = self.__dict__[var] 

        # Interpolate AOT to 550 nm if needed
        # -----------------------------------
        aot_550 = MISSING * ones(len(self.lon))

        aot_550a = aodInterpAngs(550.,self.AOT_532,self.AOT_667,532.,667.)
        aot_550b = aodInterpAngs(550.,self.AOT_532,self.AOT_675,532.,675.)
        aot_550c = aodInterpAngs(550.,self.AOT_500,self.AOT_667,500.,667.)
        aot_550d = aodInterpAngs(550.,self.AOT_500,self.AOT_675,500.,675.)
        
        aot_550 = _updAOT(aot_550,aot_550d)
        aot_550 = _updAOT(aot_550,aot_550c)
        aot_550 = _updAOT(aot_550,aot_550b)
        aot_550 = _updAOT(aot_550,aot_550a)

        aot_550 = _updAOT(aot_550,self.AOT_555) # close enough
        aot_550 = _updAOT(aot_550,self.AOT_551) # close enough

        self.AOT_550 = aot_550[:] # update undefs with interpolated values

        # Interpolate AOT to 670 nm if needed
        # -----------------------------------
        aot_670 = MISSING * ones(len(self.lon))

        aot_670a = aodInterpAngs(670.,self.AOT_667,self.AOT_675,667.,675.)
        
        aot_670 = _updAOT(aot_670,aot_670a)

        aot_670 = _updAOT(aot_670,self.AOT_667) # close enough
        aot_670 = _updAOT(aot_670,self.AOT_675) # close enough

        self.AOT_670 = aot_670[:] # update undefs with interpolated values        

        # Interpolate AOT to 470 nm if needed
        # -----------------------------------
        aot_470 = MISSING * ones(len(self.lon))

        aot_470a = aodInterpAngs(470.,self.AOT_443,self.AOT_490,443.,490.)
        aot_470b = aodInterpAngs(470.,self.AOT_440,self.AOT_490,440.,490.)
        aot_470c = aodInterpAngs(470.,self.AOT_443,self.AOT_500,443.,500.)
        aot_470d = aodInterpAngs(470.,self.AOT_440,self.AOT_500,440.,500.)
        
        aot_470 = _updAOT(aot_470,aot_470d)
        aot_470 = _updAOT(aot_470,aot_470c)
        aot_470 = _updAOT(aot_470,aot_470b)
        aot_470 = _updAOT(aot_470,aot_470a)

        self.AOT_470 = aot_470[:] # update undefs with interpolated values        
        
        # Create timedate
        # ---------------
        self.tyme = array([ isoparse('-'.join(d.split(':')[-1::-1])+'T'+t) for d,t in zip(self.Date,self.Time) ])

        # Unique list of locations
        # ------------------------
        Locations = {}
        for loc in self.Location:
            Locations[loc] = 1
        self.Stations = list(Locations.keys())

        # By default all is good if coordinates are ok
        # --------------------------------------------
        self.iValid = (abs(self.Longitude)<=180.) & \
                      (abs(self.Latitude)<=90.)

        # clean up
        # --------
        del self.formats, self.converters, self.iVars, self.columns, self.skip

#---
    def _readList(self,List):
        """
        Recursively, look for files in list; list items can
        be files or directories.
        """
        for item in List:
            if os.path.isdir(item):      self._readDir(item)
            elif os.path.isfile(item):   self._readGranule(item)
            else:
                print("%s is not a valid file or directory, ignoring it"%item)
#---
    def _readDir(self,dir):
        """Recursively, look for files in directory."""
        for item in os.listdir(dir):
            path = dir + os.sep + item
            if os.path.isdir(path):      self._readDir(path)
            elif os.path.isfile(path):   self._readGranule(path)
            else:
                print("%s is not a valid file or directory, ignoring it"%item)

#---
    def _readGranule(self,filename):
        """Reads one AERONET granule."""

        # Locate Header
        # -------------
        if self.columns == None: # once per file
            i = 0
            for line in open(filename).readlines():
                if 'AERONET_Site' in line:
                    self.skip = i + 1 # number of rows to skip for data
                    columns = line[0:-1].split(',') # remove \n at end of line
                    self.columns = [] # sanitize variable names
                    for c in columns:                    
                        self.columns += [c.replace('%','').replace('-','_').replace(' ','')\
                                          .replace('AOD','AOT').replace('nm','')\
                                          .replace('_Name','').replace('Precipitable_','')\
                                          .split('(')[0],]
                    break
                i += 1

            if self.columns == None:
                raise ValueError("Cannot find Column header")

            # Read relevant columns from AERONET granule
            # ----------------------------------------
            self.iVars = ()
            self.formats = ()
            self.converters = {}
            for name in self.Vars:
                try:
                    i = self.columns.index(name)
                except:
                    raise ValueError("cannot find <%s> in file <%s>"%(name,filename))
                self.iVars += (i,)
                if name=='Date':
                    self.formats += ('S10',)
                elif name=='Time':
                    self.formats += ('S8',)
                elif name=='AERONET_Site':
                    self.formats += ('S20',)
                else:
                    self.converters[i] = _convert2Float
                    self.formats += ('f4',)
                    
#       Read the data
#       -------------
        data = loadtxt(filename, delimiter=',',
                       dtype={'names':self.Vars,'formats':self.formats},
                       converters = self.converters,
                       skiprows=self.skip, usecols=self.iVars)
        N = len(data)

        self.nobs += N

        # Save data columns as attributes, with nicer names
        # -------------------------------------------------
        for i in range(len(self.Vars)):
            if self.formats[i]=='f4':
                v = ones(N)
                for j in range(N):
                    v[j] = data[j][i]
            else:
                v = []
                for j in range(N):
                    v.append(data[j][i])

            self.__dict__[self.Vars[i]] += [v,]

#---
    def writeNPZ(self,npzFile,I=None):
        """
        Writes out a NPZ file with the relevant variables.
        """
        Vars = dict()
        Nicknames = list(ALIAS.values())
        for name in self.__dict__:
            if name in Nicknames:
                continue # alias do not get reduced
            q = self.__dict__[name]
            if type(q) is type(self.lon):
                if len(q) == self.nobs:
                    if I is None:
                        Vars[name] = self.__dict__[name]    # all observations
                    else:
                        Vars[name] = self.__dict__[name][I] # reduced set
        savez(npzFile,**Vars)
      
#---
    def writeODS(self, syn_tyme, filename=None, dir='.', expid='aeronet', nsyn=8):
        """
        Writes the un-gridded AERONET object to an ODS file.
        """
        from pyods import ODS
        
        nymd = syn_tyme.year*10000 + syn_tyme.month*100  + syn_tyme.day
        nhms = syn_tyme.hour*10000 + syn_tyme.minute*100 + syn_tyme.second

        if filename is None:
            #filename = '%s/%s.obs.%d_%02dz.ods'%(dir,expid,nymd,nhms/10000)
            filename = '%s/%s.obs.%d.ods'%(dir,expid,nymd)

        # Stop here is no good obs available
        # ----------------------------------
        if self.nobs == 0:
            return (filename, self.nobs) # no data to work with

        # Interval for this synoptic time
        # -------------------------------
        hdt = timedelta(seconds=60*60*int(24/nsyn)/2) # usually 3/2 hours
        t1 = syn_tyme - hdt
        t2 = syn_tyme + hdt
        I = (self.tyme>=t1)&(self.tyme<t2)&(self.AOT_550>0)
        
        # Create and populated ODS object
        # -------------------------------
        lon = self.lon[I]
        nobs = len(lon)

        if nobs>0:

            ods = ODS(nobs=nobs, kx=KX, kt=KT['AOD'])

            ods.ks[:] = list(range(1,1+nobs))
            ods.lat[:] = self.lat[I]
            ods.lon[:] = self.lon[I]
            ods.qch[:] = zeros(nobs).astype('int')
            ods.qcx[:] = zeros(nobs).astype('int')
            ods.time[:] = zeros(nobs).astype('int') # fix this if needed
            ods.lev[:] = 550. * ones(nobs)

            ods.obs[:] = self.AOT_550[I]
            ods.xvec[:] = self.AOT_500[I]
            ods.xm[:] = self.AOT_675[I]
            
            if self.verb:
                print("[w] Writing file <"+filename+"> with %d observations at %dZ"%\
                   (ods.nobs,nhms/10000))

            ods.write(filename,nymd,nhms,nsyn=nsyn)

        return (filename, nobs)

#---
    def writeGridded(self, syn_tyme,
                         filename=None,dir='.',expid='aeronet',refine=8,res=None,
                         nsyn=8,doFilter=True):
       """
        Writes gridded AERONET AOD to a GFIO file.

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

       """
       from binObs_  import binobs2d, binobs3d
       from gfio     import GFIO
       
       # Interval for this synoptic time
       # -------------------------------
       if doFilter:
           hdt = timedelta(seconds=60*60*int(24/nsyn)/2) # usually 3/2 hours
           t1 = syn_tyme - hdt
           t2 = syn_tyme + hdt
           I = (self.tyme>=t1)&(self.tyme<t2)
       else:
           I = ones(self.lon.shape).astype('bool') # all that comes in
       lon = self.lon[I]
       nobs = len(lon)

       nymd = syn_tyme.year*10000 + syn_tyme.month*100  + syn_tyme.day
       nhms = syn_tyme.hour*10000 + syn_tyme.minute*100 + syn_tyme.second

       if filename is None:
           #filename = '%s/%s.sfc.%d_%02dz.nc4'%(dir,expid,nymd,nhms/10000)
           filename = '%s/%s.aod.%d.nc4'%(dir,expid,nymd)

       # Stop here is no good obs available
       # ----------------------------------
       if nobs == 0:
           return filename # no data to work with

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

       glon = linspace(-180.,180.,im,endpoint=False)
       glat = linspace(-90.,90.,jm)

       vtitle = [ 'AERONET Aerosol Optical Depth at 550nm (Level 2, interpolated)',]
       vname  = ['tau_550', ]
       vunits = [ '1',  ]
       kmvar  = [  0 ,  ]
       levs = array([550.,])
       
       title = 'Gridded AERONET Level 2 Aerosol Retrievals'
       source = 'NASA/GSFC GEOS-5 Aerosol Group (from L2 AERONET Retrievals)'
       contact = 'arlindo.dasilva@nasa.gov'

       # Create the file
       # ---------------
       if os.path.exists(filename):
           f = GFIO(filename,mode='w')
       else:
           f = GFIO()
           timinc = (24/nsyn) * 10000
           f.create(filename, vname, nymd, nhms, timinc=timinc,
                    lon=glon, lat=glat, levs=levs, levunits='nm',
                    vtitle=vtitle, vunits=vunits,kmvar=kmvar,amiss=MISSING,
                    title=title, source=source, contact=contact)

       # Grid variable and write to file
       # -------------------------------
       f.write('tau_550', nymd, nhms, 
               binobs2d(self.lon[I],self.lat[I],self.AOT_550[I],im,jm,MISSING) )

       if self.verb:
           print("[w] Wrote file "+filename+" at %02dZ"%(nhms/10000))

       return filename

#---
    def reduce(self,I):
        """
        Reduce all variables according to booleqn Index I.
        """
#        if len(I) != self.nobs:
#            raise ValueError, 'Index size <%d> must be equal to nobs <%d>'%(len(I),self.nobs)
        
        for v in self.__dict__:
            q = self.__dict__[v]
            s = shape(q)
            if len(s) == 1:
                if s[0] == self.nobs:
                    self.__dict__[v] = q[I]

        self.nobs = self.lon.size

    def bin_hourly(self,x,stn):
        """
        Plot attribute x, first binning by hourly.
        """
        I = self.Location==stn

        # create time grid
        # ----------------
        dt = timedelta(seconds=60) # one hour
        t_ = a.tyme[I]
        t0, tf = t_.min(), t_.max()
        t0 = datetime(t0.year,t0.month,t0.day,t0.hour)
        tf = datetime(tf.year,tf.month,tf.day,tf.hour)
        n = int(0.5+(tf-t0).total_seconds()/60.)
        t = array([t0 + i*dt for i in range(n)])

        print(t0, tf)
        return t

        
#---
def _convert2Float(s,missing=MISSING):
   if s in ('N/A', 'NA','','NaN',None):
      result = missing
   else:
      result = float(s or missing)
      
   return result

def _updAOT(tau1,tau2):
    """
    Replace tau1 with tau2 where tau2 is defined.
    """
    I = (tau2>0)
    tau = tau1[:]
    tau[I] = tau2[I]
    return tau

#---
def aodInterpAngs(lambda_,tau1,tau2,lambda1,lambda2):
    """                                                                                               
       Angstrom-interpolated AOD.
    """
    I = (tau1>0) & (tau2>0)
    angstrom = -log(tau1[I]/tau2[I])/log(lambda1/lambda2)
    tau = MISSING * ones(len(tau1))
    tau[I] = tau2[I] * (lambda2/lambda_)**angstrom
    return tau


#.......................................................................................
def retrieve( filename='aeronet.csv',
              request='year=2016&month=9&day=1&year2=2016&month2=9&day2=30&AOD15=1&AVG=10',
              verbose=True):
    """
    Use wget to retrieve aeronet data using the Version 3 webservice.'
    """

    cmd = 'wget --no-check-certificate -q -O %s' \
        + ' "https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3?%s&if_no_html=1"' \
        %(filename,request)

    if verbose:
        print(cmd)
        
    if os.system(cmd):
        raise ValueError("Cannot retrieve request <%cmd>")

        
#.......................................................................................
def granules( tyme,
              bracket=None,
              Version = '20',
              RootDir='/nobackup/AERONET/Level2',
              template='Y$year/M$month/aeronet_v$ver.$nymd.txt' ):
    """
    Given a date in *tyme* get files corresponding to bracketing days. On inout

     brackt = None, 'lfeft', 'right', or 'both'

    """

    oneday = timedelta(seconds=24*60*60)
    t2 = datetime(tyme.year,tyme.month,tyme.day)
    t1 = t2 - oneday
    t3 = t2 + oneday

    if bracket is None:
        Times = (t2,)
    elif bracket == 'left':
        Times = (t1,t2)
    elif bracket == 'right':
        Times = (t2,t3)
    else:
        Times = (t1,t2,t3)
    
    Files = []
    
    for t in Times:

        year = str(t.year)
        month = '%02d'%t.month
        day = '%02d'%t.day

        nymd = str(t.year*10000 + t.month*100 + t.day)

        pat = RootDir+'/'+template.replace('$year',year)\
                                  .replace('$month',month)\
                                  .replace('$ver',Version)\
                                  .replace('$nymd',nymd)

        Files += sorted(glob(pat))

    return Files
        
#..........................................................................

if __name__ == "__main__":

        a = AERONET_L2('aeronet.csv')
    
def __test__():
    
    t = datetime(2008,0o7,0o1)
    one_hour = timedelta(seconds=60*60) # synoptic dt =  3 hours

    Files = granules(t,bracket='left')

    a = AERONET_L2(Files,Verbose=True)

    for h in range(0,24,3):

        t_ = t + h * one_hour

        a.writeODS(t_)
        a.writeGridded(t_)
