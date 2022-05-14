#!/bin/env python

import os
import sys

from types     import *
from pyhdf     import SD
from netCDF4   import Dataset
from datetime  import date, datetime, timedelta
from numpy     import savez, meshgrid, array, concatenate, zeros, ones, \
                      linspace, sqrt, load, shape, random, interp
from glob      import glob

from .npz       import NPZ

SDS_ =  ( "latitude","longitude","aot1_Final","aot2_Final")

META_ = ('scan_line_time','anchor_relative_azimuth','anchor_sensor_zenith',
         'anchor_solar_azimuth','anchor_solar_zenith',
         'ch1_reflectance','ch2_reflectance', 'ch3a_reflectance')

ALIAS = dict (
                 latitude  = 'lat',
                longitude  = 'lon',
                aot1_Final = 'tau_630',
                aot2_Final = 'tau_860',
            scan_line_time = 'scanTime',
   anchor_relative_azimuth = 'SensorAzimuth',
      anchor_sensor_zenith = 'SensorZenith',
      anchor_solar_azimuth = 'SolarAzimuth',
       anchor_solar_zenith = 'SolarZenith',
           ch1_reflectance = 'ref_630',
           ch2_reflectance = 'ref_860',
          ch3a_reflectance = 'ref_1600',
    packed_pixel_meta_data = 'qc')           

MISSING = 1.E20

KX = 324 
KT = dict ( AOD = 45, )

# ..................................................................................

class AVHRR_L2B(object):

  """
      Read Level-2 AVHRR Aeorosol Optical Thickness.
  """


  def __init__ (self,Path=None,syn_time=None,nsyn=8,Verb=False,doMeta=True,
                     only_good=True,SDS=SDS_,META=META_):
     """
       Create an AVHRR object defining the attributes corresponding
       to the SDS of interest.
     """
     self.SDS  = SDS   # AOD retrievals
     self.META = META  # scan time, angles, reflectances
     self.verb = Verb
     self.only_good = only_good
     self.doMeta = doMeta
     
     if self.doMeta:
       self.Names = self.SDS + self.META + ('tyme',)
     else:
       self.Names = self.SDS
     
     # Create empty lists for SDS to be read from file
     # -----------------------------------------------
     for sds in self.Names:
         name = os.path.basename(sds) # strip the groups
         assert(name != '')
         self.__dict__[sds] = []

     # Make sure Path is a list
     # ------------------------   
     if Path is None:
            self.nobs = 0
            print("WARNING: Empty AVHRR object created")
            return         
     if type(Path) in (ListType,TupleType):
         pass # good to go
     elif type(Path) is StringType:
         if Path[-4:] == '.npz':
           self._initNPZ(Path,Verb) # reduced NPZ files
           return
         else:
           Path = [Path, ]
     else:
         raise ValueError('invalid type for Path')

     # Read each orbit, appending them to the list
     # -------------------------------------------
     self._readList(Path)

     # Make each attribute a single numpy array
     # ----------------------------------------
     for name in self.Names:
            try:
                self.__dict__[name] = concatenate(self.__dict__[name])
            except:
                print(self.__dict__[name])
                raise IndexError("Failed concatenating "+name)

     # Make aliases for compatibility with older code
     # ----------------------------------------------
     Alias = list(ALIAS.keys())
     for sds in self.Names:
         if sds in Alias:
             self.__dict__[ALIAS[sds]] = self.__dict__[sds]

     self.nobs = len(self.longitude)

  def _initNPZ(self,npzFiles,Verbose=False):
        """
        Load from NPZ file.
        """
        if type(npzFiles) == StringType:
            npzFiles = sorted(glob(npzFiles))

        if Verbose:
            print('[] Loading ', npzFiles[0]) 

        # Single file
        # -----------
        if len(npzFiles) == 1:
            f = load(npzFiles[0])
            for v in list(f.keys()):
                self.__dict__[v] = f[v]

        # Multiple files
        # --------------
        else:
            
            # For many files, process first file in list
            # ------------------------------------------
            f = load(npzFiles[0])
            V = dict()
            for v in list(f.keys()):
                if len(shape(f[v])) == 0: 
                    V[v] = [[f[v],],]
                else:
                    V[v] = [f[v],]

            # Append the other files
            # ----------------------
            for npzFile in npzFiles[1:]:
                if Verbose:
                    print('[] Loading ', npzFile)
                f = load(npzFile)
                for v in V:
                    if len(shape(f[v])) == 0: 
                        V[v].append([f[v],])
                    else:
                        V[v].append(f[v])
                f.close()

            # Concatenate
            # -----------
            for v in V:
                self.__dict__[v] = concatenate(V[v])

        # Undo erroneous scaling of metadata
        # ----------------------------------
        if self.anchor_relative_azimuth.max() > 180.:

          self.anchor_relative_azimuth = _unscale(self.anchor_relative_azimuth,0.,180.,-127,127)
          self.anchor_sensor_zenith    = _unscale(self.anchor_sensor_zenith,   0.,180.,-127,127)
          self.anchor_solar_azimuth    = _unscale(self.anchor_solar_azimuth,   0.,360.,-127,127)
          self.anchor_solar_zenith     = _unscale(self.anchor_solar_zenith,    0.,180.,-127,127)

          self.ch1_reflectance = _unscale(self.ch1_reflectance, -2.,120.,-32767,32767)
          self.ch2_reflectance = _unscale(self.ch2_reflectance, -2.,120.,-32767,32767)

        # Make aliases for compatibility with older code
        # ----------------------------------------------
        Alias = list(ALIAS.keys())
        for sds in self.Names:
          if sds in Alias:
             self.__dict__[ALIAS[sds]] = self.__dict__[sds]

        self.nobs = len(self.longitude)


  def balance(self,N):
    """
    Return indices of observations so that each species does not have more than
    N observations. This is meant to be performed with a reduced dataset.
    """
    I = zeros(self.lon.shape).astype(bool)
    random.seed(32768) # so that we get the same permutation
    for f in (self.fdu,self.fss,self.fcc,self.fsu):

      J = f>0.5                      # all obs for which species dominate
      n = len(self.lon[J])              # no. obs for this species
      P = random.permutation(n)      # randomize obs for this species
      m = min(n,N)                   # keep this many

      K = I[J]
      K[P[0:m]] = True
      I[J] = K

    return I

  def getCoxMunk(self,filename='/nobackup/NNR/Misc/avhrr.coxmunk_lut.npz',channel=630):
    """
    Returns ocean albedo as a function of wid speed.
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
    albedo = zeros(len(w10m))
    albedo[:] = interp(w10m,lut.speed,lut.albedo[:,j])

    self.albedo = albedo

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
  def _readList(self,List):
    """
      Recursively, look for files in list; list items can
      be files or directories.
      """
    for item in List:
      if   os.path.isdir(item):    self._readDir(item)
      elif os.path.isfile(item):   self._readGranule(item)
      else:
            print("%s is not a valid file or directory, ignoring it"%item)

#-------

  def _readGranule(self,filename):
      """ Read in a single granule """

      if self.verb:
        print("[] working on <%s>"%filename)

      if self.doMeta:
        Names = self.SDS + self.META
      else:
        Names = self.SDS

      g = dict()

      # Read in retrievals
      # ------------------
      f = SD.SD(filename)
      for name in self.SDS:
        
         v = f.select(name)
         data = v.get()
         Attr = v.attributes()

         rmiss = Attr['RANGE_MISSING']
         
         # Scale data if necessary
         # (see http://cimss.ssec.wisc.edu/patmosx/doc/v4_level3_format.html)
         # ------------------------------------------------------------------
         scaled = Attr['SCALED']
         if scaled:
             smiss = Attr['RANGE_MISSING']
             smin, smax = Attr['SCALED_MIN'], Attr['SCALED_MAX']
             rmin, rmax = Attr['RANGE_MIN'], Attr['RANGE_MAX']
             # print "---> scaling ", name, scaled, smin, smax, rmin, rmax
             I = (data==smiss)
             J = (data!=smiss)
             data_ = (data-smin)/float(smax-smin)
             if scaled == 1:
               data = rmin + (rmax-rmin) * data_
             elif scaled==2:
               data = 10.**(rmin + (rmax-rmin) * data_)
             elif scaled==3:
               data = rmin + (rmax-rmin) * (data_**2)
             else:
               raise ValueError("Unknow scaling algorith,")

             data[I] = rmiss # restore missing value
           
         data[data==rmiss] = MISSING # single missing value for everything
             
         g[name] = data

      # Metadata: Read in angles, time, etc.
      # ------------------------------------
      if self.doMeta:
          metafile = filename[:].replace('SFinal002_patmosx_','patmosx_sw_').replace('.hdf','.nc')
          self._readMetadata(metafile,g)
          
      # Expand lon, lat
      # ---------------
      Lon, Lat = meshgrid(g['longitude'],g['latitude'])
      g['longitude'], g['latitude'] = Lon.ravel(), Lat.ravel()

      # Quality control
      # ---------------
      if self.only_good:
        Igood = (g['aot1_Final']>-0.01) &\
                (abs(g['longitude'])<=180.) &\
                (abs(g['latitude'])<=90.)
      
      # append to object for later concatenation
      # ----------------------------------------
      for name in Names:

         data = g[name].ravel()
         
         if self.only_good:
             data = data[Igood]

         # Store data for this granule
         # ---------------------------
         self.__dict__[name].append(data)

      # Build python time for this daily granule
      # ----------------------------------------
      year = int(f.attributes()['END_YEAR'])
      doy  = int(f.attributes()['END_DAY'])
      if self.doMeta:
          scanTime = self.__dict__['scan_line_time'][-1]
          t0 = datetime(year,1,1)+(doy-1)*timedelta(seconds=24*60*60)
          tyme = array([t0 + timedelta(seconds=int(utc*60*60)) for utc in scanTime])
          self.tyme += [tyme,]
          
      # All done with file
      # ------------------
      #f.close()

#---
  def _readMetadata(self,metafile,g):
      """
      Read metafile with scan line time and angles.
      """
      nc = Dataset(metafile)
      for name in self.META:
        
         v = nc.variables[name]
         data = v[:,:]

         rmiss = v.getncattr('range_missing')
         
         # Scale data if necessary 
         # (see http://cimss.ssec.wisc.edu/patmosx/doc/v4_level3_format.html)
         # ------------------------------------------------------------------
         scaled = v.getncattr('scaled')
         scaled = 0                     # These files are not really scaled!  
         if scaled:
             smiss = v.getncattr('range_missing')
             smin, smax = v.getncattr('scaled_min'), v.getncattr('scaled_max')
             rmin, rmax = v.getncattr('range_min'), v.getncattr('range_max')
             I = (data==smiss)
             J = (data!=smiss)
             data_ = (data-smin)/float(smax-smin)
             if scaled == 1:
               data = rmin + (rmax-rmin) * data_
             elif scaled==2:
               data = 10.**(rmin + (rmax-rmin) * data_)
             elif scaled==3:
               data = rmin + (rmax-rmin) * (data_**2)
             else:
               raise ValueError("Unknow scaling algorith,")

             data[I] = rmiss # restore missing value
           
             data[data==rmiss] = MISSING # single missing value for everything
             
         g[name] = data

      # All done with file
      # ------------------
      nc.close()

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
  def writeODS(self, syn_tyme, filename=None, dir='.', expid='avhrr', nsyn=8, doNNR=False):
        """
        Writes the un-gridded AVHRR object to an ODS file.
        """
        from pyods     import ODS
        
        # Stop here is no good obs available
        # ----------------------------------
        if self.nobs == 0:
            return # no data to work with

        nymd = syn_tyme.year*10000 + syn_tyme.month*100  + syn_tyme.day
        nhms = syn_tyme.hour*10000 + syn_tyme.minute*100 + syn_tyme.second

        if filename is None:
            #filename = '%s/%s.obs.%d_%02dz.ods'%(dir,expid,nymd,nhms/10000)
            filename = '%s/%s.obs.%d.ods'%(dir,expid,nymd)

        # Interval for this synoptic time
        # -------------------------------
        dt = timedelta(seconds=60*60*int(24/nsyn)/2) # usually 3 hours
        t1 = syn_tyme - dt
        t2 = syn_tyme + dt
        I = (self.tyme>=t1)&(self.tyme<t2)
        
        # Create and populated ODS object
        # -------------------------------
        lon = self.lon[I]
        nobs = len(lon)

        ods = ODS(nobs=nobs, kx=KX, kt=KT['AOD'])

        ods.ks[:] = list(range(1,1+nobs))
        ods.lat[:] = self.lat[I]
        ods.lon[:] = self.lon[I]
        ods.qch[:] = zeros(nobs).astype('int')
        ods.qcx[:] = zeros(nobs).astype('int')
        ods.time[:] = zeros(nobs).astype('int') # fix this if needed
        if doNNR:
            ods.lev[:] = 550. * ones(nobs)
            ods.obs[:] = self.tau_550[I]
            ods.xvec[:] = self.ref_860[I]
            ods.xm[:] = self.ref_630[I]
        else:
            ods.lev[:] = 630. * ones(nobs)
            ods.obs[:] = self.tau_630[I]
            ods.xvec[:] = self.tau_860[I]
            ods.xm[:] = self.ref_630[I]
            
        if self.verb:
            print("[w] Writing file <"+filename+"> with %d observations at %dZ"%\
                   (ods.nobs,nhms/10000))

        ods.write(filename,nymd,nhms,nsyn=nsyn)

#---
  def writeGridded(self, syn_tyme,
                         filename=None,dir='.',expid='avhrr',refine=8,res=None,
                         nsyn=8,doNNR=False,doFilter=True):
       """
        Writes gridded AVHRR AOD to a GFIO file.

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
           dt = timedelta(seconds=60*60*int(24/nsyn)/2) # usually 3/2 hours
           t1 = syn_tyme - dt
           t2 = syn_tyme + dt
           I = (self.tyme>=t1)&(self.tyme<t2)
       else:
           I = ones(self.lon.shape).astype('bool') # all that comes in
       lon = self.lon[I]
       nobs = len(lon)

       # Stop here is no good obs available
       # ----------------------------------
       if nobs == 0:
           return # no data to work with

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

       nymd = syn_tyme.year*10000 + syn_tyme.month*100  + syn_tyme.day
       nhms = syn_tyme.hour*10000 + syn_tyme.minute*100 + syn_tyme.second

       vtitle = [ 'AVHRR Aerosol Optical Depth at 630nm (NOAA CDR)',]
       vname  = ['tau_630', ]
       vunits = [ '1',  ]
       kmvar  = [  0 ,  ]
       levs = array([630.,])

       if doNNR:
         vtitle += [ 'AVHRR Aerosol Optical Depth at 550nm (NASA NNR)',]
         vname  += ['tau_550', ]
         vunits += [ '1',  ]
         kmvar  += [  0 ,  ]
         levs = array([550.,])
       
       title = 'Gridded AVHRR Aerosol Retrievals'
       source = 'NASA/GSFC GEOS-5 Aerosol Group'
       contact = 'arlindo.dasilva@nasa.gov'

       if filename is None:
           #filename = '%s/%s.sfc.%d_%02dz.nc4'%(dir,expid,nymd,nhms/10000)
           filename = '%s/%s.aod.%d.nc4'%(dir,expid,nymd)

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
       f.write('tau_630', nymd, nhms, 
               binobs2d(self.lon[I],self.lat[I],self.tau_630[I],im,jm,MISSING) )
           
       if doNNR:
         f.write('tau_550', nymd, nhms, 
                 binobs2d(self.lon[I],self.lat[I],self.tau_550[I],im,jm,MISSING) )

       if self.verb:
           print("[w] Wrote file "+filename)

#---
  def reduce(self,I):
    """
    Reduce observations according to index I. 
    """
    Nicknames = list(ALIAS.values())
    for name in self.__dict__:
      if name in Nicknames:
        continue # alias do not get reduced
      q = self.__dict__[name]
      if type(q) is type(self.lon):
        if len(q) == self.nobs:
          # print "{} Reducing "+name
          self.__dict__[name] = q[I]

    Alias = list(ALIAS.keys())
    for sds in self.Names:
      if sds in Alias:
        self.__dict__[ALIAS[sds]] = self.__dict__[sds] # redefine aliases

    self.nobs = len(self.lon)

#---
  def speciate(self,aer_x,FineMode=False):
    """
    Use GAAS to derive fractional composition.
    """
    from gfio import GFIOHandle

    self.sampleFile(aer_x,onlyVars=('TOTEXTTAU',
                                    'DUEXTTAU',
                                    'SSEXTTAU',
                                    'BCEXTTAU',
                                    'OCEXTTAU',
                                    'SUEXTTAU',
                                    ))
    s = self.sample
    I = (s.TOTEXTTAU<=0)
    s.TOTEXTTAU[I] = 1.E30
    self.fdu  = s.DUEXTTAU / s.TOTEXTTAU
    self.fss  = s.SSEXTTAU / s.TOTEXTTAU
    self.fbc  = s.BCEXTTAU / s.TOTEXTTAU
    self.foc  = s.OCEXTTAU / s.TOTEXTTAU
    self.fcc  = self.fbc + self.foc
    self.fsu  = s.SUEXTTAU / s.TOTEXTTAU

    if FineMode:
      TOTEXTTAU = s.TOTEXTTAU[:]
      self.sampleFile(aer_x,onlyVars=('DUEXTTFM','SSEXTTFM'))
      self.fduf = s.DUEXTTFM / TOTEXTTAU
      self.fssf = s.SSEXTTFM / TOTEXTTAU

    del self.sample

#---
  def sampleG5(self,gas_x=None,avk_x=None,int_x=None,slv_x=None,ext_Nc=None):
    """
    Sample key parameters from GAAS files.
    """
    from gfio import GFIOHandle

    if gas_x is not None:
      self.sampleFile(gas_x,onlyVars=('AODANA',))
      self.tau_550 = self.sample.AODANA[:]

    if avk_x is not None:
      tyme = self.tyme[:]
      self.tyme = getSyn(tyme)
      self.sampleFile(avk_x,onlyVars=('AOD',))
      self.avk = self.sample.AOD[:]
      self.tyme[:] = tyme[:]

    if int_x is not None:
      try:
        self.sampleFile(int_x,onlyVars=('TQV',)) # As in file spec
        self.tpw = self.sample.TQV[:]
      except:
        self.sampleFile(int_x,onlyVars=('TPW',)) # Larry's name
        self.tpw = self.sample.TPW[:]

    if slv_x is not None:
      self.sampleFile(slv_x,onlyVars=('U10M','V10M'))
      self.wind = sqrt(self.sample.U10M[:]**2 + self.sample.V10M[:]**2)

    if ext_Nc is not None:
      self.sampleFile(ext_Nc,onlyVars=('taod',))
      self.tau_660 = self.sample.taod[:,5] # 660


    del self.sample

#---
  def sampleFile(self, inFile, npzFile=None, onlyVars=None, Verbose=False):
        """
        Interpolates all variables of inFile and optionally
        save them to file *npzFile*
        """
        from gfio import GFIO, GFIOctl, GFIOHandle

        # Instantiate grads and open file
        # -------------------------------
        name, ext = os.path.splitext(inFile)
        if ext in ( '.nc4', '.nc', '.hdf'):
          fh = GFIO(inFile)     # open single file
          if fh.lm == 1:
            timeInterp = False    # no time interpolation in this case
          else:
            raise ValueError("cannot handle files with more tha 1 time, use ctl instead")
        else:
          fh = GFIOctl(inFile)  # open timeseries
          timeInterp = True     # perform time interpolation

        self.sample = GFIOHandle(inFile)
        if onlyVars is None:
            onlyVars = fh.vname

        nt = self.lon.shape
        tymes = self.tyme
        lons = self.lon
        lats = self.lat

        # Loop over variables on file
        # ---------------------------
        for v in onlyVars:
            if Verbose:
                print("<> Sampling ", v)
            if timeInterp:
              var = fh.sample(v,lons,lats,tymes,Verbose=Verbose)
            else:
              var = fh.interp(v,lons,lats)
            if len(var.shape) == 1:
                self.sample.__dict__[v] = var
            elif len(var.shape) == 2:
                var = var.T # shape should be (nobs,nz)
                self.sample.__dict__[v] = var
            else:
                raise IndexError('variable <%s> has rank = %d'%(v,len(var.shape)))

        if npzFile is not None:
            savez(npzFile,**self.sample.__dict__)            

       
  def sampleLoadz(self,npzFile):
        """
        Loads sample from npz file.
        """
        from grads.gahandle import GaHandle
        self.sample = GaHandle(npzFile)
        npz = load(npzFile)
        for v in list(npz.keys()):
            self.sample.__dict__[v] = npz[v]
                
#.......................................................................................
def granules(sat, orb, tyme,
              bracket=True,
              RootDir='/Users/adasilva/data/AVHRR/Level2B',
              template='$year/$doy/SFinal002_patmosx_$sat_$orb_$year_$doy_v05r02.hdf'):
    """
    Given a date in *tyme* get files corresponding to bracketing days.
    """

    oneday = timedelta(seconds=24*60*60)
    t2 = datetime(tyme.year,tyme.month,tyme.day)
    t1 = t2 - oneday
    t3 = t2 + oneday

    if bracket:
        Times = (t1,t2,t3)
    else:
        Times = (t2,)
    
    Files = []
    
    for t in Times:

        dt = (t - datetime(t.year,1,1))
        doy = '%03d'%int(1 + dt.total_seconds() / oneday.total_seconds())
        year = str(t.year)
        pat = RootDir+'/'+template.replace('$year',year)\
                                  .replace('$doy',doy)\
                                  .replace('$sat',sat)\
                                  .replace('$orb',orb)

        Files += sorted(glob(pat))

    return Files
        
#----------

def _count_des():

  RootDir = '/nobackup/AVHRR/Level2/PATMOSX'

  Files = sorted(glob(RootDir+'/????/???/*_des_*.hdf'))

  f = open('des_inventory.txt','w')
  for fname in Files:
    a = AVHRR_L2B(fname,doMeta=False,Verb=False)
    if a.nobs>0:
       tokens = os.path.basename(fname).split("_")
       sat, orb, year, doy = tokens[2:6]
       line = "%s %s %s %s ... %5d AOT retrievals"%(orb, sat, year, doy,a.nobs)
       print(line)
       f.write(line+'\n')

  f.close()

def getSyn(tyme, nsyn=8):
  """
  Return synoptic time. Can be optimized.
  """
  dtsyn = 24/nsyn
  dt_secs = dtsyn * 60. * 60.
  oneday = timedelta(seconds=24*60*60)
  syn_time = []
  for t in tyme:
      sod = t.hour*3600 + t.minute * 60 + t.second
      isyn = int(0.5+sod/dt_secs)
      if isyn<nsyn:
         t_ = datetime(t.year,t.month,t.day,dtsyn*isyn)
      elif isyn==nsyn:
         t_ = datetime(t.year,t.month,t.day) + oneday
      else:
         raise ValueError('invalid isyn')
      syn_time += [t_,]

  return array(syn_time)

def _unscale(x,rmin,rmax,smin,smax):
  """
  Undo linear scaling.
  """
  r = (smax-smin)/(rmax-rmin)
  x_ = smin + r * (x-rmin)
  return x_

#---
if __name__ == "__main__":

    a = AVHRR_L2B('/nobackup/AVHRR/Level2/NPZ/2008/*.npz',Verb=True)

def xxx():

    # _count_des()

    RootDir = '/nobackup/AVHRR/Level2/PATMOSX'
    MAerDir = '/nobackup/MERRAero'
    MDir = '/nobackup/MERRA'

    gas_x = MAerDir + '/inst2d_gaas_x.ddf'
    aer_x = MAerDir + '/tavg2d_aer_x.ddf'
    avk_x = MAerDir + '/inst2d_avk_x.ddf'
    ext_Nc = MAerDir + '/ext_Nc.ddf'
    int_x = MDir + '/int_Nx'
    slv_x = MDir + '/slv_Nx'

    #tyme = datetime(2008,6,9)
    tyme = datetime(1981,11,9)
    Files = granules('n??','asc',tyme,RootDir=RootDir,bracket=False)

    print("Files: ", Files)

    #tyme = datetime(2008,6,9)
    #Files += granules('n??','asc',tyme,RootDir=RootDir,bracket=True)
    
    a = AVHRR_L2B(Files,Verb=True)
    #a.speciate(aer_x)
    #a.sampleG5(gas_x,avk_x,int_x)
    #a.sampleG5(slv_x=slv_x)

def later():

    for syn_hour in range(0,24,3):
        syn_tyme = tyme + timedelta(seconds=syn_hour*60*60)
        #a.writeODS(syn_tyme)
        a.writeGridded(syn_tyme)


    a.sampleFile(gas_x)
    #b = AVHRR_L2B(Files,Verb=True,doMeta=False)
    for syn_hour in range(0,24,3):
        syn_tyme = tyme + timedelta(seconds=syn_hour*60*60)
        #a.writeODS(syn_tyme)
        a.writeGridded(syn_tyme)
