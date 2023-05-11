"""
Reads MODIS Level 2 MOD06/MYD06 Cloud Product granules and returns
a single object with the relevant data. Can read a single granule,
a directory of granules (such as for a whole day of data), or a list
of either or both of the above. Directory reads proceed recursively
to subdirectories.

This software is hereby placed in the public domain.
Arlindo.daSilva@nasa.gov

08/??/16 PMN General improvements
09/26/16 PMN 'f' grid capability
09/27/16 PMN Collection 6 version, also uses masked array intermediates
10/24/16 PMN timinc for GFIO, im and jm options to write
"""

import os
import sys
import numpy as np
import numpy.ma as ma

from types     import *
from datetime  import datetime, timedelta
from glob      import glob

from pyhdf.SD import SD, HDF4Error

#from bits import BITS

#---  

# Time Reference used by MODIS granules
DATE_START = datetime(1993,1,1,0,0,0)

# Scientific Data Sets 
# COP = Cloud Optical Properties
SDS = dict (
      META = ('Longitude', 'Latitude', 'Scan_Start_Time',
#             'Sensor_Zenith', 'Sensor_Azimuth',
#             'Solar_Zenith',  'Solar_Azimuth',
            ),
      DATA = ( 'Cloud_Top_Pressure', 
               'Cloud_Top_Temperature',
#               'Brightness_Temperature',
#               'Cloud_Fraction',
##              'Cloud_Phase_Infrared',	# byte
##              'Cloud_Height_Method',	# byte
##              'Cloud_Mask_5km', 	# byte	# 3D, 3rd dimension is #bytes in mask
##              'Surface_Type', 		# not in c006, in Cloud_Mask???
            ),
      COP = ( 
#              'Cloud_Effective_Radius',
#              'Cloud_Optical_Thickness',
#              'Cloud_Water_Path',
##             'Cloud_Phase_Optical_Properties',		# byte
##             'Cloud_Multi_Layer_Flag',			# byte
##             'Cloud_Effective_Radius_Uncertainty',
##             'Cloud_Mask_1km',				# byte  # 3D, 3rd dimension is #bytes in mask
##             'Quality_Assurance_1km',			# byte  # 3D, 3rd dimension is #bytes in flag
            ),
      )

ALIAS = dict (              Longitude = 'lon',
                             Latitude = 'lat',
                        Sensor_Zenith = 'SensorZenith',
                       Sensor_Azimuth = 'SensorAzimuth',
                         Solar_Zenith = 'SolarZenith',
                        Solar_Azimuth = 'SolarAzimuth',
                     Scattering_Angle = 'ScatteringAngle',
                   Cloud_Top_Pressure = 'CTP',  
                Cloud_Top_Temperature = 'CTT',
               Brightness_Temperature = 'BT',
                       Cloud_Fraction = 'F',
               Cloud_Effective_Radius = 'RE',
              Cloud_Optical_Thickness = 'TAU',
                     Cloud_Water_Path = 'CWP',
       Cloud_Phase_Optical_Properties = 'PHASE',
               Cloud_Multi_Layer_Flag = 'MLAYER',
   Cloud_Effective_Radius_Uncertainty = 'eRE',
             )

BAD, MARGINAL, GOOD, BEST = ( 0, 1, 2, 3 )  # QA marks

MISSING = -999.999

#...........................................................................
class MxD06Handle(object):
    """
    Generic container for MxD06.
    """
    def __init__(self,name):
        self.name = name
        self.nobs = 0
        
class MxD06Error(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
  
class MxD06_L2(object):
    """
    This class implements the MODIS Level 2 CLOUD product, usually
    referred to as MOD06 (TERRA satellite) and MYD06 (AQUA satellite).
    """

    def __init__ (self,Path,decimate='thinning',Verb=0,only_good=True):
       """
       Reads all Level 2 MOD06/MYD06 granules present on a given *Path*
       and returns a single object with all data concatenated for 5 km
       and 1 km. On input, 

       Required parameters:
         Path -- can be a single file, a single directory, of a list
                 of files and directories. Directories are transversed
                 recursively. If a non MOD06/MYD06 Level 2 file is
                 encountered, it is simply ignored.

       Optional parameters:
         decimate  --- decimate 1 km data, 'binning' or 'thinning' it to 5 km
         only_good --- keep only *good* observations
         Verb      --- Verbose level:
                       0 - really quiet (default)
                       1 - basic commentary, warns if invalid file is found
                       2 - detailed commentary
       """

       self.verb = Verb
       self.nobs = 0
       self.sat  = None # Satellite name
       self.coll = None # collection, e.g., 005
       self.SDS = SDS['META'] + SDS['DATA'] + SDS['COP']
       self.decimate = decimate
       self.iFilter = None  # used to filter I/O
       self.missing = MISSING

       # Resolution dependent containers
       # -------------------------------
       self.COP = MxD06Handle('1km')

       # Create an empty list for each SDS
       # ---------------------------------
       for name in self.SDS:
           self._setName(name,[])

       # Read each granule, appending to the SDS lists
       # ---------------------------------------------
       if type(Path) is ListType:
           if len(Path) == 0:
               print("WARNING: Empty MxD06_L2 object created")
               return
       else:
           Path = [Path, ]

       # read all the granules,
       #   (doing decimation in the process, if requested)
       self._readList(Path)

       # Concatenate each SDS list along swath to a single numpy array
       # -------------------------------------------------------------
       for sds in self.SDS:
           try:
               self._catName(sds)
               if self.verb>1:
                   print("() Concatenated <%s>"%sds)
           except:
               print("WARNING: Failed concatenating "+sds)

       # If decimating to 5km, provide access in 5km section
       # ---------------------------------------------------
       if decimate is not None:
             for name in SDS['COP']:
                   self.__dict__[name] = self.COP.__dict__[name] 
                   if self.verb>1:
                        print("/\ Moved <%s> up"%name)
       
       # Make aliases for compatibility with older code 
       # ----------------------------------------------
       Alias = list(ALIAS.keys())
       for sds in self.SDS:
           if sds in Alias:
               if self.verb>1:
                        print(">< Aliasing %s --> %s"%(sds,ALIAS[sds]))
               self._aliasName(sds)

       # Create corresponding python time
       # --------------------------------
       im, jm = self.Scan_Start_Time.shape
       T = self.Scan_Start_Time[:,0]
       Time = np.array([DATE_START+timedelta(seconds=s) for s in T])
       self.Time = Time.repeat(jm).reshape((im,jm)) # just like Scan_Start_Time

#---
    def _setName(self,name,value):
        """Assign the value to SDS 'name'."""

        if name in SDS['COP']:
              self.COP.__dict__[name] = value
        else:
              self.__dict__[name] = value

#---
    def _appendName(self,name,value):
        """Append the value to SDS 'name'."""

        if name in SDS['COP']:
              self.COP.__dict__[name].append(value)
        else:
              self.__dict__[name].append(value)

#---
    def _catName(self,name):
        """Concatenate the granules in 'name' along swath to a single numpy array."""

        # reference the SDS 'name' as V
        # V is a list of the SDS for many granules
        if name in SDS['COP']:
              b = self.COP
        else:
              b = self
        V = b.__dict__[name]
 
# needs fixing for some multi-byte flags, where multi-byte dimension comes last
# would be better to detect along swath dimension from dimension attributes

        # find the *along* swath axis of the SDS
        rank = V[0].ndim
        if rank==2:
              # (along,across)
              along = 0
        elif rank==3:
              # ("channel",along,across)
              along = 1
        else:
              raise MxD06Error('Invalid SDS rank=%d for <%s> '%(rank,name))

        # concatenate along swath to a single numpy.ma array
        b.__dict__[name] = ma.concatenate(V,axis=along)

#---
    def _aliasName(self,name):
        if self.decimate is None and name in SDS['COP']:
              self.COP.__dict__[ALIAS[name]] = self.COP.__dict__[name] 
        else:
              self.__dict__[ALIAS[name]] = self.__dict__[name] 

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
                print("%s is not a valid file or directory, ignoring it" % item)

#---
    def _readDir(self,dir):
        """Recursively, look for files in directory."""
        for item in os.listdir(dir):
            path = dir + os.sep + item
            if   os.path.isdir(path):    self._readDir(path)
            elif os.path.isfile(path):   self._readGranule(path)
            else:
                print("%s is not a valid file or directory, ignoring it" % item)

#---
    def _readGranule(self,filename):
        """Reads one MOD06/MYD06 granule with Level 2 cloud data."""

        # Don't fuss if the file cannot be opened
        # ---------------------------------------
        try:
            if self.verb > 0:
                print(" [] Working on " + filename)
            hfile = SD(filename)
        except HDF4Error:
            if self.verb > 0:
                print("- %s: not recognized as an HDF file" % filename)
            return 

        # default 5km sampling if require decimation
        # ------------------------------------------
        if self.decimate is not None:
          sampFound = False
          lSamp = [3, 2028, 5]
          xSamp = [3, 1348, 5]

        # Read selected variables
        # -----------------------
        for sds in self.SDS:
            if self.verb > 1: print(' <> Doing %s' % sds)

            # read and scale the current SDS
            # (uses masked array as best way to preserve data type)
            v = hfile.select(sds).get()
            a = hfile.select(sds).attributes()
            v = ma.masked_equal(v,a['_FillValue'])
            if a['scale_factor'] != 1 or a['add_offset'] != 0:
                if self.verb > 1: print(' <> Scaling %s' % sds)
                v = a['scale_factor'] * (v - a['add_offset'])
                # PS: I know that looks strange.
                #   The add_offset is what was used to get the HDF stored values within range.
                #   So in backing out the physical value, it must actually be subtracted!

            # append the current SDS to its list of granules
            self._appendName(sds,v) 

            # preparation for decimation
            if self.decimate is not None:

              # recorded SDS sampling
              lSamp_ = a['Cell_Along_Swath_Sampling']
              xSamp_ = a['Cell_Across_Swath_Sampling']

              # try to detect a common 5-km sampling
              if sds in SDS['META'] or sds in SDS['DATA']:
                if lSamp_[2] == 5 and xSamp_[2] == 5:
                  if not sampFound:
                    # overwrite default
                    lSamp = lSamp_
                    xSamp = xSamp_
                    sampFound = True
                  else:
                    if lSamp_ != lSamp or xSamp_ != xSamp:
                      print('  WARNING: mixed 5-km Sampling found')
                      print('     along:', lSamp, lSamp_)
                      print('    across:', xSamp, xSamp_)
                      print('    Will use first')

              # also verify that COP 1-km sampling starts from 1
              #   otherwise it will mess up decimation assumptions
              if sds in SDS['COP']:
                if lSamp_[0] != 1 or lSamp_[2] != 1:
                  raise MxD06Error('Bad 1-km  along sampling for ' + sds + lSamp_)
                if xSamp_[0] != 1 or xSamp_[2] != 1:
                  raise MxD06Error('Bad 1-km across sampling for ' + sds + xSamp_)

        # Decimate 1km data to lower 5km resolution (Needs better Q/C)
        # ------------------------------------------------------------
        # Note: [-1] below gets the current granule (that we just appended above)
        # TODO: Only 2D SDSs are currently implemented for decimation. This should
        #   be all we need, but otherwise add a 3D case with a "channel" loop.

        if self.decimate is not None:

          # report 5km sampling
          if self.verb > 1:
            print(' 5km sampling:', end=' ')
            if sampFound: print('detected:', end=' ')
            else:         print('default:', end=' ')
            print('along', lSamp, 'across', xSamp)

          # verify assumed 5-km sampling, since binning
          #   currently hardwired for -2:+2 for simplicity
          if lSamp[2] != 5 or xSamp[2] != 5:
            raise MxD06Error('Not 5-km sampling as assumed!')

          # adjust to python zero-based indexing
          lSamp[0] -= 1; lSamp[1] -= 1
          xSamp[0] -= 1; xSamp[1] -= 1

          if self.decimate == 'thinning':
            # thin 1-km data to 5-km
            for name in SDS['COP']:
              v = self.COP.__dict__[name][-1]
              if v.ndim != 2:
                raise MxD06Error('Only 2D SDSs implemented for decimation')
              self.COP.__dict__[name][-1] = v[lSamp[0]:lSamp[1]+1:5,
                                              xSamp[0]:xSamp[1]+1:5]

          elif self.decimate == 'binning':
            # average 5x5 1-km pixel blocks
            for name in SDS['COP']:
              v = self.COP.__dict__[name][-1]
              if v.ndim != 2:
                raise MxD06Error('Only 2D SDSs implemented for decimation')
              # calculate number of 5x5 boxes
              nl = (lSamp[1]-lSamp[0])//5+1
              nx = (xSamp[1]-xSamp[0])//5+1
              # assume all 5x5 boxes are complete
              v = v[lSamp[0]-2:lSamp[1]+3,xSamp[0]-2:xSamp[1]+3]
              # mean over 5x5 boxes
              self.COP.__dict__[name][-1] = \
                v.reshape(nl,5,nx,5).swapaxes(1,2).reshape(nl,nx,-1).mean(axis=-1)

        # Core Metadata
        # -------------
        cm = hfile.attributes()['CoreMetadata.0']

        # Satellite name
        # --------------
        sat = cm.split('ASSOCIATEDPLATFORMSHORTNAME')[1].split('\n')[3].split('=')[1]
        sat = sat.lstrip().replace('"','')
        if self.sat is None:
          self.sat = sat
        else:
          if sat != self.sat:
            raise MxD06Error('Mixed satellites encountered')
 
        # Collection
        # ----------
        coll = int(cm.split('COLLECTION')[1].split('VERSIONID')[1].split('\n')[2].split('=')[1])
        coll = "%03d"%coll
        if self.coll is None:
          self.coll = coll
        else:
          if coll != self.coll:
            raise MxD06Error('Mixed collections encountered')
        
#---
    def write(self, filename, syn_time, dt=timedelta(hours=6),
              dir='.', iFilter=None, refine=8, res=None, im=None, jm=None,
              Verb=1):
       """
        Writes gridded MODIS measurements to file.

         filename -- full path of filename to write to,
                       or None to use the standardized filename:
                           <dir>/<prod>_c<coll>.cld_Nx.yyyymmdd_hhmmz.nc4
                       In that case, provide the directory <dir> (default '.').
                       <prod> = [MOD06|MYD06] is obtained from the object's recorded satellite.
                       The collection <coll> is the object's recorded collection.

         syn_time -- python datetime appropriate for the MXD06_L2 object
         dt       -- a timedelta used to set GFIO timinc

         refine  -- refinement level for a base 4x5 GEOS-5 grid
                       refine=1  produces a  4     x 5      grid
                       refine=2  produces a  2     x 2.5    grid
                       refine=4  produces a  1     x 1.25   grid
                       refine=8  produces a  0.5   x 0.625  grid
                       refine=16 produces a  0.25  x 0.3125 grid
                       refine=32 produces a  0.125 x 0.125  grid <- NB!
                   NOTE: starting refine=32 the base grid is a 4x4.

        Alternatively, one can specify the grid resolution with a
        single letter:

         res     -- single letter denoting GEOS-5 resolution,
                       res='a'  produces a  4     x 5      grid
                       res='b'  produces a  2     x 2.5    grid
                       res='c'  produces a  1     x 1.25   grid
                       res='d'  produces a  0.5   x 0.625  grid
                       res='e'  produces a  0.25  x 0.3125 grid
                       res='f'  produces a  0.125 x 0.125  grid <- NB!

                   NOTE: *res*, if specified, supersedes *refine*.

         Alternatively, one can specify the grid resolution directly
         with explicit im and jm. In that case, both must be present,
         and *both res and refine are superceeded*.

         Verb -- Verbose level:
                 0 - really quiet
                 1 - basic commentary (default)

       """
       from gfio import GFIO

       if self.sat is None:
         # nothing to write
         return

       # Optional filter
       # ---------------
       self.iFilter = iFilter
       
       # Output grid resolution
       # ----------------------
       if im is None and jm is None:

         # neither im or jm provided, so use refine, 
         # possibly superceeded by res.

         if res is not None:
           if res=='a': refine = 1 
           if res=='b': refine = 2
           if res=='c': refine = 4
           if res=='d': refine = 8
           if res=='e': refine = 16
           if res=='f': refine = 32

         if refine >= 32:
           dx = 4. / refine
           dy = 4. / refine
         else:
           dx = 5. / refine
           dy = 4. / refine

         im = int(360. / dx)
         jm = int(180. / dy + 1)

       else:

         # both im and jm must be explicitly provided
         if im is None or jm is None:
           raise MxD06Error('Must provide both im and jm is either is provided!')

       # Lat lon grid
       # ------------
       glon = np.linspace(-180.,180.,im,endpoint=False)
       glat = np.linspace(-90.,90.,jm)

       # file details
       # ------------
       nch = 7
       levs = np.linspace(1,7,7) # for now   ## <<<< needs fixing and levunits below

       t = syn_time
       Y, M, D, h, m, s = (t.year,t.month,t.day,t.hour,t.minute,t.second)
       nymd = Y*10000+M*100+D
       nhms = h*10000+m*100+s

       ds = dt.seconds
       dh = ds//3600; ds -= dh*3600
       dm = ds//60;   ds -= dm*60
       timinc = dh*10000+dm*100+ds

       vtitle = [ 'Cloud_Top_Pressure', 
                  'Cloud_Top_Temperature',
                ]
#                 'Brightness_Temperature',
#                 'Cloud_Fraction',
#                 'Cloud_Effective_Radius',
#                 'Cloud_Optical_Thickness',
#                 'Cloud_Water_Path',
#                ]

       vname  = ['CTP', 'CTT',] #  'BT', 'F', 'RE',     'TAU', 'CWP'   ]
       vunits = ['Pa',  'K',  ] #  'K',  '1', 'micron', '1',   'g m-2' ]
       kmvar  = [ 0,     0,   ] #  nch,   0,   0,        0,     0      ]
#should use dictionaries

       if self.sat.lower() == 'terra':
         prod = 'MOD06'
       elif self.sat.lower() == 'aqua':
         prod = 'MYD06'
       else:
         raise MxD06Error('Unknown satellite <%s>' % self.sat)

       if self.coll is None: raise MxD06Error('Missing collection')

       title = 'Gridded MODIS Cloud Retrievals'
       source = '%s_L2 collection %s gridded at NASA/GSFC/GMAO' % (prod, self.coll)
       contact = 'peter.norris@nasa.gov'

       # filename
       # --------
       if filename is None:
           filename = os.sep.join((dir,
             '%s_c%s.cld_Nx.%08d_%04dz.nc4' % (prod,self.coll,nymd,nhms/100)))

       if Verb > 0:
           print("[w] Writing", filename)

       # Create the file
       # ---------------
       f = GFIO()
       f.create(filename, vname, nymd, nhms, timinc=timinc,
                lon=glon, lat=glat, levs=levs, levunits='nm',
                vtitle=vtitle, vunits=vunits,kmvar=kmvar,amiss=MISSING,
                title=title, source=source, contact=contact)

       # Grid variables and write to file
       # -------------------------------
       for oVar in ('CTP','CTT',): #'BT','F','RE','TAU','CWP'):
         if Verb > 1: print(' [] writing variable', oVar)
         f.write (oVar, nymd, nhms, self._binObs(self.__dict__[oVar],im,jm))

       f.close()

       # Reset filter
       # ------------
       self.iFilter = None

    def _binObs(self,q,im,jm):
      """
      Bins observations. Assumes a global GEOS-5 A-Grid: Lons in [-180,10) Lats in [-90,90].
      q can be 2D or 3D. If 3D, channel is first dimension and each channel must be valid to include in binning.
      """

      # detect 3D channel variable
      if q.ndim == 3:
        nc = q.shape[0]
      else:
        nc = None

      # remove any masked values or out-of-range values
      # a. acceptable lat, lon
      mask = self.lon.mask | self.lat.mask
      mask |= np.abs(self.lon)>180.
      mask |= np.abs(self.lat)> 90.
      # b. apply optional pass filter
      if self.iFilter is not None:
        mask |= np.logical_not(self.iFilter)
      # c. mask out missing obs
      if nc is None:
        mask |= q.mask
      else:
        # all chamnels must be valid to pass
        for n in range(nc):
          mask |= q[n].mask
      # d. mask out and make 1D
      lon = ma.array(self.lon, mask=mask, keep_mask=False).compressed()
      lat = ma.array(self.lat, mask=mask, keep_mask=False).compressed()
      if nc is None:
        obs = ma.array(q, mask=mask, keep_mask=False).compressed()
      else:
        # form a list of the masked channels, each 1D
        obs = []
        for n in range(nc):
          obs.append(ma.array(q[n], mask=mask, keep_mask=False).compressed())
        # put 1D index first, channels second
        obs = np.array(obs).T

      # bins
      lon[lon>=180.] -= 360.
      dLon = 360. / im
      dLat = 180. / (jm - 1)
      ivals = np.rint((lon + 180.) / dLon)
      jvals = np.rint((lat +  90.) / dLat)
      ivals[ivals>=im] -= im
      ivals[ivals< 0 ] += im

      # do the binnning
      cnts = np.zeros((im,jm),dtype=int)
      if nc is None:
        gObs = np.zeros((im,jm))
        for i,j,v in zip(ivals,jvals,obs):
          gObs[i,j] += v
          cnts[i,j] += 1
      else:
        gObs = np.zeros((im,jm,nc))
        for i,j,channels in zip(ivals,jvals,obs):
          gObs[i,j] += channels
          cnts[i,j] += 1

      # normalize & return
      hasData = (cnts > 0)
      if nc is None:
        gObs[hasData] /= cnts[hasData]
        gObs[~hasData] = self.missing
      else:
        gObs = np.transpose(gObs,(2,0,1)) # now (nc,im,jm)
        for n in range(nc):
          gObs[n][hasData] /= cnts[hasData]
          gObs[n][~hasData] = self.missing
        gObs = np.transpose(gObs,(1,2,0)) # now (im,jm,nc)
      return gObs

#---
    def addVar(self,ctlfile,vname,**kwds):
        """
        Interpolates variable <vname> in <ctlfile> to pixel locations
        and saves it as an attribute named <vname>.
        """
        from gfio import GFIOctl
        f = GFIOctl(ctlfile)
        self.__dict__[vname] = f.interpXY(vname,self.lon,self.lat,**kwds)
        f.close()
        
#............................................................................

def granules ( path, prod, syn_time, coll='006', nsyn=8 ):
    """
    Returns a list of MxD06 granules for a given product at given synoptic time.
    On input,

    path      ---  mounting point for the MxD06 Level 2 files
    prod      ---  either MOD06 or MYD06
    syn_time  ---  synoptic time (timedate format)

    coll      ---  collection: 005, 051, 006 (optional)
    nsyn      ---  number of synoptic times per day (optional)

    """

    # Determine synoptic time range
    # -----------------------------
    dt = timedelta(seconds = 12. * 60. * 60. / nsyn)  # half of time window
    t1, t2 = (syn_time-dt,syn_time+dt)  # time window boundaries

    # Find MODIS granules in synoptic time range
    # ------------------------------------------
    dtGranule = timedelta(minutes=5)
    t = datetime(t1.year,t1.month,t1.day,t1.hour,0,0)
    Granules = []
    while t < t2:
        if t >= t1:
            doy = t.timetuple()[7]
            basen = "%s/%s/%04d/%03d/%s_L2.A%04d%03d.%02d%02d.%s.*.hdf"\
                     %(path,prod,t.year,doy,prod,t.year,doy,t.hour,t.minute,coll)
            try:
                filen = glob(basen)[0]
                Granules += [filen,]
            except:
                pass
        t += dtGranule

    if len(Granules) == 0:
        print("WARNING: no %s collection %s granules found for time" % (prod, coll), syn_time)

    return Granules

#............................................................................

if __name__ == "__main__":

      Collection = '006'
      syn_time = datetime(2011,7,1,12,0,0)
      L2mount = os.sep.join((os.environ['NOBACKUP'],'MODIS',Collection,'Level2'))
      Files = granules(L2mount,'MOD06',syn_time,nsyn=nsyn,coll=Collection)
      m = MxD06_L2(Files,Verb=2,decimate='binning')
      m.write(None,syn_time,res='f',Verb=2)

def hold():

      m = MxD06_L2('/nobackup/MODIS/051/Level2/MOD06/2012/228',Verb=1)

      Files = sorted(glob('/nobackup/MODIS/051/Level2/MOD06/2012/228/*.hdf'))
      Files = '/nobackup/MODIS/051/Level2/MOD06/2012/228/MOD06_L2.A2012228.1200.051.2012228213410.hdf'
      m = MxD06_L2(Files,decimate='thinning',Verb=1)
      m = MxD06_L2(Files,decimate='binning',Verb=1)

      aer_Nv = '/nobackup/fp/opendap/aer_Nv.ddf'
      m.addVar(aer_Nv,'DU001',kbeg=36,kount=36)
#     m.write(None,syn_time)

    
