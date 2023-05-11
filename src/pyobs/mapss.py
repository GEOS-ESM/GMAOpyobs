"""
   Classes for reading CSV files produced by MAPSS.
"""

import os
import sys
from   types    import *
from   numpy    import loadtxt, ones, zeros,  savez, pi, log, interp, \
                       concatenate, arange, savez, shape, array
from   datetime import datetime as TIME

from .npz       import NPZ
from .mcd43gf   import MCD43GF
from .igbp      import IGBP

MISSING = -999.

BAD, MARGINAL, GOOD, BEST = list(range(4))

VARS = dict (
       META = ( 'Date', 'Time', 'Longitude', 'Latitude', 'Location' ),
       ANET = ( 'mean_AOD1640', 'mean_AOD1020', 'mean_AOD0870', 'mean_AOD0670',  
                'mean_AOD0550', 'mean_AOD0530', 'mean_AOD0500', 'mean_AOD0490', 
                'mean_AOD0440', 'mean_AOD0410', 
                'mean_AOD0380', 'mean_AOD0340' ),
   DEEP_AOD = ( 'SolarZenith', 'SolarAzimuth', 'SensorZenith',
                'SensorAzimuth', 'ScatteringAngle',
                'mean_AOD0412dpbl-l',
                'mean_AOD0470dpbl-l',
                'mean_AOD0660dpbl-l',
                'QAdpbl-l'),
  DEEP_MREF = ( 'SolarZenith', 'SolarAzimuth', 'SensorZenith',
                'SensorAzimuth', 'ScatteringAngle',
                'mean_mref0412dpbl-l',
                'mean_mref0470dpbl-l',
                'mean_mref0660dpbl-l',
                'QA-l'),
  DEEP_SREF = ('mean_surfre0412dpbl-l',
               'mean_surfre0470dpbl-l',
               'mean_surfre0660dpbl-l',
               'QA-l',),
  MISR_GEOM = ( 'irowc','icolc',
                'SolarZenith','SensorZenith1','SensorZenith2','SensorZenith3','SensorZenith4',
                'SensorZenith5','SensorZenith6','SensorZenith7','SensorZenith8','SensorZenith9',
                'RelativeAzimuth1','RelativeAzimuth2','RelativeAzimuth3','RelativeAzimuth4',
                'RelativeAzimuth5','RelativeAzimuth6','RelativeAzimuth7','RelativeAzimuth8',
                'RelativeAzimuth9',
                'ScatteringAngle1','ScatteringAngle2','ScatteringAngle3','ScatteringAngle4',
                'ScatteringAngle5','ScatteringAngle6','ScatteringAngle7','ScatteringAngle8',
                'ScatteringAngle9',
                'GlitterAngle1','GlitterAngle2','GlitterAngle3','GlitterAngle4',
                'GlitterAngle5','GlitterAngle6','GlitterAngle7','GlitterAngle8','GlitterAngle9',),
  MISR_MREF = ('irowc','icolc',
               'mean_Refl0446eq1', 'mean_Refl0446eq2', 'mean_Refl0446eq3',
               'mean_Refl0446eq4', 'mean_Refl0446eq5', 'mean_Refl0446eq6',
               'mean_Refl0446eq7', 'mean_Refl0446eq8',
               'mean_Refl0446eq9', 'mean_Refl0558eq1', 'mean_Refl0558eq2',
               'mean_Refl0558eq3', 'mean_Refl0558eq4', 'mean_Refl0558eq5',
               'mean_Refl0558eq6', 'mean_Refl0558eq7', 'mean_Refl0558eq8',
               'mean_Refl0558eq9',
               'mean_Refl0672eq1', 'mean_Refl0672eq2', 'mean_Refl0672eq3',
               'mean_Refl0672eq4', 'mean_Refl0672eq5', 'mean_Refl0672eq6',
               'mean_Refl0672eq7', 'mean_Refl0672eq8', 'mean_Refl0672eq9',
               'mean_Refl0866eq1', 'mean_Refl0866eq2', 'mean_Refl0866eq3',
               'mean_Refl0866eq4', 'mean_Refl0866eq5', 'mean_Refl0866eq6',
               'mean_Refl0866eq7', 'mean_Refl0866eq8', 'mean_Refl0866eq9',),

   MISR_AOD = ('irowc','icolc',
               'mean_AOD0446b', 'mean_AOD0558b', 
               'mean_AOD0672b', 'mean_AOD0866b'),

               )
              
ALIAS = dict (             Longitude = 'lon',
                            Latitude = 'lat',
                        mean_AOD1640 = 'tau1640',
                        mean_AOD1020 = 'tau1020',
                        mean_AOD0870 = 'tau870',
                        mean_AOD0670 = 'tau670',
                        mean_AOD0550 = 'tau550',
                        mean_AOD0530 = 'tau530',
                        mean_AOD0500 = 'tau500',
                        mean_AOD0490 = 'tau490',
                        mean_AOD0440 = 'tau440',
                        mean_AOD0410 = 'tau410',
                        mean_AOD0380 = 'tau380',
                        mean_AOD0340 = 'tau340',
                        mean_AOD0446b = 'tau446',
                        mean_AOD0558b = 'tau558',
                        mean_AOD0672b = 'tau672',
                        mean_AOD0866b = 'tau866',
              )

ALIAS['mean_AOD0412dpbl-l'] = 'tau412'
ALIAS['mean_AOD0470dpbl-l'] = 'tau470'
ALIAS['mean_AOD0660dpbl-l'] = 'tau660'
ALIAS['QAdpbl-l']           = 'qa_flag'

ALIAS['mean_mref0412dpbl-l'] = 'dRef412'
ALIAS['mean_mref0470dpbl-l'] = 'dRef470'
ALIAS['mean_mref0660dpbl-l'] = 'dRef660'

ALIAS['mean_surfre0412dpbl-l'] = 'sRef412'
ALIAS['mean_surfre0470dpbl-l'] = 'sRef470'
ALIAS['mean_surfre0660dpbl-l'] = 'sRef660'

ALIAS['QA-l']                = 'qa_flag'

for v in VARS['MISR_MREF']:
    if 'mean_Refl' in v:
        ALIAS[v] = 'mRef'+v[10:]

#---- 
class MAPSS(MCD43GF,IGBP):
    """Base class for Multi-sensor Aerosol Products Subset Statistics (MAPSS)"""

    def __init__ (self,Path,xVars=(),Verbose=False):
        """
        Base class for generic MAPSS dataset.
        """

        self.verb = Verbose
        self.columns = None
        self.nobs = 0
        self.Vars = VARS['META'] + xVars
        self.sample = None

        # Past is string or list
        # ----------------------
        if type(Path) is ListType:
            if len(Path) == 0:
                print("WARNING: Empty MAPSS object created")
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

        # Create timedate
        # ---------------
        self.time = []
        self.julian = ones(self.Date.size).astype('int')
        self.minutes = ones(self.Date.size).astype('int')
        for i in range(self.Date.size):
            yy, mm, dd = self.Date[i].split('-')
            h, m = self.Time[i].split(':')
            time = TIME(int(yy),int(mm),int(dd),int(h),int(m),0)
            self.time.append(time)
            self.julian[i] = time.toordinal()
            self.minutes[i] = int(h) * 60 + int(m)

        self.time = array(self.time)

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

#---
    def colocate(self,that,dt_min=120.):
        """
        Collocate another MAPSS object.
        """
        this = self
        Index = arange(this.Location.size)
        Jndex = arange(that.Location.size)
        Cndex = - ones(this.Location.size).astype('int')
        for station in self.Stations:
            if self.verb:
                print("Working on <%s>"%station)
            for julian in range(self.julian.min(),self.julian.max()+1):
                I = (this.iValid) & (this.Location==station) & (this.julian==julian) 
                J = (that.iValid) & (that.Location==station) & (that.julian==julian)
                if any(I) and any(J):
                    for i in Index[I]:
                        dt = abs(that.minutes[J]-this.minutes[i])
                        if dt.min()<=dt_min:
                            k = dt.argmin()
                            Cndex[i] = Jndex[J][k]
                            if self.verb:
                                print("- Found match on ", this.time[i], " at <%s>"%station)
        return Cndex

    collocate = colocate
    
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
        """Reads one MAPSS granule."""

        # Locate Header
        # -------------
        if self.columns == None:
            i = 0
            for line in open(filename).readlines(): 
                if line[0:9] == 'Date,Time':
                    self.columns = line[0:-1].split(',') # remove \n at end of line
                    self.skip = i + 1 # number of rows to skip for data
                i += 1


            if self.columns == None:
                raise ValueError("Cannot find Column header")

            # Read relevant columns from MAPSS granule
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
                    self.formats += ('S5',)
                elif name=='Location':
                    self.formats += ('S20',)
                else:
                    self.converters[i] = lambda s: float(s or MISSING)
                    self.formats += ('f4',)

#       Read the data
#       -------------
        print(filename)
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
    def collapse(self,I):
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

    colapse = collapse # early mispelling
    
    def savez(self,outFile):
        """
        Save all atributes to npz file.
        """
        xxx = self.converters
        self.converters = None
        savez(outFile,**self.__dict__)
        self.converters = xxx

#---
    def speciate(self,aer_x,FineMode=False):
        """
        Use GAAS to derive fractional composition.
        """

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

        if self.sample==None:
            self.sample = GFIOHandle(inFile)
        if onlyVars is None:
            onlyVars = fh.vname

        nt = self.lon.shape
        tymes = self.time
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
        albedo = zeros(len(w10m))
        albedo[:] = interp(w10m,lut.speed,lut.albedo[:,j])

        self.ocnAlbedo = albedo

#..........................................................................
class ANET(MAPSS):

    def __init__ (self, Path, xVars = VARS['ANET'], Verbose=False):
        """
        Given a directory/file list returns a MAPSS object for AERONET 
        measurements.
        """
        MAPSS.__init__(self,Path,xVars=xVars,Verbose=Verbose)

        # Interpolate to 550 nm, if missing
        # ---------------------------------
        i = (self.tau550<0) & (self.tau500>=0) & (self.tau670>=0)
        alpha = ( log(550.) - log(500.) ) / ( log(670.) - log(500.) )
        self.tau550[i] = (1.-alpha)*self.tau500[i] + alpha * self.tau670[i]

        self.iValid = self.iValid & (self.tau550>-0.01)

#..........................................................................

class DEEP_AOD(MAPSS):

    def __init__ (self, Path, xVars = VARS['DEEP_AOD'], Verbose=False):
        """
        Given a directory/file list returns a MAPSS object for AERONET 
        measurements.
        """
        MAPSS.__init__(self,Path,xVars=xVars,Verbose=Verbose)

        # Interpolate to 550 nm, if missing
        # ---------------------------------
        self.tau550 = MISSING * ones(self.tau470.shape)
        i = (self.tau470>=0) & (self.tau660>=0)
        alpha = ( log(550.) - log(470.) ) / ( log(660.) - log(470.) )
        self.tau550[i] = (1.-alpha)*self.tau470[i] + alpha * self.tau660[i]

        self.iValid = self.iValid & (self.tau550>-0.01) & (self.qa_flag>0)

#..........................................................................
class DEEP_MREF(MAPSS):

    def __init__ (self, Path, xVars = VARS['DEEP_MREF'], Verbose=False):
        """
        Given a directory/file list returns a MAPSS object for AERONET 
        measurements.
        """
        MAPSS.__init__(self,Path,xVars=xVars,Verbose=Verbose)

        self.iValid = self.iValid & (self.dRef470>=0) # & (self.qa_flag>0)

#..........................................................................
class DEEP_SREF(MAPSS):

    def __init__ (self, Path, xVars = VARS['DEEP_SREF'], Verbose=False):
        """
        Given a directory/file list returns a MAPSS object for AERONET 
        measurements.
        """
        MAPSS.__init__(self,Path,xVars=xVars,Verbose=Verbose)

        self.iValid = self.iValid & (self.sRef470>=0) # & (self.qa_flag>0)

#..........................................................................
class MISR_GEOM(MAPSS):

    def __init__ (self, Path, xVars = VARS['MISR_GEOM'], Verbose=False):
        """
        Given a directory/file list returns a MAPSS object for MISR  
        Geometry.
        """
        MAPSS.__init__(self,Path,xVars=xVars,Verbose=Verbose)

        self.iValid = (self.SolarZenith>0) # for now
        
#..........................................................................
class MISR_AOD(MAPSS):

    def __init__ (self, Path, xVars = VARS['MISR_AOD'], Verbose=False):
        """
        Given a directory/file list returns a MAPSS object for MISR AOD
        retrievals.
        """
        MAPSS.__init__(self,Path,xVars=xVars,Verbose=Verbose)

        # Interpolate to 550 nm, if missing
        # ---------------------------------
        self.iValid = self.iValid & (self.tau558>-0.01)

#..........................................................................
class MISR_MREF(MAPSS):

    def __init__ (self, Path, xVars = VARS['MISR_MREF'], Verbose=False):
        """
        Given a directory/file list returns a MAPSS object for MISR  
        Geometry.
        """
        MAPSS.__init__(self,Path,xVars=xVars,Verbose=Verbose)

        self.iValid = (self.mRef446eq1>0) # for now
        
#..........................................................................

if __name__ == "__main__":

    tdir = '/nobackup/MAPSS'
    adir = tdir+'/AERONET/AerosolOpticalDepth/2008'
    gdir = tdir+'/MISR/Geometry/2008'
    rdir = tdir+'/MISR/RegEqRefl/2008'
    mdir = tdir+'/MISR/RegBestEstimateSpectralOptDepth/2008'

    #    tdir = '/Users/adasilva/workspace/Data_Analysis/DeepBlue'

#    a = ANET(adir,Verbose=True)
#    g = MISR_GEOM(gdir,Verbose=True)
#    r = MISR_MREF(rdir,Verbose=True)
    m = MISR_AOD(mdir,Verbose=True)


def deep():
    m = MAPSS(tdir+'/Aeronet/AOD/2008')
    d = DEEP_AOD(tdir+'/DeepBlue_AOD_Land/2008',Verbose=True)
    r = DEEP_MREF(tdir+'/DeepBlue_Reflectance_Land/2008',Verbose=True)
    s = DEEP_SREF(tdir+'/DeepBlue_Surface_Reflectance_Land/2008',Verbose=True)

