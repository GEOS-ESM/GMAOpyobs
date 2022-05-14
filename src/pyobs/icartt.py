"""
   Class for reading ICARTT ASCII files.
"""

from datetime import datetime, timedelta
from numpy    import loadtxt, ones, NaN, concatenate, array, pi, cos, sin, arccos, zeros
from glob     import glob
import gzip
import collections

from .config import Config

MISSING = -99999.0
ULOD    = -77777.0 # Upper Limit of Detection (LOD) flag
LLOD    = -88888.0 # Lower Limit of Detection (LOD) flag

class ICARTT(object):
    """Reads ICARTT text files into Numpy arrays"""

    def __init__ (self,Filenames,Alias=None,FixUnits=True,Verbose=False):
        """
        Loads one or more ICART text files, creating an ICART object.
        When entering many files, make sure they are in chronological
        order.
        
        FixUnits --- Convert some enginering units to MKS
        """

        self.verb = Verbose

        # Many files
        # ----------
        if (type(Filenames) is list) or type(Filenames) is tuple:
            self._readManyFiles (Filenames,Alias)

        # Or glob
        # -------
        else:
            if Filenames.count('*')>0 or \
               Filenames.count('[')>0 or \
               Filenames.count('?')>0:
                   
               Filenames_ = glob(Filenames)
               Filenames_.sort()
               self._readManyFiles (Filenames_,Alias)

            else:
                self._readOneFile (Filenames,Alias)

        # If variable name in spreadsheet does not match names in units...
        # This is a kludge to cope with weird DC8 NAV file in ARCTAS
        # ----------------------------------------------------------------
        for var in self.vIndex:
            if var not in self.Vars:
                self.__dict__[var] = self.__dict__[self.Vars[self.vIndex[var]]]
                
        # Fix Units
        # ---------
        if FixUnits:
            for var in self.Units:
                if self.Units[var].upper() == "FEET":
                    self.Units[var] = 'm'
                    self.__dict__[var] = 0.3048 * self.__dict__[var]
                if self.Units[var].upper() == "KILOMETER" or \
                   self.Units[var].upper() == "KILOMETERS" or \
                   self.Units[var].upper() == "KM":
                    self.Units[var] = 'm'
                    self.__dict__[var] = self.__dict__[var] / 1000.
                if self.Units[var].upper() == "KNOT" or \
                   self.Units[var].upper() == "KNOTS":
                    self.Units[var] = 'm s-1'
                    self.__dict__[var] = 0.514444444444444 * self.__dict__[var] 
                if self.Units[var].upper() == "C" or\
                   self.Units[var].upper() == "CELSIUS": 
                    self.Units[var] = 'K'
                    self.__dict__[var] = self.__dict__[var] + 273.15 
                if self.Units[var].upper() == "F" or\
                   self.Units[var].upper() == "FAHRENHEIT":
                    self.Units[var] = 'K'
                    self.__dict__[var] = self.__dict__[var]/1.8 + 255.372

        # Standardize Navigation information
        # Supported Aircrafts: UMd, P3, UC12, NCAR C130
        # ---------------------------------------------
        self.Nav = dict ( Time=self.tyme, Longitude=None, Latitude=None, Altitude=None, Pressure=None )
        for var in self.Vars:
            VAR = var.upper()
            if VAR in ('LONGITUDE', 'LONGITUDE_YANG', 'LONGITUDE_DEG','FMS_LON', 'GPS_LON', 'LON', 'GGLON'):
                self.Nav['Longitude'] = self.__dict__[var]
            if VAR in ('LATITUDE', 'LATITUDE_YANG', 'LATITUDE_DEG','FMS_LAT', 'GPS_LAT', 'LAT', 'GGLAT' ):
                self.Nav['Latitude'] = self.__dict__[var]
            if VAR in ('GPSALT', 'MSL_GPS_ALTITUDE_YANG', 'GPSALT_M', 'FMS_ALT_PRES', 'GPS_ALT', 'GGALT' ):
                self.Nav['Altitude'] = self.__dict__[var]
            if VAR in ('PRESSURE', 'PRESSURE_YANG', 'C_STATICPRESSURE', 'PSXC'):
                self.Nav['Pressure'] = self.__dict__[var]

        # Navigation shorthands
        # ---------------------
        self._shorthands()
                
#--
    def _readManyFiles (self,Filenames,Alias=None):
        """
        Loads multiple ICART text files. It is assumed that these are
        conformant files with same variables, etc., only the times
        and data values are different.
        """

        # Read first file
        # ---------------
        if self.verb:
            print("[] Processing ", Filenames[0])
        self._readOneFile(Filenames[0],Alias)

        AllVars = self.Vars + ['tyme',]

        # Container for holding data in each file
        # ---------------------------------------
        VARS = {}
        for var in AllVars:
            VARS[var] = [self.__dict__[var],]
            
        # Loop over remaining files
        # -------------------------
        for filename in Filenames[1:]:

            if self.verb:
                print("[] Processing ", filename)
                
            # Read this file
            # --------------
            self._readOneFile(filename,Alias)

            # Save variables for final merge
            # ------------------------------
            for var in AllVars:
                VARS[var] += [self.__dict__[var],]

        # Concatenate variables
        # ---------------------
        for var in AllVars:
            self.__dict__[var] = concatenate(VARS[var])
            
#--
    def _readOneFile (self,filename,Alias=None):
        """
        Loads a single ICART text file.
        """

        # Top header
        # ----------
        if filename[-3:] == '.gz':
            f = gzip.open(filename, mode='r')
        else:
            f = open(filename,mode='r')
        header = f.readline().replace('\n','').replace('\r','')
        if ',' in header:
            tokens = header.split(',')
            delim = ','
        else:
            tokens = header.split()
            delim = None
        self.n_header, self.FFI = int(tokens[0]), int(tokens[1])
        if self.FFI != 1001:
            raise ValueError("Only ICARTT File Format Index 1001 currently support; got %d"%self.FFI)

        self.PI = f.readline().replace('\n','').replace('\r','')
        self.INSTITUTION = f.readline().replace('\n','').replace('\r','')
        self.PRODUCT = f.readline().replace('\n','').replace('\r','')
        self.CAMPAIGN = f.readline().replace('\n','').replace('\r','')

        # Date range on file
        # ------------------
        f.readline() # ignore file volume number
        tokens = f.readline().replace('\n','').replace('\r','')
        if ',' in tokens:
            tokens = tokens.replace(' ','').split(',')  
        else:
            tokens = tokens.split()
        t = [ int(i) for i in tokens ]
        self.beg_date = datetime(t[0],t[1],t[2])
        self.rev_date = datetime(t[3],t[4],t[5])

        # Units
        # -----
        f.readline() # skip this for now
        f.readline() # skip this for now
        self.nVars = int(f.readline().replace('\n','').replace('\r',''))
        tokens = f.readline().replace(' ','').replace('\n','').replace('\r','').split(',')
        self.scale_factors = [1.,] + [ float(x) for x in tokens ]
        tokens = f.readline().replace('\n','').replace('\r','')
        if ',' in tokens:
            tokens = tokens.replace(' ','').split(',')  
        else:
            tokens = tokens.split()
        self.missing_values = [MISSING,] + [ float(x) for x in tokens ]
        self.Units = {}
        self.Long = {}
        self.vIndex = {}
        for i in range(self.nVars):
            tokens = f.readline().replace('\n','').replace('\r','')
            if ',' in tokens:
                tokens = tokens.replace(' ','').split(',')  
            else:
                tokens = tokens.split(' ')
            name, units = tokens[0], tokens[1]
            name = name.replace('.','')
            if len(tokens)>2:
                longn = tokens[2]
            else:
                longn = name
            self.Units[name] = units
            self.Long[name] = longn
            self.vIndex[name] = i+1 # because UTC is first
        self.nVars += 1 # to include Start_UTC
            
        # Variable Names
        # --------------
        for i in range(self.n_header-13-self.nVars+1): f.readline()
        if delim == ',':
            self.Vars = f.readline().replace(' ','').replace('\n','').replace('\r','').replace('.','').split(',')
        else:
            self.Vars = f.readline().replace('\n','').replace('\r','').replace('.','').split(delim)
        f.close()

        # Remove duplicates
        # -----------------
        for dup in [item for item, count in list(collections.Counter(self.Vars).items()) if count > 1]:
            self.Vars.remove(dup)
        
#       Use Config to load other attributes
#       -----------------------------------
        cf = config.Config(filename)
        for rc in list(cf.keys()):
            self.__dict__[rc.upper()] = cf(rc)
        
        # Read relevant columns from MAPSS granule
        # ----------------------------------------
        formats = ()
        converters = {}
        i = 0
        for name in self.Vars:
            converters[i] = lambda s: float(s or MISSING)
            formats += ('f4',)
            i += 1

        # Read the data
        # -------------
        data = loadtxt(filename, delimiter=delim,
                       dtype={'names':self.Vars,'formats':formats},
                       converters = converters,
                       skiprows=self.n_header)
 
        try:
            N = len(data)
        except:
            N = 0 # should be 1, but don't bother with these
        
#       Save data columns as attributes
#       -------------------------------
        try:
            ULOD = float(self.ULOD_FLAG)
        except:
            ULOD = -77777.0
        try:
            LLOD = float(self.LLOD_FLAG)
        except:
            LLOD = -88888.0
        BAD_ = -980. # Hack for inconsistent files
        for i in range(len(self.Vars)):
            MISSING_ = self.missing_values[i]
            v = ones(N)
            for j in range(N):
                v[j] = data[j][i]
            bad = (v==MISSING_)|(v==LLOD)|(v==ULOD)|(v<=BAD_)
            v[bad] = NaN              # use NaN for bad data
            self.__dict__[self.Vars[i]] = v

        # For some merged files there is no mid time, so create one
        # ---------------------------------------------------------
        if (self.Vars[0].upper()=='TIME_START') and (self.Vars[1].upper()=='TIME_STOP'):
            self.Time_Mid = (self.__dict__[self.Vars[0]]+self.__dict__[self.Vars[1]])/2.
            self.Vars = self.Vars + ['Time_Mid',]

        # Find time variables
        # -------------------
        Tvar = []
        for var in self.Vars:
            if var.upper().count('UTC_')>0   or \
               var.upper().count('_UTC')>0   or \
               var.upper() == 'GPS_TIME'     or \
               var.upper() == 'TIME_START'   or \
               var.upper() == 'TIME_STOP'    or \
               var.upper() == 'TIME_MID'     or \
               var.upper() == 'UTC':
               Tvar += [var,]

        if len(Tvar)==0:
            Tvar = [self.Vars[0],]   # safe guard

        # Express time variables in python native format
        # ----------------------------------------------
        T0 = self.beg_date
        for tvar in Tvar:
            self.__dict__[tvar] = array([ T0+timedelta(seconds=t) for t in self.__dict__[tvar] ])

            
        # Find representative time variable
        # ---------------------------------
        tvar = self.Vars[0]
        for tvar_ in ( 'MID_UTC', 'UTC_MID', 'TIME_MID'):
            for var in self.Vars:
                if var.upper() == tvar_: # case insensitive
                    tvar = var
                    
        self.tyme = self.__dict__[tvar] # just an alias

#--
    def _shorthands(self):
        """
        Add navigation short hands: lon, lat, etc.
        """
        self.lon = self.Nav['Longitude']
        self.lat = self.Nav['Latitude']
        self.prs = self.Nav['Pressure']
        self.alt = self.Nav['Altitude']

        #--
    def fixNav(self, nav):
        """
        Given a possibly longer ICARTT object *nav* with (Longitude,Latitude,Altitude)
        properly defined, complete the navigation of this object. This function is needed
        because some instrument teams do not include geo-location on their ICARTT files.
        """

        n_this = len(self.tyme)
        n_nav = len(nav.tyme)

        # Consistency tests
        # -----------------
        if n_this > n_nav:
            raise ValueError('Navigation object is not long enough.')
        if self.tyme[0] < nav.tyme[0]:
            raise ValueError('Navigation object start after this object')
        if self.tyme[-1] > nav.tyme[-1]:
            raise ValueError('Navigation object end before this object')
        
        # Compute proper time offsets
        # ---------------------------
        i_off = int((self.tyme[0]-nav.tyme[0]).total_seconds())
        j_off = int((nav.tyme[-1]-self.tyme[-1]).total_seconds())

        # More consistency tests
        # ----------------------
        if self.tyme[0] != nav.tyme[i_off]:
            raise ValueError('inconsistent start times')
        if self.tyme[-1] != nav.tyme[n_nav-j_off-1]:
            raise ValueError('inconsistent end times')
        if (n_nav-j_off-i_off) != n_this:
            raise ValueError('inconsistent sizes')
        
        # Complete navigation of this object where needed
        # -----------------------------------------------
        for coord in ( 'Longitude', 'Latitude', 'Altitude', 'Pressure' ):
            if self.Nav[coord] is None:
                self.Nav[coord] = nav.Nav[coord][i_off:n_nav-j_off]

        # Navigation shorthands
        # ---------------------
        self._shorthands()
        
#--
    def getSpeed(self, metric=True,skip=60):
        """
        Return plane speed in km/h or nm/h depending on metric)
        """
        # coarsen trajectory
        # ------------------
        lon = self.Nav['Longitude'][::skip]
        lat = self.Nav['Latitude'][::skip]
        t = self.tyme[::skip]
        
        dst = _getDist(lon,lat) # meters
        if metric:
            dst = dst / 1000.               # km
        else:
            dst = 0.000539957 * dst # nautical miles

        # Fix units
        # ---------
        dt = t[1:]-t[:-1]
        dt_hour = array([dt_.total_seconds()/3600. for dt_ in dt]) # hours
        
        spd = dst / dt_hour

        return spd
        
def _gcDist(x1,y1,z1,x2,y2,z2):
    """
    Return great circle distance.
    """
    a = (6378137.00+6356752.3142)/2. # mean earth radius
    cosa = x1*x2 + y1*y2 + z1*z2
    cosa[cosa>1.] = 1.
    return a * arccos(cosa)
    
def _getDist(lon,lat):
    """
    Return distance along trajectory.
    """
    d2r = pi / 180.
    rlon, rlat = (d2r * lon, d2r * lat)
    x, y, z = (cos(rlat)*cos(rlon),cos(rlat)*sin(rlon),sin(rlat))
    dist = _gcDist(x[:-1],y[:-1],z[:-1],
                   x[1:], y[1:], z[1:])
    return dist

        
#....................................................................

if __name__ == "__main__":

    x = ICARTT('SEAC4RS-mrg60-dc8_merge_20130923_R7.ict')

