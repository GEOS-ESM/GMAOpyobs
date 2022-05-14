"""
   Classes for reading PM2.5 AQS data including IMPROVE files.
"""

import os
import sys
from types    import *
from numpy    import loadtxt, ones, savez, pi, log, concatenate, arange
from datetime import timedelta,datetime as TIME

MISSING = -999.

BAD, MARGINAL, GOOD, BEST = list(range(4))

VARS = dict(
       PM25 = ('State_Code','County_Code','Site_ID','Parameter',
               'Sample_Duration','Start_Time','Unit','Date',
               'Sample_Value','Uncertainty'),
       SITES = ('State_Code','County_Code','Site_ID',
                'Latitude','Longitude','Land_Use', 'Location_Setting'),
            )
ALIAS = dict(Sample_Value = 'pm25',
             Uncertainty = 'error')
   

#---- 

class AQS(object):
    """Base class for EPA data AQS PM2.5"""

    def __init__ (self,Path,xVars=(),SitesFile=None,Verbose=False):
        """
        Instantiate a generic EPA AQS dataset. On input,

        Path       --- single file name or least of files.
        SitesFile  --- file name with des ription of AQS sites.

        """

        self.verb = Verbose
        self.columns = None
        self.nobs = 0
        self.Vars = VARS['PM25']

        # Create empty lists for SDS to be read from file
        # -----------------------------------------------
        for name in self.Vars:
            self.__dict__[name] = []

        # Read each granule, appending them to the list
        # ---------------------------------------------
        if type(Path) is ListType:
            if len(Path) == 0:
                print("WARNING: Empty AQS object created")
                return
        else:
                Path = [Path, ]
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
        self.julian = ones(self.Date.shape[0]).astype('int')
        self.minutes = ones(self.Date.shape[0]).astype('int')
        for i in range(self.Date.size):
            yy = self.Date[i][0:4]
            mm = self.Date[i][4:6]
            dd = self.Date[i][6:8]
            h, m = self.Start_Time[i].split(':')
            time = TIME(int(yy),int(mm),int(dd),int(h),int(m),0) 
            dt = timedelta(seconds=12*60*60) # add 12 hours           
            self.time.append(time+dt)
            self.julian[i] = time.toordinal()
            self.minutes[i] = int(h) * 60 + int(m)

        # Unique list of locations
        # ------------------------
        Locations = {}
        for st in self.State_Code:
            Locations[st] = 1
        self.Stations_st = list(Locations.keys())

        # If AQS Sites file specified, geo-locate observations
        # ----------------------------------------------------
        if SitesFile is not None:
           sites = AQS_SITES(SitesFile)
            


        # By default all is good if coordinates are ok
        # --------------------------------------------
#        self.iGood = (abs(self.Longitude)<=180.) & \
#                     (abs(self.Latitude)<=90.)


#---
    def _readList(self,List):
        """
        Recursively, look for files in list; list items can
        be files or directories.
        """
        for item in List:
            if os.path.isdir(item):      self._readDir(item)
            elif os.path.isfile(item):   self._readEPAdata(item)
            else:
                print("%s is not a valid file or directory, ignoring it"%item)
#---
    def _readDir(self,dir):
        """Recursively, look for files in directory."""
        for item in os.listdir(dir):
            path = dir + os.sep + item
            if os.path.isdir(path):      self._readDir(path)
            elif os.path.isfile(path):   self._readEPAdata(path)
            else:
                print("%s is not a valid file or directory, ignoring it"%item)

#---
    def _readEPAdata(self,filename):
        """Reads EPA AQS data file."""

        # Locate Header
        # -------------
        if self.columns == None:
            i = 0
            for line in open(filename).readlines(): 
                if line[0:4] == '# RD':
                    self.columns = line[0:-1].replace('\r','')\
                    .replace(' ','_')\
                    .replace('(','')\
                    .split('|')     #  split the header 
                    self.skip = i + 1 # number of rows to skip for data
                if line[0:4] == '# RC': # header different  
                    self.skip = i + 1 # number of rows to skip for data
                i += 1


            if self.columns == None:
                raise ValueError("Cannot find Column header")

            # Read relevant columns
            # ----------------------------------------
            self.iVars = ()
            self.formats = ()
            self.converters = {}
            for name in self.Vars:
                try:
                    i = self.columns.index(name)
                    
                except:
                    raise ValueError("cannot find <%s> in file"%name)
                self.iVars += (i,)
                if name=='Date':
                    self.formats += ('S8',)
                elif name=='Start_Time':
                    self.formats += ('S5',)
                elif name=='State_Code':
                    self.formats += ('S2',)
                elif name=='County_Code':
                    self.formats += ('S3',) 
                elif name=='Site_ID':
                    self.formats += ('S4',) 
                else:
                    self.converters[i] = lambda s: float(s or MISSING)
                    self.formats += ('f4',)

#       Read the data
#       -------------
       
        data = loadtxt(filename, delimiter='|',
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

        self.sitemap = {}
        self.sitemap = ((self.State_Code,self.County_Code,self.Site_ID))
        self.sitemap_= list(zip(*self.sitemap))

#---
class AQS_SITES(object):
    """Base class for EPA data AQS lat lon sites"""

    def __init__ (self,Path,xVars=(),Verbose=False):
        """
        Base class for generic EPA AQS dataset.
        """

        self.verb = Verbose
        self.columns = None
        self.nobs = 0
        self.Vars = VARS['SITES']

        # Create empty lists for SDS to be read from file
        # -----------------------------------------------
        for name in self.Vars:
            self.__dict__[name] = []

        # Read each granule, appending them to the list
        # ---------------------------------------------
        if type(Path) is ListType:
            if len(Path) == 0:
                print("WARNING: Empty MxD04_L2 object created")
                return
        else:
                Path = [Path, ]
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

#---
    def _readList(self,List):
        """
        Recursively, look for files in list; list items can
        be files or directories.
        """
        for item in List:
            if os.path.isdir(item):      self._readDir(item)
            elif os.path.isfile(item):   self._readLatLon(item)
            else:
                print("%s is not a valid file or directory, ignoring it"%item)
#---
    def _readDir(self,dir):
        """Recursively, look for files in directory."""
        for item in os.listdir(dir):
            path = dir + os.sep + item
            if os.path.isdir(path):      self._readDir(path)
            elif os.path.isfile(path):   self._readLatLon(path)
            else:
                print("%s is not a valid file or directory, ignoring it"%item)

#---
    def _readLatLon(self,filename):
        """Reads EPA AQS data file."""
        # Locate Header
        # -------------
        if self.columns == None:
            i = 0
            for line in open(filename).readlines(): 
                if line[0:5] == 'State':
                    self.columns = line[0:-1].replace('\r','')\
                    .replace(' ','_')\
                    .split('|')     #  split the header 
                    self.skip = i + 1 # number of rows to skip for data
                i += 1


            if self.columns == None:
                raise ValueError("Cannot find Column header")

            # Read relevant columns
            # ----------------------------------------
            self.iVars = ()
            self.formats = ()
            self.converters = {}
            for name in self.Vars:
                try:
                    i = self.columns.index(name)
                    
                except:
                    raise ValueError("cannot find <%s> in file"%name)
                self.iVars += (i,)
               
                
                if name=='Location_Setting':
                    self.formats += ('S21',)
                elif name=='Land_Use':
                    self.formats += ('S21',)
                elif name=='State_Code':
                    self.formats += ('S2',)
                elif name=='County_Code':
                    self.formats += ('S3',) 
                elif name=='Site_ID':
                    self.formats += ('S4',) 
                else:
                    self.converters[i] = lambda s: float(s or MISSING)
                    self.formats += ('f4',)

#       Read the data
#       -------------
       
        data = loadtxt(filename, delimiter='|',
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


#----
def _site(self,filename):

        self.sitemap = {}
        self.sitemap = ((self.State_Code,self.County_Code,self.Site_ID))
        self.sitemap_= list(zip(*self.sitemap))


#..........................................................................

if __name__ == "__main__":

    m = IMPROVE('/nobackup/5/AQS/pm2.5_local/Y2011/RD_501_88101_2011-0.txt')
    g = SITE_MAP('/nobackup/5/SITE_MON_5_13/SITE_MON_5_13_AA_select.txt')
    SITE_MAP._site(g,'/nobackup/5/SITE_MON_5_13/SITE_MON_5_13_AA_select.txt')

   
