"""
   Simple script to return the names and locations of AERONET sites
   that have reported data in the last N days.

   P. Castellanos June 2017
"""

import os
import sys
from   datetime import datetime, timedelta
from   dateutil.parser import parse as isoparse
from   glob     import glob
import argparse
import numpy as np

VARS = ( 'AERONET_Site',
         'Longitude',
         'Latitude',
         )

ALIAS = { 'Longitude': 'lon',
          'Latitude' : 'lat',
          'AERONET_Site'  : 'Location',
        }
MISSING = -999
#---- 
class AERONET(object):

    def __init__ (self,Path,Vars=VARS,Verbose=False):
        """
        Class for generic AERONET dataset.
        """

        self.verb = Verbose
        self.Vars = Vars
        self.columns = None
        self.nobs    = 0

        # Past is string or list
        # ----------------------
        if type(Path) is list:
            if len(Path) == 0:
                print "WARNING: NO AERONET FILES FOUND"
                return
        else:
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
                self.__dict__[var] = np.concatenate(self.__dict__[var])
            except:
                print "Failed concatenating "+var

        # Make aliases
        # ------------
        Alias = ALIAS.keys()
        for var in self.Vars:
            if var in Alias:
                self.__dict__[ALIAS[var]] = self.__dict__[var]

        # Unique list of locations
        # ------------------------
        Locations = {}
        for loc,lon,lat in zip(self.Location,self.lon,self.lat):
            Locations[loc] = (lon,lat) 
        self.Stations = Locations.keys()


        # write aeronet_stations.rc
        # -------------------------
        f = open('aeronet_stations.rc','w')
        f.write('name,lon,lat\n')
        for loc in Locations:
          lon, lat = Locations[loc]
          f.write('{},{},{}\n'.format(loc,lon,lat))

        f.close()


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
                print "%s is not a valid file or directory, ignoring it"%item
#---
    def _readDir(self,dir):
        """Recursively, look for files in directory."""
        for item in os.listdir(dir):
            path = dir + os.sep + item
            if os.path.isdir(path):      self._readDir(path)
            elif os.path.isfile(path):   self._readGranule(path)
            else:
                print "%s is not a valid file or directory, ignoring it"%item

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
                                          .replace('Site_','')\
                                          .split('(')[0],]
                    break
                i += 1

            if self.columns == None:
                raise ValueError, "Cannot find Column header"

            # Read relevant columns from AERONET granule
            # ----------------------------------------
            self.iVars = ()
            self.formats = ()
            self.converters = {}
            for name in self.Vars:
                try:
                    i = self.columns.index(name)
                except:
                    raise ValueError, "cannot find <%s> in file <%s>"%(name,filename)
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
        data = np.loadtxt(filename, delimiter=',',
                       dtype={'names':self.Vars,'formats':self.formats},
                       converters = self.converters,
                       skiprows=self.skip, usecols=self.iVars)
        N = len(data)

        self.nobs += N

        # Save data columns as attributes, with nicer names
        # -------------------------------------------------
        for i in range(len(self.Vars)):
            if self.formats[i]=='f4':
                v = np.ones(N)
                for j in range(N):
                    v[j] = data[j][i]
            else:
                v = []
                for j in range(N):
                    v.append(data[j][i])

            self.__dict__[self.Vars[i]] += [v,]

#---
def _convert2Float(s,missing=MISSING):
   if s in ('N/A', 'NA','','NaN',None):
      result = missing
   else:
      result = float(s or missing)

   return result

#.......................................................................................
def retrieve( filename='aeronet.csv',
              request='year=2016&month=9&day=1&year2=2016&month2=9&day2=30&AOD15=1&AVG=10',
              verbose=True):
    """
    Use wget to retrieve aeronet data using the Version 3 webservice.'
    """
    webportal = 'https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3?'
    cmd = 'wget --no-check-certificate -q -O %s "%s%s&if_no_html=1"'%(filename,webportal,request)

    if verbose:
        print cmd
        
    if os.system(cmd):
        raise ValueError, "Cannot retrieve request <%cmd>"

        
        
#..........................................................................

if __name__ == "__main__":

    # Defaults
    nDays = 30

#   Parse command line options
#   --------------------------

    parser = argparse.ArgumentParser()
    parser.add_argument("sdate",
                        help="search date (YYYY-MM-DD)")

    parser.add_argument("-n","--nDays", default=nDays, type=int,
                        help="number of days before search date to search for AERONET data (default=%i)"%nDays)  


    args = parser.parse_args()

    edate = isoparse(args.sdate)
    DT    = timedelta(days=args.nDays)
    sdate = edate - DT

    srequest = 'year={}&month={}&day={}&'.format(sdate.year,sdate.month,sdate.day)
    erequest = 'year2={}&month2={}&day2={}&'.format(edate.year,edate.month,edate.day)
    request  = srequest + erequest + 'AOD15=1&AVG=10'
    
    retrieve(request=request)

    a = AERONET('aeronet.csv')

    os.remove('aeronet.csv')
