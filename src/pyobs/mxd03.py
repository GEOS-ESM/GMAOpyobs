"""
Reads Level 1 MOD03/MYD03 granules, adds

This software is hereby placed in the public domain.
Arlindo.daSilva@nasa.gov

"""

import os
import sys
from types    import *

from numpy    import zeros, ones, concatenate, array, shape, arange, tile 

from datetime import date, datetime, timedelta
from glob     import glob

# from gfio      import GFIO

from pyhdf.SD import SD, HDF4Error

#---  

DATE_START = datetime(1993,1,1,0,0,0)

ALIAS = dict ( Longitude = 'lon', Latitude = 'lat', Height = 'ZS' )

MISSING = 999.999

#...........................................................................
class MxD03Handle(object):
    """
    Generic container for MxD03.
    """
    def __init__(self,name):
        self.name = name

class MxD03(object):
    """
    This class implements the MODIS Level 1 MOD03 interface.
    """

    def __init__ (self,Path,Verb=0):
       """
       Reads individual granules or a full day of Level 1 MOD03/MYD03 files
       present on a given *Path* and returns a single object with
       all data concatenated for a given algorithm. On input, 

       Required parameters:
         Path -- can be a single file, a single directory, of a list
                 of files and directories.  Directories are
                 transversed recursively. If a non MOD03/MYD03 Level 1
                 file is encountered, it is simply ignored.
       Optional parameters:
         Verb      -- Verbose level:
                 0 - really quiet (default)
                 1 - Warns if invalid file is found
                 2 - Prints out non-zero number of aerosols in each file.

       """

#      Initially are lists of numpy arrays for each granule
#      ------------------------------------------------
       self.verb = Verb
       self.sat  = None # Satellite name
       self.col  = None # collection, e.g., 005
       self.sample = None # will hold sampled variables later
       self.ica = None # will hold ICA indices later
       self.SDS = ('Longitude', 'Latitude', 'EV start time', 'Height')
       self.path = Path

       # Create empty lists for SDS to be read from file
       # -----------------------------------------------
       for sds in self.SDS:
           sds_ = sds.replace(' ','_')
           self.__dict__[sds_] = []

       # Read each granule, appending them to the list
       # ---------------------------------------------
       if type(Path) is ListType:
           if len(Path) == 0:
               self.nobs = 0
               print("WARNING: Empty Mxd03 object created")
               return
       else:
           Path = [Path, ]
       self._readList(Path)

       # Make each attribute a single numpy array
       # ----------------------------------------
       for sds in self.SDS:
           sds_ = sds.replace(' ','_')
           try:
               self.__dict__[sds_] = concatenate(self.__dict__[sds_])
           except:
               print("Failed concatenating "+sds)

       # Make aliases for convenience
       # ----------------------------
       Alias = list(ALIAS.keys())
       for sds in self.SDS:
           sds_ = sds.replace(' ','_')
           if sds_ in Alias:
               self.__dict__[ALIAS[sds_]] = self.__dict__[sds_] 

       # Expand scan time
       # ----------------
       t = ones(self.lon.shape[0])
       for i in range(self.EV_start_time.size):
           t[10*i:10*i+10] = self.EV_start_time[i]
           
       # Create corresponding python time
       # --------------------------------
       self.time = array([DATE_START+timedelta(seconds=int(s)) for s in t])
       self.gtime = DATE_START+timedelta(seconds=int(self.EV_start_time.mean())) # mean granule time

#---
    def linearSampleFile(self, inFile, onlyVars=None):
        """
        Interpolates some or all variables of inFile.
        """
        from gfio import GFIOctl

        # Open file
        # ---------
        f = GFIOctl(inFile)
        if self.sample == None:
            self.sample = MxD03Handle(inFile)
        if onlyVars is None:
            onlyVars = f.vname

        # Flatten coordinates
        # -------------------
        nt, nr = self.lon.shape
        tyme = tile(self.time,(nr,1)).T
        lons = self.lon.ravel()
        lats = self.lat.ravel()
        tymes = tyme.ravel()
        tymes[:] = self.gtime # use mean granule time
             
        # Loop over variables on file
        # ---------------------------
        for v in onlyVars:
            if self.verb:
                print(" <> Linear sampling ", v)
            var = f.sample(v,lons,lats,tymes,Verbose=self.verb)
            if len(var.shape) == 1:
                self.sample.__dict__[v] = var.reshape((nt,nr))
            elif len(var.shape) == 2:
                var = var.T # shape should be (nobs,nz)
                self.sample.__dict__[v] = var.reshape((nt,nr,-1))
            else:
                raise IndexError('variable <%s> has rannk = %d'%len(var.shape))

#---
    def nearestSampleFile(self, inFile, onlyVars=None):
        """
        Sample some or all variables of inFile. This simply returns the
        values of the gridbox the pixel fall into (nearest neighbor).
        """
        from gfio import GFIOctl

        # Instiate grads and open file
        # ----------------------------
        f = GFIOctl(inFile)
        if self.sample == None:
            self.sample = MxD03Handle(inFile)
        if onlyVars is None:
            onlyVars = f.vname

        # Flatten coordinates
        # -------------------
        nt, nr = self.lon.shape
        tyme = tile(self.time,(nr,1)).T
        lons = self.lon.ravel()
        lats = self.lat.ravel()
        tymes = tyme.ravel()
        tymes[:] = self.gtime # use mean granule time
             
        # Loop over variables on file
        # ---------------------------
        for v in onlyVars:
            if self.verb:
                print(" <> Nearest sampling ", v)
            var = f.sample(v,lons,lats,tymes,
                           algorithm='nearest',Verbose=self.verb)
            self.var = var
            if len(var.shape) == 1:
                self.sample.__dict__[v] = var.reshape((nt,nr))
            elif len(var.shape) == 2:
                var = var.T # shape should be (nobs,nz)
                self.sample.__dict__[v] = var.reshape((nt,nr,-1))
            else:
                raise IndexError('variable <%s> has rannk = %d'%len(var.shape))

#---
    def getICAindx(self,inFile):
        """
        Get gridded related indices for ICA calculations.  It calculates
        grid-box indices useful to implement Independent COlumn
        Approximation (ICA) type of algorithms.
        """
        from gfio import GFIOctl

        # Instiate grads and open file
        # ----------------------------
        f = GFIOctl(inFile)
        if self.ica == None:
            self.ica = MxD03Handle(inFile)

        # Handle grid indices for ICA algorithms
        # --------------------------------------
        nt, nr = self.lon.shape
        if self.verb:
            print(" <> Performing ICA index generation")
        iCoord, jCoord = f.coordNN(self.lon.ravel(),self.lat.ravel())
        iS, jS = iCoord.astype('S5'), jCoord.astype('S5'), 
        keys = [ ii+','+jj for ii,jj in zip(iS,jS) ]
        Indices = dict()
        for n in range(iCoord.size):
            Indices[keys[n]] = (iCoord[n],jCoord[n])

        # save his for later
        # ------------------
        self.ica.iCoord = iCoord.reshape((nt,nr))
        self.ica.jCoord = jCoord.reshape((nt,nr))
        self.ica.Indices = Indices
        
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
        """Reads one MOD03/MYD03 granule."""

        # Don't fuss if the file cannot be opened
        # ---------------------------------------
        try:
            if self.verb:
                print("[] Working on "+filename)
            hfile = SD(filename)
        except HDF4Error:
            if self.verb > 2:
                print("- %s: not recognized as an HDF file"%filename)
            return 

        # Read select variables (reshape to allow concatenation later)
        # ------------------------------------------------------------
        for sds in self.SDS:
            sds_ = sds.replace(' ','_')
            v = hfile.select(sds).get()
            #a = hfile.select(sds).attributes()
            #if a['scale_factor']!=1.0 or a['add_offset']!=0.0:
            #    v = a['scale_factor'] * v + a['add_offset']
            self.__dict__[sds_].append(v) 

        # Core Metadata
        # -------------
        cm = hfile.attributes()['CoreMetadata.0']

#       Satellite name
#       --------------
        if self.sat is None:
            sat = cm.split('ASSOCIATEDPLATFORMSHORTNAME')[1].split('\n')[3].split('=')[1]
            self.sat = sat.lstrip().replace('"','')
            
#       Collection
#       ----------
        if self.col is None:
            col = int(cm.split('COLLECTION')[1].split('VERSIONID')[1].split('\n')[2].split('=')[1])
            self.col = "%03d"%col

#............................................................................

if __name__ == "__main__":

    mod03 = '/nobackup/MODIS/Level1/MYD03/sample/MYD03.A2012228.1200.005.2012229145830.hdf'
    asm_Nx = '/nobackup/fp/opendap/asm_Nx.ctl'
    asm_Nv = '/nobackup/fp/opendap/asm_Nv.ctl'

    m = MxD03(mod03, Verb=1)
    m.getICAindx(asm_Nx)
    
def hold():

    m.linearSampleFile(asm_Nx,  onlyVars=('PS', 'T2M', 'SLP', 'U10M', 'V10M'))
    m.nearestSampleFile(asm_Nx, onlyVars=('QV2M',))

    m.linearSampleFile(asm_Nv, onlyVars=('T', 'DELP'))
    m.nearestSampleFile(asm_Nv, onlyVars=('O3','RH'))
    
                        
