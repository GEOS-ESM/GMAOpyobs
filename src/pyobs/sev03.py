"""
Reads Level 1 SEVIRI pixel location/geometry and surface height files
"""

import os
import sys
from types    import *

from pyhdf.SD import SD, HDF4Error

import numpy as np

#---  

ALIAS = dict ( Longitude = 'lon', Latitude = 'lat', P001_MSG_GLOBE_DEM = 'ZS' )

OFFVIEW = -999
OCEAN = -500

#...........................................................................
class SEV03Handle(object):
    """
    Generic container for SEV03.
    """
    def __init__(self,name):
        self.name = name

class SEV03(object):
    """
    This class implements the SEVIRI Level 1 SEV03 interface.
    """

    def __init__ (self,fgeom,fhgt,Verb=0):
       """
       Reads Level 1 SEVIRI pixel location/geometry and surface height files.
       On input, 

       Required parameters:
         fgeom -- filename of pixel location/geometry file
         fhgt  -- filename of pixel surface height file

       Optional parameters:
         Verb      -- Verbose level:
                 0 - really quiet (default)
                 1 - Warns if invalid file is found
                 2 - Prints out non-zero number of aerosols in each file.
       """

#      Initially are lists of numpy arrays for each granule
#      ------------------------------------------------
       self.verb	= Verb
       self.sat		= 'MSG+0000' # Satellite name
       self.sample	= None # will hold sampled variables later
       self.ica		= None # will hold ICA indices later
       self.SDS_geom	= ('Longitude', 'Latitude')
       self.SDS_hgt	= ('P001_MSG_GLOBE_DEM',)
       self.fgeom	= fgeom
       self.fhgt	= fhgt

       # Read each file
       # --------------
       self._read(fgeom,fhgt)

       # Make aliases for convenience
       # ----------------------------
       for sds in self.SDS_geom + self.SDS_hgt:
           sds_ = sds.replace(' ','_')
           if sds_ in ALIAS:
               self.__dict__[ALIAS[sds_]] = self.__dict__[sds_] 

#---
    def linearSampleFile(self, inFile, tyme, onlyVars=None):
        """
        Interpolates some or all variables of inFile.
        """
        from gfio import GFIOctl

        # constant time
        tymes = np.repeat(tyme,self.nobs)
             
        # Open file
        # ---------
        f = GFIOctl(inFile)
        if self.sample == None:
            self.sample = SEV03Handle(inFile)
        if onlyVars is None:
            onlyVars = f.vname

        # Loop over variables on file
        # ---------------------------
        for v in onlyVars:
            if self.verb:
                print(" <> Linear sampling ", v)
            var = f.sample(v,self.lon,self.lat,tymes,Verbose=self.verb)
            if var.ndim == 1:
                self.sample.__dict__[v] = var
            elif var.ndim == 2:
                self.sample.__dict__[v] = var.T   # -> (nobs,nz)
            else:
                raise IndexError('variable <%s> has rank = %d' % var.ndim)

#---
    def nearestSampleFile(self, inFile, tyme, onlyVars=None):
        """
        Sample some or all variables of inFile. This simply returns the
        values of the gridbox the pixel fall into (nearest neighbor).
        """
        from gfio import GFIOctl

        # constant time
        tymes = np.repeat(tyme,self.nobs)
             
        # Instiate grads and open file
        # ----------------------------
        f = GFIOctl(inFile)
        if self.sample == None:
            self.sample = SEV03Handle(inFile)
        if onlyVars is None:
            onlyVars = f.vname

        # Loop over variables on file
        # ---------------------------
        for v in onlyVars:
            if self.verb:
                print(" <> Nearest sampling ", v)
            var = f.sample(v,self.lon,self.lat,tymes,
                           algorithm='nearest',Verbose=self.verb)
            if var.ndim == 1:
                self.sample.__dict__[v] = var
            elif var.ndim == 2:
                self.sample.__dict__[v] = var.T   # -> (nobs,nz)
            else:
                raise IndexError('variable <%s> has rank = %d' % var.ndim)

#---
    def getICAindx(self,inFile):
        """
        Get grid related indices for ICA calculations.  It calculates
        grid-box indices useful to implement Independent Column
        Approximation (ICA) type of algorithms.
        """
        from gfio import GFIOctl

        # Instiate grads and open file
        # ----------------------------
        f = GFIOctl(inFile)
        if self.ica == None:
            self.ica = SEV03Handle(inFile)

        # Handle grid indices for ICA algorithms
        # --------------------------------------
        if self.verb:
            print(" <> Performing ICA index generation")
        # GEOS-5 coordinates of each pixel
        iCoord, jCoord = f.coordNN(self.lon,self.lat)
        # dictionary containing pixels for each gridcolumn
        Indices = {}
        for n in range(iCoord.size):
          Indices.setdefault((iCoord[n],jCoord[n]),[]).append(n)

        # save for later
        # --------------
        self.ica.iCoord = iCoord
        self.ica.jCoord = jCoord
        self.ica.Indices = Indices
        
#---
    def _read(self,fgeom,fhgt):
        """Reads location/geometry and height files"""

        try:
            if self.verb:
                print("[] Working on "+fgeom)
            hfile = SD(fgeom)
        except HDF4Error:
            if self.verb > 2:
                print("- %s: not recognized as an HDF file"%fgeom)
            return 
        for sds in self.SDS_geom:
            sds_ = sds.replace(' ','_')
#T          self.__dict__[sds_] = hfile.select(sds).get()
#           self.__dict__[sds_] = hfile.select(sds).get()[:1856,:1856] #00 #T too big as is
#           self.__dict__[sds_] = hfile.select(sds).get()[:1856,1856:] #01 #T too big as is
#           self.__dict__[sds_] = hfile.select(sds).get()[1856:,:1856] #10 #T too big as is
            self.__dict__[sds_] = hfile.select(sds).get()[1856:,1856:] #11 #T too big as is
        hfile.end()

        try:
            if self.verb:
                print("[] Working on "+fhgt)
            hfile = SD(fhgt)
        except HDF4Error:
            if self.verb > 2:
                print("- %s: not recognized as an HDF file"%fhgt)
            return 
        for sds in self.SDS_hgt:
            sds_ = sds.replace(' ','_')
#T          self.__dict__[sds_] = hfile.select(sds).get()
#           self.__dict__[sds_] = hfile.select(sds).get()[:1856,:1856] #00 #T too big as is
#           self.__dict__[sds_] = hfile.select(sds).get()[:1856,1856:] #01 #T too big as is
#           self.__dict__[sds_] = hfile.select(sds).get()[1856:,:1856] #10 #T too big as is
            self.__dict__[sds_] = hfile.select(sds).get()[1856:,1856:] #11 #T too big as is
        hfile.end()

        # convert "ocean height" to zero
        self.P001_MSG_GLOBE_DEM[self.P001_MSG_GLOBE_DEM == OCEAN] = 0

        # store mask of off-view points
        self.offview = np.logical_or(
          self.Longitude == OFFVIEW, 
          self.Latitude  == OFFVIEW)
        self.offview = np.logical_or(self.offview, 
          self.P001_MSG_GLOBE_DEM == OFFVIEW)

        # store unflattened shape
        self.orgshape = self.offview.shape

        # store flattened, in-view SDS
        for sds in self.SDS_geom + self.SDS_hgt:
            sds_ = sds.replace(' ','_')
            self.__dict__[sds_] = self.__dict__[sds_][~self.offview]

        # number of obs
        self.nobs = self.Longitude.size

def showdict(d,name):
    print('<<<<<< dictionary', name, 'begin >>>>>>')
    for key in d:
      t = type(d[key])
      print('>>', key, t)
      if t in (int, float, str):
        print('  scalar:', d[key])
      elif t in (tuple, list):
        print('  sequence:', d[key])
      elif t is np.ndarray:
        v = d[key]
        print('  numpy.ndarray: shape', v.shape)
        if v.ndim == 0:
          print('  scalar:', v)
        elif v.ndim == 1:
          print('  vect: first 5 vals:', v[:min(5,v.size)])
        elif v.ndim == 2:
          print('  array: first column:', v[0])
      print()
    print('<<<<<< dictionary', name, 'end >>>>>>')
    print()
    print()

#............................................................................

if __name__ == "__main__":

    sev03  = '/home/gwind/MSG+0000.3km.hdf'
    sevhgt = '/home/gwind/seviri_height.hdf'
    asm_Nx = '/nobackup/fp/opendap/asm_Nx.ctl'
    asm_Nv = '/nobackup/fp/opendap/asm_Nv.ctl'

    m = SEV03(sev03, sevhgt, Verb=1)
    m.getICAindx(asm_Nx)
    
def hold():

    m.linearSampleFile(asm_Nx,  onlyVars=('PS', 'T2M', 'SLP', 'U10M', 'V10M'))
    m.nearestSampleFile(asm_Nx, onlyVars=('QV2M',))

    m.linearSampleFile(asm_Nv, onlyVars=('T', 'DELP'))
    m.nearestSampleFile(asm_Nv, onlyVars=('O3','RH'))
    
                        
