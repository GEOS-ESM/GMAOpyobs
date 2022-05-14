"""
Reads Geostationary pixel location/geometry file
"""

import os
import sys
from types import *

from netCDF4 import Dataset

from datetime import timedelta

import numpy as np

#---  

ALIAS = dict ( clon = 'lon', clat = 'lat' )

#...........................................................................
class GCS03Handle(object):
    """
    Generic container for GCS03.
    """
    def __init__(self,name):
        self.name = name

class GCS03(object):
    """
    This class implements the Geostationary GCS03 interface.
    """

    def __init__ (self,fgeom,varTime=False,Verb=0):
       """
       Reads Geostationary pixel location/geometry file.

       On input, 

       Required parameters:
         fgeom -- filename of pixel location/geometry file

       Optional parameters:
         varTime --- if True, read scanTime from file 
           (to allow intra-scan time interpolation)
           
         Verb      -- Verbose level:
                 0 - really quiet (default)
                 1 - Warns if invalid file is found
                 2 - Prints out non-zero number of aerosols in each file.
       """

       # Initially are lists of numpy arrays for each granule
       # ----------------------------------------------------
       self.verb	= Verb
       self.sat_lat     = None
       self.sat_lon     = None
       self.sample	= None # will hold sampled variables later
       self.ica		= None # will hold ICA indices later
       self.SDS_geom	= ('clon', 'clat')
       self.fgeom	= fgeom
       self.varTime	= varTime

       # Read file
       # ---------
       self._read(fgeom)

       # Make aliases for convenience
       # ----------------------------
       for sds in self.SDS_geom:
         sds_ = sds.replace(' ','_')
         if sds_ in ALIAS:
           self.__dict__[ALIAS[sds_]] = self.__dict__[sds_] 

#---
    def _read(self,fgeom):
        """Reads location/geometry file"""

        print("[] Working on "+fgeom)
        hf = Dataset(fgeom,'r')
        for sds in self.SDS_geom:
          sds_ = sds.replace(' ','_')
          self.__dict__[sds_] = hf.variables[sds][:]
        self.sat_lat = float(hf.sat_lat)
        self.sat_lon = float(hf.sat_lon)
        if self.varTime:
          # scanTime is a 1D array varying EW
          # units are seconds since start of hour
          scanTime = hf.variables['scanTime'][:]
          # replicate in NS direction
          self.scanTime = np.tile(scanTime, (self.clat.shape[0], 1))
        hf.close()

        # need this in case doesnt automatically pick up missing values
        self.clon = np.ma.masked_outside(self.clon, -180, 180)
        self.clat = np.ma.masked_outside(self.clat,  -90,  90)

        # store mask of off-view points
        self.offview = (np.ma.getmaskarray(self.clon) | np.ma.getmaskarray(self.clat))

        # store unflattened shape
        self.orgshape = self.offview.shape

        # store flattened, in-view SDS
        for sds in self.SDS_geom:
          sds_ = sds.replace(' ','_')
          self.__dict__[sds_] = self.__dict__[sds_][~self.offview]
        if self.varTime:
          self.scanTime = self.scanTime[~self.offview]

        # number of obs
        self.nobs = self.clon.size

#---
    def getICAindx(self,inFile):
        """
        Get grid related indices for ICA calculations. Calculates grid-box
        indices useful to implement Independent Column Approximation (ICA)
        type of algorithms.
        """
        from gfio import GFIOctl

        # Instiate grads and open file
        # ----------------------------
        f = GFIOctl(inFile)
        if self.ica is None:
            self.ica = GCS03Handle(inFile)

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
    def linearSampleFile(self, inFile, tyme, onlyVars=None):
        """
        Interpolates some or all variables of inFile.
        """
        from gfio import GFIOctl

        if self.varTime:
          # tyme is hour mark, and scanTime displaces from this in seconds
          tymes = np.array([tyme+timedelta(seconds=int(s)) for s in self.scanTime])
        else:
          # constant time
          tymes = np.repeat(tyme,self.nobs)
             
        # Open file
        # ---------
        f = GFIOctl(inFile)
        if self.sample == None:
            self.sample = GCS03Handle(inFile)
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

        if self.varTime:
          # tyme is hour mark, and scanTime displaces from this in seconds
          tymes = np.array([tyme+timedelta(seconds=int(s)) for s in self.scanTime])
        else:
          # constant time
          tymes = np.repeat(tyme,self.nobs)

        # Instiate grads and open file
        # ----------------------------
        f = GFIOctl(inFile)
        if self.sample == None:
            self.sample = GCS03Handle(inFile)
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

    gcs03  = '/discover/nobackup/pcastell/tempo_lg1_cld/tempo.lg1.cld.invariant.410.nc4'
    asm_Nx = '/nobackup/fp/opendap/asm_Nx.ctl'
    asm_Nv = '/nobackup/fp/opendap/asm_Nv.ctl'

    g = GCS03(gcs03, Verb=9)
#   g.getICAindx(asm_Nx)
