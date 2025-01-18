"""

   Bin observations on LatLon grid. 

"""

import numpy as np
import bisect
import xarray as xr

from .constants import MAPL_UNDEF

class llBinError(Exception):
    """
    Defines general exception errors.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class LLBIN(object):
    
    def __init__(self, **kwargs):
        """
        Supporting keywords: GridName,lon_bnds, lat_bnds.
        GridName: It should have the standard format : POLEiiiXjjj-DATELINE
           POLE      : PE or PC
           iii       : number of lon
           x         : the separator
           jjj       : number of lat
           DATELINE  : DE or DC

           For example, PE360x180-DE represents a lat-lon grid with 360 lons, 180 lats, pole on edge and Dateline on edge

        lon_bnds: longitude bounds, an ordered array from -180. to 180. 
        lat_bnds: latitude bounds, an ordered array from -90.  to 90.
        """

        kw_ = {k.upper():v for k,v in kwargs.items()}        

        if 'GRIDNAME' in kw_.keys():
           grid = kw_['GRIDNAME'].upper()
           pole = grid[:2]
           x    = grid.find('X')
           dateline = grid[-2:]
           nlon = int(grid[2:x])
           nlat = int(grid[x+1:-3])
           self.pole= pole
           self.dateline=dateline
           self.Xdim = nlon   
           self.Ydim = nlat
        elif ( 'LAT_BNDS' in kw_.keys() and 'LON_BNDS' in kw_.keys()):
           self.lon_bnds = lon_bnds
           self.lat_bnds = lat_bnds
           self.pole='unknown'
           self.dateline='unknown'
           self.Xdim = len(self.lon_bnds)-1   
           self.Ydim = len(self.lat_bnds)-1
        else:
           raise llBinError('Init keywork is not correct')

    def setIndices (self, lons, lats):
        """
        Given a list of longitude and latitudes in degree (lons,lats), 
        find the cell indices.
        """
        I = np.zeros(len(lons))
        J = np.zeros(len(lats))
        # get I indices
        if (self.dateline =='DE'):
           dlon = 360./self.Xdim
           I    = np.floor((lons+180.)/dlon).astype(int)
        elif (self.dateline =='DC'):
           dlon = 360./self.Xdim
           lons = lons + dlon/2.
           lons = np.where(lons >= 180.0, lons - 360., lons)
           I    = np.floor((lons+180.0)/dlon).astype(int)
        elif (self.dateline =='unknown'):
           I    = np.array([ bisect.bisect(self.lon_bnds, lon)  for lon in lons]) - 1
        # get I indices
        if (self.pole == 'PE'):
           dlat = 180./self.Ydim
           J    = np.floor((lats + 90.)/dlat).astype(int)
        elif (self.pole == 'PC'):
           dlat = 180./(self.Ydim-1)  
           J    = np.floor((lats + 90.+ dlat/2.)/dlat).astype(int)
        elif (self.pole == 'unknown'):
           J    = np.array([ bisect.bisect(self.lat_bnds, lat)  for lat in lats]) - 1
        # m
        I = np.where( I == self.Xdim, 0,   I)
        J = np.where( J == self.Ydim, J-1, J)

        self.IJ = I + J*self.Xdim
 
        return 
    
    def binObs ( self, obs, average = True ):
        """
        Given a list of longitude and latitudes in (lons,lats),
        bin list of observations *obs* on the lat-lon grid,
        returning gridded observations *llObs*.
        """
        
        llObs     = np.zeros((self.Xdim, self.Ydim), dtype=float)
        grid_size = self.Xdim * self.Ydim
        aObs      = np.zeros(grid_size, dtype=float)
        nObs      = np.zeros(grid_size, dtype=int)
       
        if (len(obs.shape) > 1):
           obs = np.ravel(obs) 

           # Accumulator, counter
           # --------------------
           
           # Accumulate observations lying on this grid
           # ------------------------------------------------------
           try:
               IJ = self.IJ
           except:
               raise llBinError('Indices (F,IJ) not yet set.')

           np.add.at(nObs, IJ, 1)
           np.add.at(aObs, IJ, obs)
           
           # Normalize
           # ---------
           gObs = MAPL_UNDEF * np.ones(grid_size, dtype=float)
           K = nObs>0

           if (average) :
              gObs[K] = aObs[K] / nObs[K]
           else :
              gObs[K] = aObs[K]

           llObs = np.reshape(gObs, (self.Xdim, self.Xdim))
        # All done
        # --------
        return llObs    

#...........................................................................
if __name__ == '__main__':

   llbin = LLBIN(GridName='PE360x180-DE')
   nlon = 360
   nlat = 180
   dlat =  180./nlat
   dlon =  360./nlon

   lats  = xr.DataArray([  -90.])
   lons  = xr.DataArray([ -180.])
   llbin.setIndices(lons, lats)
   print(llbin.IJ[0] == 0)

   lats  = xr.DataArray([  -90. + 1.5*dlat ])
   lons  = xr.DataArray([ -180.])
   llbin.setIndices(lons, lats)
   print(llbin.IJ[0] == 360)

   lats  = xr.DataArray([  -90. +1.5*dlat ])
   lons  = xr.DataArray([ -180. +1.5*dlon ])
   llbin.setIndices(lons, lats)
   print(llbin.IJ[0] == 361)

   lats  = xr.DataArray([  -90. +10.5*dlat ])
   lons  = xr.DataArray([ -180. +10.5*dlon ])
   llbin.setIndices(lons, lats)
   print(llbin.IJ[0] == 3610)

   lon_bnds = xr.DataArray([ -180.0+i*dlon for i in range(nlon+1)])
   lat_bnds = xr.DataArray([ -90.0 +i*dlat for i in range(nlat+1)])

   llbin = LLBIN(lon_bnds=lon_bnds, lat_bnds=lat_bnds)

   lats  = xr.DataArray([  -90. + 0.1*dlat])
   lons  = xr.DataArray([ -180. + 0.1*dlon])
   llbin.setIndices(lons, lats)
   print(llbin.IJ[0] == 0 )

   lats  = xr.DataArray([  -90. +10.5*dlat ])
   lons  = xr.DataArray([ -180. +10.5*dlon ])
   llbin.setIndices(lons, lats)
   print(llbin.IJ[0] == 3610)
