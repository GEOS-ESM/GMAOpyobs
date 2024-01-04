"""
Reads Level 3 grid MXD43 BRDF files.

"""

import os
import sys
import numpy       as np
import cartopy.crs as ccrs
import xarray      as xr

from datetime import date, datetime, timedelta
from glob     import glob
#from pyhdf.SD import SD, HDF4Error

MISSING = -99999


SDS = dict (
      LAND = ('BRDF_Albedo_Parameters_Band1','BRDF_Albedo_Parameters_Band2',
              'BRDF_Albedo_Parameters_Band3','BRDF_Albedo_Parameters_Band4',
              'BRDF_Albedo_Parameters_Band5','BRDF_Albedo_Parameters_Band6',
              'BRDF_Albedo_Parameters_Band7'),

      QUAL = ('BRDF_Albedo_Quality',
              'Snow_BRDF_Albedo',
              'BRDF_Albedo_Ancillary', )
           )

ALIAS = dict ( BRDF_Albedo_Parameters_Band1 = 'BRDF_b1_645',
               BRDF_Albedo_Parameters_Band2 = 'BRDF_b2_856',
               BRDF_Albedo_Parameters_Band3 = 'BRDF_b3_465',
               BRDF_Albedo_Parameters_Band4 = 'BRDF_b4_553',
               BRDF_Albedo_Parameters_Band5 = 'BRDF_b5_1241',
               BRDF_Albedo_Parameters_Band6 = 'BRDF_b6_1629',
               BRDF_Albedo_Parameters_Band7 = 'BRDF_b7_2114',
             )

#...........................................................................

modis  = ccrs.Sinusoidal.MODIS # Built in
latlon = ccrs.PlateCarree()

# 4 points West, East, South, North
# ---------------------------------
nH = 36
nV = 18

x_180, dummy = modis.transform_point(-180,0,latlon)   
x180,  dummy = modis.transform_point(180,0,latlon)   
dummy,  y_90 = modis.transform_point(0,-90,latlon)   
dummy,  y90  = modis.transform_point(0,90,latlon)   
       
dH = (x180-x_180)/nH     # delta-h, size of each horizontal tile
dV = (y90-y_90)/nV       # delta-v, size of each vertical tile

    
#...........................................................................
  
def _ll2sinu(lon,lat):
    """
    Given (lon,lat) coordinates, returns (x,y) coordinates of MODIS sinusoidal
    transform. On input, (lon,lat) can be scalars or arrays.
    """

    # Input are arrays
    # ----------------
    if isinstance(lon,np.ndarray):
        z = modis.transform_points(latlon, lon, lat)
        return (z[:,0], z[:,1]) # x, y coordinartes

    # Input are scalars
    # -----------------
    else:
        return modis.transform_point(lon,lat, latlon)

def _tn2bbox(TileName):
    """
    Given a tile name, return bounding boxes.
    TileName example: h32v08
    """
    h = int(TileName[1:3])
    v = nV - 1 - int(TileName[4:6]) # vertical tiles are upside down!

    xs, ys = np.zeros(2), np.zeros(2)
     
    xs[0] = x_180 + h*dH
    xs[1] = xs[0] + dH
    
    ys[0] = y_90  + v*dV
    ys[1] = ys[0] + dV

    return (xs,ys)
    

#...........................................................................

class McD43(object):
    """
    This class implements the MODIS LAND BRDF daily Level 3 products, MCD43B1 (1000m horz res), 
    """

    def __init__ (self,Path,lon,lat,Verb=1):
       """
       Reads individual tile files for full day of Level 3 MCD43 
       present on a given *Path* and returns a objects with
       all 3 kernels coeff.  

       On input,
         Path -- can be a single file, a single directory, of a list
                 of files and directories. 
         
         lon, lat -- coordinates to sample BRDF kernels
       
       """

       # List of HDF files for a given directory,
       # typically corresponding to a date.
       #-----------------------------------------
       self.verb = Verb
       self.SDS = SDS['LAND'] 

       # Lazyload each MCD43 file, retrieving bounding boxes
       # and assigning missing coordinate variables
       # ---------------------------------------------------
       self._lazyload(Path)

       # From a list of lat and lon, return the tile numbers v(vertical), h(horiz)
       # global coordinates (x,y) inside each tile
       # -------------------------------------------------------------------------  
       self.nobs = len(lon)
       self._findTile(lon,lat)
       
       # Read BRDF kernels at each observation location
       # ----------------------------------------------
       self.interpBRDF()
      
          
#---
    def _lazyload(self,Path):
        """
        Lazy load each tile file, retrieve bounding boxes
        and add coordinate variables to xarray object.
        """
        TileFileNames = glob(Path + '*.hdf')
        self.Tiles = dict()
        for fn in TileFileNames:
            tn = os.path.basename(fn).split('.')[2]
            ds = xr.open_dataset(fn,engine='netcdf4')
            
            ## xs, ys = np.zeros(2), np.zeros(2)
            ## meta = ds.attrs['ArchiveMetadata.0']
            ## lons = np.array((float(meta.split('WESTBOUNDINGCOORDINATE')[1].split('=')[2].split('\n')[0]),
            ##                  float(meta.split('EASTBOUNDINGCOORDINATE')[1].split('=')[2].split('\n')[0]) ))
            ## lats = np.array(( float(meta.split('SOUTHBOUNDINGCOORDINATE')[1].split('=')[2].split('\n')[0]),
            ##                   float(meta.split('NORTHBOUNDINGCOORDINATE')[1].split('=')[2].split('\n')[0]) ))
            
            ## meta = ds.attrs['StructMetadata.0']
            ## xs[0],ys[1] = np.array(meta.split('UpperLeftPointMtrs')[1].split('\n')[0].replace('=(','').replace(')','').split(',')).astype('float')
            ## xs[1],ys[0] = np.array(meta.split('LowerRightMtrs')[1].split('\n')[0].replace('=(','').replace(')','').split(',')).astype('float')

            ## self.Tiles[tn] = dict ( ds=ds, lons=lons, lats=lats, xs=xs, ys=ys )

            # Give dimensions sensible names
            # ------------------------------
            ds = ds.rename({'XDim:MOD_Grid_BRDF': 'x',
                            'YDim:MOD_Grid_BRDF': 'y',
                            'Num_Parameters:MOD_Grid_BRDF': 'k'})
            self.Tiles[tn] = ds

        # Add coordinate variables
        # ------------------------
        self.nx, self.ny = ds.dims['x'], ds.dims['y']
        for tn in self.Tiles:
            xs, ys = _tn2bbox(tn)
            x = np.linspace(xs[0],xs[1],self.nx,endpoint=True)
            y = np.linspace(ys[1],ys[0],self.ny,endpoint=True)  # Vertical gris is North to South
            self.Tiles[tn] = ds.assign_coords({'x': x, 'y':y})  # dataset now has sinusoidal coordinates with sensible names
            
    def _findTile(self,lon,lat):

       """
          Given a list of lat, lon, find tile coordinates (v,h) and the
       unique set of tiles that contain those observations.
       """

       self.lon, self.lat = lon, lat
       x, y = _ll2sinu(lon,lat)  # corresponding sinusoidal coordinates

       # Make these DataArrays with fake dimension
       # -----------------------------------------
       self.x = xr.DataArray(x, dims='nobs')
       self.y = xr.DataArray(y, dims='nobs')

       v  = nV - 1 - (self.y-y_90)//dV     # tile vertical coordinate v in [0,18]
       h  = (self.x-x_180)//dH             # tile horizontal coordinate in [0,36]

       self.obsTileIndex = dict()        # index of those observations on a given tile
       T1D = v*nH + h
       for t1d in np.unique(T1D):        # loop over unique set of tiles
           v = t1d // nH
           h = t1d % nH
           tn = 'h%02dv%02d'%(h,v)            # tile name, e.g., h00v18
           self.obsTileIndex[tn] = T1D==t1d   # save obs indices corresponding to this tile
       
#---
    def interpBRDF(self):
        """
        Sample BRDF kernels from MCD43A1/B1 at observation locations.
        """

        for sds in self.SDS:
           self.__dict__[sds] = xr.DataArray(np.ones((self.nobs,3)),
                                             coords = { 'lon': ('nobs',self.lon),
                                                        'lat': ('nobs',self.lat)},
                                             dims=('nobs','k'))
           #print(sds,self.__dict__[sds])
    
        for tn in self.obsTileIndex:

           I  = self.obsTileIndex[tn]  # obs indices on this tile
           ds = self.Tiles[tn]         # xarray corresponding to tile
           
           x, y = self.x[I], self.y[I] # obs coords on this tile

           print(tn,'lon = ',self.lon[I], 'lat = ',self.lat[I])
           
           for sds in self.SDS:
               self.__dict__[sds][I,:] = ds[sds].sel(x=x, y=y,method='nearest').values

        # Convenient aliases
        for sds in self.SDS:
           #print(sds,self.__dict__[sds])
           if sds in list(ALIAS.keys()):
               self.__dict__[ALIAS[sds]] = self.__dict__[sds] 
         

#............................................................................

def fluxnet(filen):

    import pandas as pd
    
    flx = xr.open_dataset(filen)

    Lon  = flx['Longitude']
    Lat  = flx['Latitude']
    Name = flx['Name']

    Stations = np.unique(Name)

    lons, lats = [], []
    for stn in Stations:

        I = Name==stn

        lons += [Lon[I][0],]
        lats += [Lat[I][0],]

    variables = dict( lons=np.array(lons), lats=np.array(lats) )
    return pd.DataFrame( variables,index=Stations)     
    

if __name__ == "__main__":

      fluxnet_fn = '/Users/adasilva/data/brdf/onefluxnet_daily_mcd43_c61.nc'
      mcd43a1_dn = '/Users/adasilva/data/brdf/2023/204/'

      stations = fluxnet(fluxnet_fn)

      print(stations)

      lon = stations['lons'].values
      lat = stations['lats'].values
      
      brdf = McD43(mcd43a1_dn,lon,lat)


def xxx():
      
      for tn in brdf.Tiles:
            ds = brdf.Tiles[tn]
            xs, ys = _tn2bbox(tn)
            x, y = ds['x'], ds['y']
            xs = np.array((x.min(), x.max()))
            ys = np.array((y.min(), y.max()))
            lons, lats = np.zeros(2), np.zeros(2)
            lons[0],lats[0] = latlon.transform_point(xs[0],ys[0],modis)
            lons[1],lats[1] = latlon.transform_point(xs[1],ys[1],modis)
            print(tn)
            print('lons = ',lons, 'lats = ',lats)
            #print('  xs = %20.3f %20.3f %20.3f'%(tile['xs'][0],tile['xs'][1],tile['xs'][1]-tile['xs'][0]))
            print('  XS = %20.3f %20.3f %20.3f'%(xs[0],xs[1],xs[1]-xs[0]))
            print('  YS = %20.3f %20.3f %20.3f'%(ys[0],ys[1],ys[1]-ys[0]))

      print()
      print('x limits = ',x_180,x180,dH)
      print('y limits = ',y_90,y90,dV)
      
