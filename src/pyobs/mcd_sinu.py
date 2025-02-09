"""
Reads Level 3 MODIS data on a sinosoidal tile grid as described at

    https://modis-land.gsfc.nasa.gov/MODLAND_grid.html

and provides functionality to sample the data at a list of (lon,lat).
It is assumed that all tile files for a give day/year are available on a single
directory.

"""

import os
import numpy       as np
import cartopy.crs as ccrs
import xarray      as xr
import pandas      as pd
    
from glob     import glob


#...........................................................................

# The 2 relevant transforms
# -------------------------
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

__DEBUG__ = False

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

def _sinu2ll(x,y):
    """
    Given (x,y) MODIS sinusoidal coordinates, returns (lon,lat) coordinates.
    On input, (x,y) can be scalars or arrays.
    """

    # Input are arrays
    # ----------------
    if isinstance(lon,np.ndarray):
        z = latlon.transform_points(modis, x, y)
        return (z[:,0], z[:,1]) # x, y coordinartes

    # Input are scalars
    # -----------------
    else:
        return latlon.transform_point(x, y, modis)

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

class NodataError(Exception):
    """
    Defines Nodata exceptions.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class MCD_SINU(object):
    """
    This class implements the MODIS LAND Level 3 products on the tiled sinosoidal grid. 
    """

    def __init__ (self,Path,lon,lat,addLatLon=False, Alias=None):
       """
       Reads individual tile files for full day of Level 3 
       tiled files on a simosoidal grid. Files assumed to be 
       present on a given *Path* and returns a objects with
       all 3 kernels coeff.  

       On input,
         Path : directory where to find tile files for a given day         
         lon, lat : coordinates to sample variables on
         addLatLon : whether to add 2D lat/lon coordinates to tiles;
                     this is only necessary for testing interpolation,
                     default is False.
         Alias : an optional dictionary for renaming variables on file.
       
       """

       # From a list of lat and lon, return the tile numbers v(vertical), h(horiz)
       # global coordinates (x,y) inside each tile
       # -------------------------------------------------------------------------  
       self.nobs = len(lon)
       self._findTile(lon,lat)
       
       # Lazy load each needed MCD_SINU file, retrieving bounding boxes
       # and assigning missing coordinate variables
       # -----------------------------------------------------------
       self._lazyLoadTiles(Path,addLatLon=addLatLon,Alias=Alias)


#---
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

       self.obsTileIndex = dict()          # index of those observations on a given tile
       T1D = v*nH + h
       for t1d in np.unique(T1D):          # loop over unique set of tiles
           v = t1d // nH
           h = t1d % nH
           tn = 'h%02dv%02d'%(h,v)             # tile name, e.g., h00v18
           self.obsTileIndex[tn] = (T1D==t1d)  # save obs indices corresponding to this tile
           #print(tn,t1d,any(self.obsTileIndex[tn]))
        
#---
    def _lazyLoadTiles(self,Path,addLatLon=False,Alias=None):
        """
        Lazy load each tile file, retrieve bounding boxes
        and add coordinate variables to xarray object.
        """
        neededTiles = list(self.obsTileIndex.keys()) 
        #print(neededTiles)
        TileFileNames = glob(Path + '/*.hdf')
        self.Tiles = dict()
        for fn in TileFileNames:
            tn = os.path.basename(fn).split('.')[2]
            if tn not in neededTiles: continue # only load what is needed
            
            ds = xr.open_dataset(fn,engine='netcdf4')
            
            # Give dimensions sensible names
            # ------------------------------
            ren_ = dict()
            for dn in ds.dims:
                if 'XDim' in dn: ren_[dn] = 'x'
                if 'YDim' in dn: ren_[dn] = 'y'
                if 'Num_Parameters:MOD_Grid_BRDF' in dn: ren_[dn] = 'k'

            if Alias is not None:
                ren_.update(Alias)

            ds = ds.rename(ren_)

            self.Tiles[tn] = ds 

        if len(self.Tiles) == 0:
            raise NodataError("No valid MCD_SINU tiles for this day.")
            
        # Add coordinate variables
        # ------------------------
        self.nx, self.ny = ds.dims['x'], ds.dims['y']
        for tn in self.Tiles:
            xs, ys = _tn2bbox(tn)
            x = np.linspace(xs[0],xs[1],self.nx,endpoint=True)
            y = np.linspace(ys[1],ys[0],self.ny,endpoint=True)  # Vertical gris is North to South

            ds = self.Tiles[tn]
            
            coords_ = dict()
            if addLatLon:

                X, Y = np.meshgrid(x,y)
                Lon, Lat = _sinu2ll(X.flatten(),Y.flatten())
                Lon = xr.DataArray(Lon.reshape(X.shape),dims=('x','y'))
                Lat = xr.DataArray(Lat.reshape(Y.shape),dims=('x','y'))

                coords_ = {  'x': x, 'y':y,          # sinusoidal coords
                             'lon':Lon, 'lat':Lat}   # and 2D lat/lon coordinates
                                                   
            else:
                coords_ = {  'x': x, 'y':y }    # sinusoidal coords

            self.Tiles[tn] = ds.assign_coords(**coords_)  # assign coordinates

        # Convenience list of variables, coordinates and dimensions
        # ---------------------------------------------------------
        self.coords, self.variables, self.dims = list(ds.coords), list(ds.data_vars), list(ds.dims)
        self.shapes = dict()
        for vn in self.variables:
            self.shapes[vn] = ds[vn].shape
#---
    def interp(self, vname, method='nearest'):
        """
        Given variable name 'vname', sample said variable on
        observations locations.

        Parameters
        ----------
        vname  : str, variable name
        method : {"linear", "nearest", "zero", "slinear", "quadratic", "cubic", "polynomial"}, default: "nearest"
        
        """

        # Container for output
        # --------------------
        var = xr.DataArray(np.ones((self.nobs,3))+np.nan,
                                     coords = { 'lon': ('nobs',self.lon),
                                                'lat': ('nobs',self.lat)},
                                     dims=('nobs','k') )
           
        for tn in self.obsTileIndex:

           I  = self.obsTileIndex[tn]  # obs indices on this tile

           x, y = self.x[I], self.y[I] # obs coords on this tile

           try:
               ds = self.Tiles[tn]         # xarray corresponding to tile

               if __DEBUG__:
                   X, Y = ds.coords['x'], ds.coords['y']
                   if x.min()<X.min() or x.max()>X.max():
                       raise ValueError('x out of range for '+tn)
                   if y.min()<Y.min() or y.max()>Y.max():
                       raise ValueError('y out of range for '+tn)

               if method == 'nearest':
                   vinterp = ds[vname].sel(x=x, y=y,method=method)
               else:               
                   vinterp = ds[vname].interp(x=x, y=y,method=method)

               if len(vinterp.shape) == 1:
                   var[I,0] = vinterp          
               elif len(vinterp.shape) == 2:   
                   var[I,:] = vinterp
               else:
                   raise Warning('Internal error, this should never happen!')

           except:
               pass  # values will remain underf
            
        if len(self.shapes[vname]) == 2:
           return var[:,0]               
        else:
           return var

#---
    def interp_many(self,Variables=None,Index=None):
        """
        Sample all variables on file, returning an xarray Data DataFrame with all
        interpolated variables. On input,

        Variables : list of variables to interpolate (default: all variables)
        Index : Pandas DataFrame index; if not specified, an Xarray dataset is returned.
      
        """
        if Variables is None:
            Variables = self.variables
        variables = dict()
        if Index is not None:       
            variables['lon'], variables['lat'] = self.lon, self.lat
        for vname in Variables:
            variables[vname] = self.interp(vname)
     
        if Index is None:       
            return xr.Dataset(variables)
        else:
            return pd.DataFrame(variables,index=Index)
#---
    def _getBoundsFromMetatada(self,ds):
        """
        Given a MCD_SINU granule xarray dataset, returns bounds
        found in the metadata. Useful to verify correctness of
        _tn2bbox(), not really used otherwise.
        """

        meta = ds.attrs['ArchiveMetadata.0']
        lons = np.array((float(meta.split('WESTBOUNDINGCOORDINATE')[1].split('=')[2].split('\n')[0]),
                         float(meta.split('EASTBOUNDINGCOORDINATE')[1].split('=')[2].split('\n')[0]) ))
        lats = np.array((float(meta.split('SOUTHBOUNDINGCOORDINATE')[1].split('=')[2].split('\n')[0]),
                         float(meta.split('NORTHBOUNDINGCOORDINATE')[1].split('=')[2].split('\n')[0]) ))

        meta = ds.attrs['StructMetadata.0']
        xs, ys = np.zeros(2), np.zeros(2)
        xs[0],ys[1] = np.array(meta.split('UpperLeftPointMtrs')[1].split('\n')[0].replace('=(','').replace(')','').split(',')).astype('float')
        xs[1],ys[0] = np.array(meta.split('LowerRightMtrs')[1].split('\n')[0].replace('=(','').replace(')','').split(',')).astype('float')

        return (lons, lats, xs, ys)

