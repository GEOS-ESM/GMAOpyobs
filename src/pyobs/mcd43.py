"""
Reads Level 3 grid MXD43 BRDF files.

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

class McD43(object):
    """
    This class implements the MODIS LAND BRDF daily Level 3 products, MCD43B1 (1000m horz res), 
    """

    def __init__ (self,Path,lon,lat,addLatLon=False):
       """
       Reads individual tile files for full day of Level 3 MCD43 
       present on a given *Path* and returns a objects with
       all 3 kernels coeff.  

       On input,
         Path : directory where to find tile files for a given day         
         lon, lat : coordinates to sample variables on
         addLatLon : whether to add 2D lat/lon coordinates to tiles;
                     this is only necessary for testing interpolation,
                     default is False.
       
       """

       # From a list of lat and lon, return the tile numbers v(vertical), h(horiz)
       # global coordinates (x,y) inside each tile
       # -------------------------------------------------------------------------  
       self.nobs = len(lon)
       self._findTile(lon,lat)
       
       # Lazy load each needed MCD43 file, retrieving bounding boxes
       # and assigning missing coordinate variables
       # -----------------------------------------------------------
       self._lazyLoadTiles(Path,addLatLon=addLatLon)


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
           tn = 'h%02dv%02d'%(h,v)            # tile name, e.g., h00v18
           self.obsTileIndex[tn] = T1D==t1d   # save obs indices corresponding to this tile
               
#---
    def _lazyLoadTiles(self,Path,addLatLon=False):
        """
        Lazy load each tile file, retrieve bounding boxes
        and add coordinate variables to xarray object.
        """
        neededTiles = list(self.obsTileIndex.keys()) 
        TileFileNames = glob(Path + '*.hdf')
        self.Tiles = dict()
        for fn in TileFileNames:
            tn = os.path.basename(fn).split('.')[2]
            if tn not in neededTiles: continue # only load what is needed
            
            ds = xr.open_dataset(fn,engine='netcdf4')
            
            # Give dimensions sensible names
            # ------------------------------
            if 'Num_Parameters:MOD_Grid_BRDF' in ds.dims:
                ds = ds.rename({'XDim:MOD_Grid_BRDF': 'x',
                                'YDim:MOD_Grid_BRDF': 'y',
                                'Num_Parameters:MOD_Grid_BRDF': 'k'})
            else:
                ds = ds.rename({'XDim:MOD_Grid_BRDF': 'x',
                                'YDim:MOD_Grid_BRDF': 'y'})
 
            self.Tiles[tn] = ds

        # Add coordinate variables
        # ------------------------
        self.nx, self.ny = ds.dims['x'], ds.dims['y']
        for tn in self.Tiles:
            xs, ys = _tn2bbox(tn)
            x = np.linspace(xs[0],xs[1],self.nx,endpoint=True)
            y = np.linspace(ys[1],ys[0],self.ny,endpoint=True)  # Vertical gris is North to South

            ds = self.Tiles[tn]
            
            if addLatLon:
                X, Y = np.meshgrid(x,y)
                Lon, Lat = _sinu2ll(X.flatten(),Y.flatten())
                Lon = xr.DataArray(Lon.reshape(X.shape),dims=('x','y'))
                Lat = xr.DataArray(Lat.reshape(Y.shape),dims=('x','y'))

                self.Tiles[tn] = ds.assign_coords({  'x': x,    'y':y,    # sinusoidal coords
                                                   'lon':Lon, 'lat':Lat}) # and 2D lat/lon coordinates
                                                   
            else:
                self.Tiles[tn] = ds.assign_coords({'x': x, 'y':y})        # sinusoidal coords

        # Convenience list of variables, coordinates and dimensions
        # ---------------------------------------------------------
        self.coords, self.variables, self.dims = list(ds.coords), list(ds.variables), list(ds.dims)
        for v in self.variables:
            if v in self.coords:
                self.variables.remove(v) # exclude coordinate variables
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
        var = xr.DataArray(np.ones((self.nobs,3)),
                                     coords = { 'lon': ('nobs',self.lon),
                                                'lat': ('nobs',self.lat)},
                                     dims=('nobs','k') )
           
        for tn in self.obsTileIndex:

           I  = self.obsTileIndex[tn]  # obs indices on this tile
           ds = self.Tiles[tn]         # xarray corresponding to tile
           
           x, y = self.x[I], self.y[I] # obs coords on this tile

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

        if len(vinterp.shape) == 1:
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
        Given a MCD43 granule xarray dataset, returns bounds
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

#............................................................................

if __name__ == "__main__":

      fluxnet_fn = '/Users/adasilva/data/brdf/fluxnet_stations.csv'
      mcd43a1_dn = '/Users/adasilva/data/brdf/MCD43A4/2023/204/'

      stations = pd.read_csv('/Users/adasilva/data/brdf/fluxnet_stations.csv',
                             index_col=0)

      print(stations)

      lon = stations['lons'].values
      lat = stations['lats'].values

      # BRDF object properly initialized
      # --------------------------------
      brdf = McD43(mcd43a1_dn,lon,lat)

      # Sample 1 variable at a time at obs locations
      # --------------------------------------------
      for vname in brdf.variables:
            print('[] Interpolating',vname)
            v = brdf.interp(vname)

      # Alternativaly, interpolate several variables
      # Note: omit variables to interpolate all variables
      # -------------------------------------------------
      Variables = ['Nadir_Reflectance_Band1', 'Nadir_Reflectance_Band2']
      V = brdf.interp_many(Variables,Index=stations.index)
      print('\n Dataframe:')
      print(V)


