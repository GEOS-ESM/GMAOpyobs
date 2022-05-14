"""
Reads Level 3 grid MXD43 BRDF files.

"""

import os
import sys
from numpy    import loadtxt, array, tile, where, concatenate
from numpy import ones
from datetime import date, datetime, timedelta
from glob     import glob
from pyhdf.SD import SD, HDF4Error
from .bits import BITS

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

class McD43(object):
    """
    This class implements the MODIS LAND BRDF daily Level 3 products, MCD43B1 (1000m horz res), 
    """

    def __init__ (self,Path,lon,lat,Verb=1):
       """
       Reads individual tile files for full day of Level 3 MCD43 
       present on a given *Path* and returns a objects with
       all 3 kernels coeff. On input, 

       Required parameters:
         Path -- can be a single file, a single directory, of a list
                 of files and directories. 
         
       Optional parameters:
         res  --- horizontal resolution of the file
                  by default : MCD43B files--> 1km res
                  MCD43A files --> 500 m                     
       """

       # List of HDF files for a given date
       #-----------------------------------
       self.verb = Verb
       self.SDS = SDS['LAND'] 
       self.Tfiles = glob(Path + '*.hdf')
      
    
       # From a list of lat and lon, return the tile numbers v(vertical), h(horiz),
       # and dx, dy inside the tile
       # -------------------------------------  
       self.nobs = len(lon)
       self._findTile(Path,lon,lat)
       
       # Create empty lists for SDS to be read from file
       # -----------------------------------------------
#       for name in self.SDS:
#           self.__dict__[name] = [] 

       # Read BRDF kernel in a MODIS tile
       # --------------------------------- 
       self.read_BRDF()
      
       # Result
       
            
       
#--- 
    def _findTile(self,path,lon,lat):

       """Given a list of lat, lon, return numbers to find the tile(v,h) 
          and position inside the tile (dx,dy)
       """
       from mpl_toolkits.basemap import Basemap # Basemap is deprecated, needs refatoring
       m = Basemap(projection='sinu',lon_0=0,rsphere=6371007.181,resolution='c')
       x,y=m(lon,lat)

       x_180,y_180 = m(-180,0)   # 4 points West, East, South, North
       x180,y180 = m(180,0)
       x_90,y_90 = m(0,-90)
       x90,y90 = m(0,90)
       
       dv = tile(int((y90-y_90)/18.),len(x))
      
       N = tile(18.,len(x))        # 18 vertical tiles in MODIS
       self.v = N-(y/dv)           # return number v tile for a list of lat, lon

       dh = tile(int((x180-x_180)/36.), len(x))
       self.h = x/dh
      
                       
       dirn, filen = os.path.split(self.Tfiles[0])   # verify dim of first file
       tokens = filen.split('.')    
       if tokens[0]=='MCD43B1' :
          xdim = 1200.
       else :
          print("- %s:not MCD43B file--> check resolution"%tokens[0])
       if self.verb:
          print('xdim', xdim)
       xdim_ = tile(xdim,len(x))

       int_h = [int(i) for i in self.h]
       int_v = [int(i) for i in self.v]
      
       dx = (self.h-int_h)*xdim_
       dy = (self.v-int_v)*xdim_

       self.dx = [int(i) for i in dx]
       self.dy = [int(i) for i in dy]
  
  
       self.h = int_h               # keep only real part
       self.v = int_v
       if self.verb:
          print('dx','dy', self.dx,self.dy)
      
       # create a list of tiles name associated with each (lat, lon) and (h,v)
       # -------------------
       self.Tiles = [] 
       for f in self.Tfiles :
                 dirn, filen = os.path.split(f)
                 tokens = filen.split('.')
                 for i in range(len(self.h)):
                     if (self.h[i] == int(tokens[2][1:3])) and (self.v[i] == int(tokens[2][4:6])) :
                         self.Tiles.append(f)

       # create a unique list of files
       # ----------------
       uniq = dict()
       for fn in array(self.Tiles) :
           uniq[fn] = 0 
#       self.unique_fn = uniq.keys()

       self.unique_fn = []
       # Get the index I for each file
       # ----------------------------
       for fn in list(uniq.keys()) :
           I = [array(self.Tiles) == fn]
           self.unique_fn.append((fn, I))

#---
    def get_Var(self, Vname = None):

        """Given filename dx, dy for a tile, return the 3 kernels for BRDF"""
        nch = 7
        BRDF = MISSING * ones((self.nobs,nch)) 
        for fn,I in self.unique_fn:
              self._read_BRDF(fn,I[0])
              if self.verb:
                print(I)
              BRDF
#---
    def read_BRDF(self):
       """Reads MCD43B1 one Tile file with Level 3 BRDF kernels for each MODIS band."""

       # Create empty lists for SDS to be read from file
       # -----------------------------------------------
       for name in self.SDS:
           self.__dict__[name] = [] 

       BRDF = MISSING * ones((len(self.SDS),self.nobs,3)) 

       for fn, I in self.unique_fn:
         index = I[0]
         if self.verb:
            print(index, type(index), len(index)) 
          # Don't fuss if the file cannot be opened
          # ---------------------------------------
         try:
            if self.verb:
                print("[] Working on "+fn)
            hfile = SD(fn)
         except HDF4Error:
            if self.verb > 2:
                print("- %s: not recognized as an HDF file"%filename)
            return 

          # Read select variables (reshape to allow concatenation later)
          # ------------------------------------------------------------
         for sds in self.SDS:  
            if self.verb:
              print('sds',self.SDS.index(sds))                 
            v = hfile.select(sds).get()           
            a = hfile.select(sds).attributes()
            if a['scale_factor']!=1.0 or a['add_offset']!=0.0:
                v = a['scale_factor'] * v + a['add_offset']
            if self.verb:
              print(array(self.dx)[index], BRDF.shape, BRDF[self.SDS.index(sds),index], v.shape) 

            BRDF[self.SDS.index(sds),index,:] = v[array(self.dx)[index], array(self.dy)[index], :]

       for sds in self.SDS:  
           self.__dict__[sds] = BRDF[self.SDS.index(sds),:,:]  
           if sds in list(ALIAS.keys()):
               self.__dict__[ALIAS[sds]] = self.__dict__[sds] 
         
        
#---

#............................................................................

if __name__ == "__main__":


#      path = '/nobackup/2/vbuchard/MODIS_LAND/MCD43B1/20070626/MCD43B1.A2007177.h%h2v%v2.*'
      path = '/nobackup/2/vbuchard/MODIS_LAND/MCD43B1/20070626/'
      
      lon=[-2.,-120.,15.2,17.2,170.1]
      lat=[88.,40.,-20.,-20.,-55.5]
      ex = McD43(path,lon,lat)

