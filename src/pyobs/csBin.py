"""

   Bin observations on cubed sphere grid. Based on externally generated 

"""

import numpy as np
import xarray as xr

from .constants import MAPL_UNDEF

class csBinError(Exception):
    """
    Defines general exception errors.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)



class CSBIN(object):
    
    def __init__(self,cs_filename):
        """
        Given *cs_filename* with cubed sphere grid coordinates.
        """
        
        self.ds = xr.open_dataset(cs_filename,chunks = {'nf': 1})
        
        # Center Cubed Sphere cordinates
        # ------------------------------
        self.lon_c = self.ds.coords['lons']
        self.lat_c = self.ds.coords['lats']
        self.Xdim  = self.ds.dims['Xdim']
        self.Ydim  = self.ds.dims['Ydim']
        # Coner coordinates are available as well.
        
    def set_Indices ( self, lons, lats):
        """
        Given a list of longitude and latitudes in degree (lons,lats), store
        corresponding cube sphere face F, and indices of the grid box
        where each observation falls into. From these indices, multiple
        variables can now be binned.
        """
        # some constants
        shift     = 0.174532925199433
        alpha     = 0.615479708670387
        sqr2      = 1.41421356237310
        dalpha    = 2.0*alpha/self.Xdim
        tolerance = 10**(-9)

        original_shp = lons.shape
        # Calculate (F,J,I)
        # shift the grid away from Japan Fuji Mt.
        lons_  = np.ravel(lons)/180*np.pi + shift
        lats_  = np.ravel(lats)/180*np.pi

        F = np.zeros(lons_.size, dtype=int)
        J = np.zeros(lons_.size, dtype=int)
        I = np.zeros(lons_.size, dtype=int)

        x = np.cos(lats_)*np.cos(lons_)
        y = np.cos(lats_)*np.sin(lons_)
        z = np.sin(lats_)

        max_abs = np.maximum(np.maximum(abs(x),abs(y)), abs(z))
        # project onto 2x2x2 cube
        x = x/max_abs
        y = y/max_abs
        z = z/max_abs

        f1Mask = (1-x) <= tolerance
        F[f1Mask] = 0
        I[f1Mask] = np.floor((np.arctan( y[f1Mask]/sqr2)+alpha)/dalpha)
        J[f1Mask] = np.floor((np.arctan( z[f1Mask]/sqr2)+alpha)/dalpha)

        f2Mask = (1-y) <= tolerance
        F[f2Mask] = 1
        I[f2Mask] = np.floor((np.arctan(-x[f2Mask]/sqr2)+alpha)/dalpha)
        J[f2Mask] = np.floor((np.arctan( z[f2Mask]/sqr2)+alpha)/dalpha)

        f3Mask = (1-z) <= tolerance
        F[f3Mask] = 2
        I[f3Mask] = np.floor((np.arctan(-x[f3Mask]/sqr2)+alpha)/dalpha)
        J[f3Mask] = np.floor((np.arctan(-y[f3Mask]/sqr2)+alpha)/dalpha)

        f4Mask = (x+1.0) <= tolerance
        F[f4Mask] = 3
        I[f4Mask] = np.floor((np.arctan(-z[f4Mask]/sqr2)+alpha)/dalpha)
        J[f4Mask] = np.floor((np.arctan(-y[f4Mask]/sqr2)+alpha)/dalpha)

        f5Mask = (y+1.0) <= tolerance
        F[f5Mask] = 4
        I[f5Mask] = np.floor((np.arctan(-z[f5Mask]/sqr2)+alpha)/dalpha)
        J[f5Mask] = np.floor((np.arctan( x[f5Mask]/sqr2)+alpha)/dalpha)

        f6Mask = (z+1.0) <= tolerance
        F[f6Mask] = 5
        I[f6Mask] = np.floor((np.arctan( y[f6Mask]/sqr2)+alpha)/dalpha)
        J[f6Mask] = np.floor((np.arctan( x[f6Mask]/sqr2)+alpha)/dalpha)

        self.F = F.reshape(*original_shp)
        self.I = I.reshape(*original_shp)
        self.J = J.reshape(*original_shp)

        return 
    
    def binObs ( self, obs ):
        """
        Given a list of longitude and latitudes in (lons,lats),
        bin list of observations *obs* on the cubed sphere,
        returning gridded observations *csObs*.
        """
        
        csObs = np.zeros(self.lon_c.shape, dtype=float)
        
        # One face at a time
        # ------------------
        for f in range(6):
            
            # Accumulator, counter
            # --------------------
            aObs = np.zeros(self.lon_c[0].shape, dtype=float)
            nObs = np.zeros(self.lon_c[0].shape, dtype=float)
            
            # Accumulate observations lying on this face of the cube
            # ------------------------------------------------------
            try:
                F = self.F==f
                I = self.I[F]
                J = self.J[F]
            except:
                raise csBinError('Indices (F,J,I) not yet set.')

            aObs[J,I] += obs[F]
            nObs[J,I] += 1
            
            # Normalize
            # ---------
            gObs = MAPL_UNDEF * np.ones(self.lon_c[0].shape, dtype=float)
            K = nObs>0
            gObs[K] = aObs[K] / nObs[K]
            
            csObs[f] = gObs
            
        # All done
        # --------
        return csObs    
            
#...........................................................................
if __name__ == '__main__':

   dyv2_dn  = '/css/g5nr/DYAMONDv2/'
   const_fn = dyv2_dn + '03KM/DYAMONDv2_c2880_L181/const_2d_asm_Mx/202002/DYAMONDv2_c2880_L181.const_2d_asm_Mx.20200201_0000z.nc4'
   csbin = CSBIN(const_fn)
   csbin.set_Indices(csbin.lon_c, csbin.lat_c)

   for f in range(6):#csbin.nf):
     for j in range(csbin.Ydim):
        for i in range(csbin.Xdim):
           if (csbin.F[f,j,i] != f):
              print("face is wrong", csbin.F[f,j, i], f, j, i )
           if (csbin.J[f,j,i] != j):
              print("J index is wrong", csbin.F[f,j, i], f, j, i )
           if (csbin.I[f,j,i] != i):
              print("I index is wrong", csbin.F[f,j, i], f, j, i )
