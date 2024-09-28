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

    def set_Indices ( self, lons, lats):
        """
        Given a list of longitude and latitudes in degree (lons,lats), store
        corresponding cube sphere face F, and indices of the grid box
        where each observation falls into. From these indices, multiple
        variables can now be binned.
        """
        # some constants
        shift     = 0.174532925199433
        tolerance = 10**(-6) # on 0-1 scale
        alpha     = 0.615479708670387
        dalpha    = 2.0*alpha/self.Xdim
        sqr2      = np.sqrt(2.)

        # Calculate (F,J,I)

        # shift the grid away from Japan Fuji Mt.
        lons_  = lons/180*np.pi + shift
        lats_  = lats/180*np.pi
        # some functions only work for 1-d array. ravel it first
        if len(lons.shape) > 1:
          lons_  = np.ravel(lons_)
          lats_  = np.ravel(lats_)
           
        self.F = np.zeros(lons_.size, dtype=int)
        self.J = np.zeros(lons_.size, dtype=int)
        self.I = np.zeros(lons_.size, dtype=int)

        x = np.cos(lats_)*np.cos(lons_)
        y = np.cos(lats_)*np.sin(lons_)
        z = np.sin(lats_)

        max_abs = np.maximum(np.maximum(abs(x),abs(y)), abs(z))
        # project onto 2x2x2 cube
        x = x/max_abs
        y = y/max_abs
        z = z/max_abs

        f1Mask = (1-x) <= tolerance
        self.F[f1Mask] = 0
        self.I[f1Mask] = np.floor((np.arctan( y[f1Mask]/sqr2)+alpha)/dalpha)
        self.J[f1Mask] = np.floor((np.arctan( z[f1Mask]/sqr2)+alpha)/dalpha)

        f2Mask = (1-y) <= tolerance
        self.F[f2Mask] = 1
        self.I[f2Mask] = np.floor((np.arctan(-x[f2Mask]/sqr2)+alpha)/dalpha)
        self.J[f2Mask] = np.floor((np.arctan( z[f2Mask]/sqr2)+alpha)/dalpha)

        f3Mask = (1-z) <= tolerance
        self.F[f3Mask] = 2
        self.I[f3Mask] = np.floor((np.arctan(-x[f3Mask]/sqr2)+alpha)/dalpha)
        self.J[f3Mask] = np.floor((np.arctan(-y[f3Mask]/sqr2)+alpha)/dalpha)

        f4Mask = (x+1.0) <= tolerance
        self.F[f4Mask] = 3
        self.I[f4Mask] = np.floor((np.arctan(-z[f4Mask]/sqr2)+alpha)/dalpha)
        self.J[f4Mask] = np.floor((np.arctan(-y[f4Mask]/sqr2)+alpha)/dalpha)

        f5Mask = (y+1.0) <= tolerance
        self.F[f5Mask] = 4
        self.I[f5Mask] = np.floor((np.arctan(-z[f5Mask]/sqr2)+alpha)/dalpha)
        self.J[f5Mask] = np.floor((np.arctan( x[f5Mask]/sqr2)+alpha)/dalpha)

        f6Mask = (z+1.0) <= tolerance
        self.F[f6Mask] = 5
        self.I[f6Mask] = np.floor((np.arctan( y[f6Mask]/sqr2)+alpha)/dalpha)
        self.J[f6Mask] = np.floor((np.arctan( x[f6Mask]/sqr2)+alpha)/dalpha)

        self.IJ = self.I + self.J*self.Xdim

        return 
    
    def binObs ( self, obs ):
        """
        Given a list of longitude and latitudes in (lons,lats),
        bin list of observations *obs* on the cubed sphere,
        returning gridded observations *csObs*.
        """
        
        csObs = np.zeros(self.lon_c.shape, dtype=float)
       
        if (len(obs.shape) > 1):
           obs = np.ravel(obs) 
        # One face at a time
        # ------------------
        for f in range(6):
            
            # Accumulator, counter
            # --------------------
            aObs = np.zeros(self.lon_c[0].size, dtype=float)
            nObs = np.zeros(self.lon_c[0].size, dtype=float)
            
            # Accumulate observations lying on this face of the cube
            # ------------------------------------------------------
            try:
                F  = self.F==f
                IJ = self.IJ[F]
            except:
                raise csBinError('Indices (F,J,I) not yet set.')

            np.add.at(nObs, IJ, 1)
            np.add.at(aObs, IJ, obs[F])
            
            # Normalize
            # ---------
            gObs = MAPL_UNDEF * np.ones(self.lon_c[0].size, dtype=float)
            K = nObs>0
            gObs[K] = aObs[K] / nObs[K]
            
            csObs[f] = np.reshape(gObs, (self.Ydim, self.Xdim))
        # All done
        # --------
        return csObs    

    def checkGrid(self):
        """
        Verify the grid is valid by two steps: 
          1)Calculate the lon a corner and compare it with the value in the file
          2)Calculate the lats on an edge and comapre them with the values in the file
        If the values agree, the grid prpbably is fine
        """
        if 'corner_lats' in self.ds and 'corner_lons' in self.ds:
           # 1) calculate the corner_lons and corner_lats
           shift     = 0.174532925199433
           tolerance = 10**(-4) # on -180 -- 180 scale
           alpha     = 0.615479708670387
           dalpha    = 2.0*alpha/self.Xdim
           lon_calculated = (1.750*np.pi - shift)/np.pi*180
           lon_in_file    = self.ds['corner_lons'].values[0,0,0]
           # 2) make sure the grid is rotated
           assert abs(lon_calculated-lon_in_file) < tolerance, "Error: Grid should have pi/18 Japan Mount shift"

           J = np.arange(self.Ydim+1)
           lats_calculated = (-alpha+J*dalpha)/np.pi*180
           lats_in_file    = self.ds['corner_lats'].values[0,:,0]
           # 3) compare calculated_lats  and lats in the file, make sure they are the same
           assert all(abs(lats_calculated-lats_in_file) < tolerance),  "Error: cannot handle this grid"
           print("The grid in the file seems fine")
        else:
           print('Not enough information to verify the grid')
            
#...........................................................................
if __name__ == '__main__':

   dyv2_dn  = '/css/g5nr/DYAMONDv2/'
   const_fn = dyv2_dn + '03KM/DYAMONDv2_c2880_L181/const_2d_asm_Mx/202002/DYAMONDv2_c2880_L181.const_2d_asm_Mx.20200201_0000z.nc4'
   csbin = CSBIN(const_fn)
   csbin.checkGrid()
   csbin.set_Indices(csbin.lon_c, csbin.lat_c)
   obs = csbin.binObs(csbin.ds['PHIS'][0].values)
   for f in range(6):
     F = np.reshape(csbin.F, (6, csbin.Ydim, csbin.Xdim))
     assert (F[f,:,:] == f).all(), "face is wrong"
   for j in range(csbin.Ydim):
     J = np.reshape(csbin.J, (6, csbin.Ydim, csbin.Xdim))
     assert (J[:,j,:] == j).all(), "Ydim is wrong"
   for i in range(csbin.Xdim):
     I = np.reshape(csbin.I, (6, csbin.Ydim, csbin.Xdim))
     assert (I[:,:,i] == i).all(), "Xdim is wrong"
   print("The coordinates' indices are consistent with the indices calculated by set_Indices")
