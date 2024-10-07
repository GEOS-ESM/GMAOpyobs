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
    
    def __init__(self, **kwargs):
        """
        Supporting keywords: xdim, stretch_factor, target_lon, target_lat, and filename.
        If Filename is given, all the other keywords are not necessary and ignored
        """

        kw_ = {k.upper():v for k,v in kwargs.items()}        

        self.stretched = False
        if 'FILENAME' in kw_.keys():
           with xr.open_dataset(kw_['FILENAME']) as ds:
              self.Xdim  = ds.dims['Xdim']
              self.stretch_factor = 1.0
              self.target_lon     = 0.0
              self.target_lat     = -np.pi/2.0

              if 'STRETCH_FACTOR' in ds.attrs :
                 self.stretch_factor = ds.attrs['STRETCH_FACTOR']
                 self.target_lon     = ds.attrs['TARGET_LON']/180*np.pi
                 self.target_lat     = ds.attrs['TARGET_LAT']/180*np.pi
        else:
           self.Xdim = kw_['XDIM']
           self.stretch_factor = kw_.get('STRETCH_FACTOR',    1.0)
           self.target_lon     = kw_.get('TARGET_LON',        0.0)
           self.target_lat     = kw_.get('TARGET_LAT', -np.pi/2.0)

        if (abs(self.stretch_factor-1) > 1.0e-4):
           self.stretched    = True
                
    @staticmethod            
    def _reverse_rotate(lon_p, lat_p, lons_, lats_):
      """
       Rotate (lons_, lats_) back to its original position before Schmidt stretch
      """
    
      half_pi= np.pi/2
      two_pi = 2*np.pi
      x = np.cos(lats_)*np.cos(lons_-lon_p)
      y = np.cos(lats_)*np.sin(lons_-lon_p)
      z = np.sin(lats_)
      X = np.sin(lat_p)*x - np.cos(lat_p)*z
      Y = -y
      Z = -np.cos(lat_p)*x - np.sin(lat_p)*z
    
      n_s = (1. - abs(Z)) < 10**(-7)
    
      lons  = np.where(n_s, 0, np.arctan2(Y,X))
      lats  = np.where(n_s, half_pi*np.copysign(1,Z), np.arcsin(Z))
    
      lons = np.where(lons < 0,       lons + two_pi, lons)
      lons = np.where(lons >= two_pi, lons - two_pi, lons)
    
      return lons, lats

    @staticmethod 
    def _reverse_schmidt(c, lon_p, lat_p, lons_, lats_):
       """
         Perform reverse Schmidt transform
         c: stretch factor
         lon_p, lat_p : target point
         lons_, lats_ : coordinate after Schmidt transform
       """
       c2p1 = 1 + c*c
       c2m1 = 1 - c*c
    
       # 1) rotate back to the pole
       lons, lats = CSBIN._reverse_rotate(lon_p, lat_p, lons_, lats_)
    
       #  2) reverse stretch
       if (abs(c2m1) > 10**(-7)): 
         lats  = np.arcsin( (c2m1-c2p1*np.sin(lats))/(c2m1*np.sin(lats)-c2p1))
    
       return lons, lats

    @staticmethod
    def _checkGridInFile(filename):
        """
        Verify the grid in a file is valid by two steps: 
          1)Calculate the lon a corner and compare it with the value in the file
          2)Calculate the lats on an edge and comapre them with the values in the file
        If the values agree, the grid prpbably is fine
        """
        
        ds   = xr.open_dataset(filename)
        Xdim = ds.dims['Xdim']

        if 'corner_lats' in ds and 'corner_lons' in ds:
           # 1) calculate the corner_lons and corner_lats
           shift     = 0.174532925199433
           tolerance = 10**(-5) # on 0--pi
           alpha     = 0.615479708670387
           dalpha    = 2.0*alpha/Xdim
           dimC = ds['corner_lons'][0,:,0].size
           lons_calculated = np.ones(dimC, dtype=float)*(1.750*np.pi - shift)

           lons = ds['corner_lons'].values[0,:,0]/180*np.pi
           lats = ds['corner_lats'].values[0,:,0]/180*np.pi
           stretch_factor = 1.0
           if 'STRETCH_FACTOR' in ds.attrs :
              stretch_factor = ds.attrs['STRETCH_FACTOR']
              target_lon     = ds.attrs['TARGET_LON']/180*np.pi
              target_lat     = ds.attrs['TARGET_LAT']/180*np.pi
           stretched = False
           if (abs(stretch_factor-1) > 1.0e-4):
             stretched    = True

           if stretched :
              lons_calculated += shift # add shift back to make sure no shift for stretched grid
              lons, lats = CSBIN._reverse_schmidt(stretch_factor, target_lon, target_lat, lons , lats)

           J = np.arange(dimC)
           lats_calculated = (-alpha+J*dalpha)
           # 2) compare calculated corners with corners in the file, make sure they are the same
           assert all(abs(lats_calculated-lats) < tolerance),  "Error lats: cannot handle this grid"
           assert all(abs(lons_calculated-lons) < tolerance),  "Error lons: cannot handle this grid"
           print("The grid in the file {} seems fine".format(filename))
        else:
           print('Not enough information to verify the grid')

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
        lons  = lons/180*np.pi
        lats  = lats/180*np.pi

        if self.stretched:
           lons, lats = CSBIN._reverse_schmidt(self.stretch_factor, self.target_lon, self.target_lat, lons, lats)
           lons  -=  shift

        # shift the grid away from Japan Fuji Mt.
        lons_  = lons + shift
        lats_  = lats
        # some functions only work for 1-d array. ravel it first
        if len(lons.shape) > 1:
          lons_  = np.ravel(lons_)
          lats_  = np.ravel(lats_)
           
        self.F = np.zeros(lons_.size, dtype=int)
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
        self.F[f1Mask] = 0
        I[f1Mask] = np.floor((np.arctan( y[f1Mask]/sqr2)+alpha)/dalpha)
        J[f1Mask] = np.floor((np.arctan( z[f1Mask]/sqr2)+alpha)/dalpha)

        f2Mask = (1-y) <= tolerance
        self.F[f2Mask] = 1
        I[f2Mask] = np.floor((np.arctan(-x[f2Mask]/sqr2)+alpha)/dalpha)
        J[f2Mask] = np.floor((np.arctan( z[f2Mask]/sqr2)+alpha)/dalpha)

        f3Mask = (1-z) <= tolerance
        self.F[f3Mask] = 2
        I[f3Mask] = np.floor((np.arctan(-x[f3Mask]/sqr2)+alpha)/dalpha)
        J[f3Mask] = np.floor((np.arctan(-y[f3Mask]/sqr2)+alpha)/dalpha)

        f4Mask = (x+1.0) <= tolerance
        self.F[f4Mask] = 3
        I[f4Mask] = np.floor((np.arctan(-z[f4Mask]/sqr2)+alpha)/dalpha)
        J[f4Mask] = np.floor((np.arctan(-y[f4Mask]/sqr2)+alpha)/dalpha)

        f5Mask = (y+1.0) <= tolerance
        self.F[f5Mask] = 4
        I[f5Mask] = np.floor((np.arctan(-z[f5Mask]/sqr2)+alpha)/dalpha)
        J[f5Mask] = np.floor((np.arctan( x[f5Mask]/sqr2)+alpha)/dalpha)

        f6Mask = (z+1.0) <= tolerance
        self.F[f6Mask] = 5
        I[f6Mask] = np.floor((np.arctan( y[f6Mask]/sqr2)+alpha)/dalpha)
        J[f6Mask] = np.floor((np.arctan( x[f6Mask]/sqr2)+alpha)/dalpha)

        self.IJ = I + J*self.Xdim

        return 
    
    def binObs ( self, obs ):
        """
        Given a list of longitude and latitudes in (lons,lats),
        bin list of observations *obs* on the cubed sphere,
        returning gridded observations *csObs*.
        """
        
        csObs = np.zeros((6, self.Xdim, self.Xdim), dtype=float)
       
        if (len(obs.shape) > 1):
           obs = np.ravel(obs) 
        # One face at a time
        # ------------------
        face_size = self.Xdim*self.Xdim

        for f in range(6):
            # Accumulator, counter
            # --------------------
            aObs = np.zeros(face_size, dtype=float)
            nObs = np.zeros(face_size, dtype=float)
            
            # Accumulate observations lying on this face of the cube
            # ------------------------------------------------------
            try:
                F  = self.F==f
                IJ = self.IJ[F]
            except:
                raise csBinError('Indices (F,IJ) not yet set.')

            np.add.at(nObs, IJ, 1)
            np.add.at(aObs, IJ, obs[F])
            
            # Normalize
            # ---------
            gObs = MAPL_UNDEF * np.ones(face_size, dtype=float)
            K = nObs>0
            gObs[K] = aObs[K] / nObs[K]
            
            csObs[f] = np.reshape(gObs, (self.Xdim, self.Xdim))
        # All done
        # --------
        return csObs    




#...........................................................................
if __name__ == '__main__':

   dyv2_dn  = '/css/g5nr/DYAMONDv2/'
   const_fn = dyv2_dn + '03KM/DYAMONDv2_c2880_L181/const_2d_asm_Mx/202002/DYAMONDv2_c2880_L181.const_2d_asm_Mx.20200201_0000z.nc4'
   #const_fn = "./v12-stock-2024Sep26-1day-c540-NoPoints-RunTimeFix.geosgcm_prog_nat.20200415_0600z.nc4"
   csbin = CSBIN(filename = const_fn)

   CSBIN._checkGridInFile(const_fn)

   ds   = xr.open_dataset( const_fn)
   lon_c = ds['lons']
   lat_c = ds['lats']
   Xdim  = ds.dims['Xdim'] 
   csbin.set_Indices(lon_c, lat_c)
   #obs = csbin.binObs(csbin.ds['PHIS'][0].values)
   ij =np.array([i+j*Xdim for j in range(Xdim) for i in range(Xdim)]).reshape((Xdim,Xdim))

   for f in range(6):
     F = np.reshape(csbin.F, (6, csbin.Xdim, csbin.Xdim))
     assert (F[f,:,:]  == f).all(), "face is wrong"

     face = csbin.F == f
     IJ = np.reshape(csbin.IJ[face], (csbin.Xdim, csbin.Xdim))
     assert (IJ == ij).all(), " IJ index is wrong"

   print("The coordinates' indices are consistent with the indices calculated by set_Indices")
