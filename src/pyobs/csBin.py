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
        
        # Coner coordinates are available as well.
        
    def set_Indices ( self, lons, lats):
        """
        Given a list of longitude and latitudes in (lons,lats), store
        corresponding cube sphere face F, and indices of the grid box
        where each observation falls into. From these indices, multiple
        variables can now be binned.
        """
        
        self.F = np.zeros(lons.shape, dtype=int)
        self.J = np.zeros(lons.shape, dtype=int)
        self.I = np.zeros(lons.shape, dtype=int)
        
        # Calculate (F,J,I)
        
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
            
            
