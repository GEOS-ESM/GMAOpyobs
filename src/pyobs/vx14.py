"""
      VIIRS Level2 Fire Products access via xarray.

"""

import os
import sys

import xarray as xr
from glob import glob
from dateutil.parser import parse as isoparser
      
class VX14(object):
    
    def __init__ (self, path, verbose=False):
        """
        Lazy loads mutiple granules of VIIRS fire products, concatenating
        into a single dataset.
        
        path:    str, list - glob style file pattern or list of fles
        
        """
        
        self.verb = verbose
        
        if isinstance(path,str):
            Files = sorted(glob(path))
        else:
            Files = path
        
        DS = []
        nobs = 0
        for f in Files:
            
            ds = xr.open_dataset(f,drop_variables=('fire mask', 'algorithm QA')) # For now
            
            n_fires = ds['FP_power'].size
            
            if n_fires > 0:
                
               nobs += n_fires
                
               ds = ds.rename({'phony_dim_0':'fire'})
            
               # Time coordinates
               # ---------------- 
               t0, t1 = isoparser(ds.PGE_StartTime), isoparser(ds.PGE_EndTime)
               tm = t0 + (t1-t0)/2 # mid time of granule
               time = xr.DataArray(n_fires*[tm,],dims='fire')
                
               ds = ds.assign_coords({'FP_time':      time, 
                                      'FP_longitude': ds.FP_longitude, 
                                      'FP_latitude':  ds.FP_latitude})
             
               DS += [ds,]
          
        # Concatenate it all in a single dataset
        # --------------------------------------
        self.ds = xr.concat(DS,dim='fire')
        
        if self.verb:
             print('[] Total number of fires:', nobs)
                
    @property
    def traj(self): 
        coords = self.ds.coords
        return (coords['FP_time'],coords['FP_longitude'],coords['FP_latitude'])
        
