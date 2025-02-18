"""
   Handling of IODA NetCDF files using Xarray DataTrees.
"""

import numpy  as np
import xarray as xr


def open_ioda(filename, **kwargs):
    """
    Loads an ioda file, tweaking the coordinates   
    """

    dt = xr.open_datatree(filename, **kwargs)

    m = dt['/MetaData']
    
    dt.coords['time'] = m.dateTime
    dt.coords['lat']  = m.latitude
    dt.coords['lon']  = m.longitude
        
    return dt

    
    
         