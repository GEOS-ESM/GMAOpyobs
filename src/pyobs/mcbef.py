"""

    Implements base class for McBEF products.

"""

import xarray as xr
import numpy as np
from glob import glob

from . import sampler as sp

class MCBEF(object):

    def __init__(self,paths,**kwargs):
        """
        Loads McBEF properties for one or more days.
        """

        self.ds = xr.open_mfdataset(paths,drop_variables=('crs',),combine='nested',concat_dim='fire',**kwargs)
        self.ds.coords['FP_Time'] = self.ds.FP_Time

    def sample(self, dataset, Variables=None, parallel=True, method='linear'):
        """
        Sample variables on the time and location of each file.
        """

        trj = sp.TRAJECTORY(self.ds.FP_Time.values, 
                            self.ds.FP_Longitude.values, 
                            self.ds.FP_Latitude.values, 
                            dataset, parallel=parallel)

        return trj.sample(Variables=Variables,method=method)



#-------------------------------------------------------------------

if __name__ == "__main__":

    ds = MCBEF('/Users/adasilva/data/mcbef/2019/VNP47MCBEF.State.A201900*.nc')

    
