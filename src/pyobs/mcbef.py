"""

    Implements base class for McBEF products.

"""

import xarray as xr
import numpy as np
from glob import glob

from . import sampler as sp

QA = dict (Uniphasic                           = 1,
           Biphasic                            = 2,
           Biphasic2Uniphasic                  = 3,
           UniphasicWithCLTBackground          = 11,
           BiphasicWithCLTBackground           = 12, 
           Biphasic2UniphasicWithCLTBackground = 13,
           Bowtie                              = 100, 
           BadBackgroundSignal                 = 101, 
           BadFireSignal                       = 102, 
           FailUniphasic                       = 103,
           FailBiphasic                        = 104,
          )

STATS = dict ( lower=0, mode=1, upper=2, mean=3, stdv=4 )

class McBEF_Dataset(xr.Dataset):
    """
    Subclassing xarray Dataset with some McBEF specific methods.
    """
    __slots__ = ()
    def Hello(self):
        print("Hello, world!")
        
        
    def Sample(self, dataset, Variables=None, parallel=True, method='linear'):
        """
        Sample variables on the time and location of each file.
        """

        trj = sp.TRAJECTORY(self.FP_Time.values, 
                            self.FP_Longitude.values, 
                            self.FP_Latitude.values, 
                            dataset, parallel=parallel)

        return trj.sample(Variables=Variables,method=method)
    

class MCBEF(object): # Deprecated

    def __init__(self,paths,**kwargs):
        """
        IMPORTANT: this is being debrecated.
        
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

def open_dataset(*args, 
                 OnlyBiphasic=False,
                 OnlyUniphasic=False,
                 OnlyGood=False,
                 OnlyBad=False,
                 **kwargs):
    """
    Wrapper of xr.open_dataset(), with some optional filtering:

    OnlyBiphasic  --- only fires with QA_flag = 2
    OnlyUniphasic --- only fires with QA_flag = 1
    OnlyGood      --- only fires with QA_flag < 100
    OnlyBad       --- only fires with QA_fkag >= 100
        
    """
    ds = xr.open_dataset(*args,**kwargs)
    ds.__class__ = McBEF_Dataset

    if OnlyBiphasic:
        I = ds.QA_flag==QA['Biphasic']
    elif OnlyUniphasic:
        I = ds.QA_flag==QA['Uniphasic']
    elif OnlyGood:
        I = ds.QA_flag<QA['Bowtie']
    elif OnlyBad:
        I = ds.QA_flag>=QA['Bowtie']
    else:
        return ds  # no filtering

    J = np.arange(ds.dims['fire'],dtype=int)[I] # fires to retain.
    
    return ds.isel(fire=J)
    
def open_mfdataset(*args, **kwargs):
    """
    Wrapper of xr.open_mfdataset().
    """
    ds = xr.open_mfdataset(*args,**kwargs)
    ds.__class__ = McBEF_Dataset
    return ds

#-------------------------------------------------------------------

if __name__ == "__main__":

    ds = MCBEF('/Users/adasilva/data/mcbef/2019/VNP47MCBEF.State.A201900*.nc')

    
