"""

   Classes implementing station and trajectgory samplers.

"""

import os

import numpy  as np 
import xarray as xr
import pandas as pd

from datetime import datetime, timedelta
from glob import glob

from . import xrctl as xc

os.environ['HDF5_USE_FILE_LOCKING']='FALSE'

#............................................................
class SamplerError(Exception):
    """
    Defines NC4ctl general exception errors.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class STATION(object):

    def __init__(self, stations, lons, lats,
                 dataset, time_range=None, verbose=False):
        """
        Specifies dataset to be sampled at obs location.
        On input,

        stations: station names (labels)
        lons, lats: cooridnates of each station
        dataset: the input dataset, it can be one of these
                 xr.Dataset: an xarray dataset
                 string : either a GrASDS-style control file
                          (must have extension .ctl or .xdf)
                          or a glob template (e.g., *.nc)
                 list,tuple: a list of file names
        time_range: when using a GrADS templates, the time interval
              to generate a list of files.

        """

        self.verb = verbose
        
        # If dataset is an xarray dataset we are good to go
        # -------------------------------------------------
        if isinstance(dataset,xr.Dataset):
            self.ds = dataset # we are good to go...

        # If dataset is a list of files...
        # --------------------------------
        elif isinstance(dataset,(list,tuple)):
            self.ds = xr.open_mfdataset(dataset,parallel=True)

        # If datatset is a string it is either a GrADS-style ctl or
        # a glob type of template
        # ---------------------------------------------------------
        elif isinstance(dataset,str):
            # Special handles GrADS-style ctl if found
            self.ds = xc.open_mfdataset(dataset,time_range=time_range,parallel=True)

        else:
            raise SamplerError("Invalid dataset specification.")
            
        # Save coordinates
        # ----------------
        self.stations = xr.DataArray(stations, dims='station')
        self.lons = xr.DataArray(lons, dims='station',attrs=self.ds.coords['lon'].attrs)
        self.lats = xr.DataArray(lats, dims='station',attrs=self.ds.coords['lat'].attrs)

        # TO DO: when using xESMF for regridding, pre-compute transforms here
        # -------------------------------------------------------------------

        
    #--
    def sample(self,Variables=None,method='linear'):
        """
        Sample variables and pre-determined obs locations

        """
        if Variables is None:
            Variables = list(self.ds.data_vars)

        elif isinstance(Variables,str):
            Variables = [Variables,]

        sampled = dict()

        for vn in Variables:
            if self.verb: print('[ ] sampling ',vn)
            sampled[vn] = self.ds[vn].interp(lon=self.lons,lat=self.lats,method=method)

        return xr.Dataset(sampled).assign_coords({'station': self.stations})

#......................................................................................

class TRAJECTORY(object):

    def __init__(self, times, lons, lats, dataset, verbose=False):
        """
        Specifies dataset to be sampled at obs location.
        On input,

        times, lons, lats: trajectory coordinates
        dataset: the input dataset, it can be one of these
                 xr.Dataset: an xarray dataset
                 string : either a GrASDS-style control file
                          (must have extension .ctl or .xdf)
                          or a glob template (e.g., *.nc)
                 list,tuple: a list of file names
        time_range: when using a GrADS templates, the time interval
              to generate a list of files.

        """

        self.verb = verbose
        self.times = xr.DataArray(times,dims='time')
        time_range = times.min(), times.max()
        
        # If dataset is an xarray dataset we are good to go
        # -------------------------------------------------
        if isinstance(dataset,xr.Dataset):
            self.ds = dataset # we are good to go...

        # If dataset is a list of files...
        # --------------------------------
        elif isinstance(dataset,(list,tuple)):
            self.ds = xr.open_mfdataset(dataset,parallel=True)

        # If datatset is a string it is either a GrADS-style ctl or
        # a glob type of template
        # ---------------------------------------------------------
        elif isinstance(dataset,str):
            # Special handles GrADS-style ctl if found
            # ----------------------------------------
            self.ds = xc.open_mfdataset(dataset,parallel=True) # special handles GrADS-style ctl if found

        else:
            raise SamplerError("Invalid dataset specification.")

        # Save coordinates with proper attributes
        # ---------------------------------------
        self.lons = xr.DataArray(lons, dims='time',attrs=self.ds.coords['lon'].attrs)
        self.lats = xr.DataArray(lats, dims='time',attrs=self.ds.coords['lat'].attrs)

        # TO DO: when using xESMF for regridding, pre-compute transforms here
        # -------------------------------------------------------------------

    #--
    def sample(self,Variables=None,method='linear'):
        """
        Sample variables and pre-determined obs locations

        """
        if Variables is None:
            Variables = list(self.ds.data_vars)

        elif isinstance(Variables,str):
            Variables = [Variables,]

        sampled = dict()

        for vn in Variables:
            if self.verb: print('[ ] sampling',vn)
            sampled[vn] = self.ds[vn].interp(time=self.times,lon=self.lons,lat=self.lats,method=method)

        return xr.Dataset(sampled).assign_coords({'time': self.times})

#......................................................................................

class TLETRAJ(TRAJECTORY):


    def __init__ (self, tleFile, t1, t2, dt, *args, **kwargs):
        """
        Generate trajectory from Two-line (TLE) file.

        t1, t2: datetime, time interval
        dt    : timedelta, timestep

        """

        from .tle import TLE
        
        # Generate coordinates
        # --------------------
        times, lons, lats = TLE(tleFile).getSubpoint(t1,t2,dt)

        # Initialize base class
        # ---------------------
        super().__init__(times, lons, lats, *args, **kwargs)
        

class WPTRAJ(TRAJECTORY):

    def __init__ (self, wpFile, plane, takeoff, *args, **kwargs):
        """
        Generate trajectory from a CSV waypoint file.

        """

        from .waypoint import WAYPOINT

        # Generate trajectory from waypoint file and takeoff time
        # -------------------------------------------------------
        traj = WAYPOINT(wpFile, plane).getTraj(takeoff)

        # Initialize base class
        # ---------------------
        times, lons, lats = traj.index.values, traj['lon'].values, traj['lat'].values 
        super().__init__(times, lons, lats, *args, **kwargs)
        

#......................................................................................

if __name__ == "__main__":

      pass
  
def test_tle():

      tleFile = '/Users/adasilva/data/tle/terra/terra.2023-04-15.tle'

      aer_Nx = '/Users/adasilva/data/merra2/ctl/tavg1_2d_aer_Nx.ctl' # GrADSctl

      t1 = datetime(2023,4,15,0,0,0)
      t2 = datetime(2023,4,15,6,0,0)
      dt = timedelta(minutes=1)
    
      wt = TLETRAJ(tleFile,t1,t2,dt,aer_Nx,verbose=True)

      ds = wt.sample()

      return ds
  
def test_waypoint():

      wpFile = '/Users/adasilva/data/wp/phillipines_waypoints.csv'

      aer_Nx = '/Users/adasilva/data/merra2/ctl/tavg1_2d_aer_Nx.ctl' # GrADSctl

      takeoff = '2023-04-15T08:00:00'       # either string or datetime
      takeoff = datetime(2023,4,15,8,0,0)
      
      wt = WPTRAJ(wpFile,'DC8',takeoff,aer_Nx,verbose=True)
      
      ds = wt.sample()

      return ds
  
def test_trajecgory():
    
      from datetime import datetime
    
      merra2_dn = '/Users/adasilva/data/merra2/Y2023/M04/'
      aer_Nx = merra2_dn + '/MERRA2.tavg1_2d_aer_Nx.????????.nc4'

      traj_fn = '/Users/adasilva/data/merra2/DC8_20230426.nc'

      c = xr.open_dataset(traj_fn)

      times, lons, lats = c['time'].values, c['lon'].values, c['lat'].values
      
      traj = TRAJECTORY(times, lons, lats, aer_Nx)
      ds = traj.sample(Variables=['DUEXTTAU', 'DUCMASS'])

      print(ds)
      
def test_stations():
      
      fluxnet_fn = '/Users/adasilva/data/brdf/fluxnet_stations.csv'

      stations = pd.read_csv('/Users/adasilva/data/brdf/fluxnet_stations.csv',
                             index_col=0)

      print(stations)

      lons = stations['lons'].values
      lats = stations['lats'].values


      # Using file lists
      # ----------------
      stn = STATION(stations.index,lons,lats,aer_Nx,verbose=1)
      ds = stn.sample(Variables=['DUEXTTAU', 'DUCMASS'])
      print(ds)

      # GrADS-style ctl
      # ---------------
      ctlfile = '/Users/adasilva/data/merra2/ctl/tavg1_2d_aer_Nx.ctl'
      tbeg, tend = datetime(2023,4,7,0,30), datetime(2023,4,15,23,30)
      stn2 = STATION(stations.index,lons,lats,ctlfile,
                     time_range=(tbeg,tend),verbose=1)
      ds2 = stn2.sample(Variables=['DUEXTTAU', 'DUCMASS'])
      print(ds2)


        
