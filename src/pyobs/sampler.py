"""

   Classes implementing station and trajectgory samplers.

"""

import os

import numpy  as np 
import xarray as xr
import pandas as pd

from datetime import datetime, timedelta
from dateutil.parser import parse as isoparser
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
        

#......................................  Station Sampler CLI ..........................................

def stn_sampler():
    
    """
    Parses command line and write files with resulting station sampling results. 
    """

    from optparse        import OptionParser

    format = 'NETCDF4'
    outFile = 'stn_sampler.nc'
    method = 'linear'
    
#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] stnFile.csv inDataset [iso_t1 iso_t2]\n"+\
                                "where: \n"+
                                "   stnFile.csv          comma separated file with (iso_time,lon,lat)\n"+\
                                "   inDataset            GrADS-style ctl or a shell-style wildcard string\n"+\
                                "   iso_t1,iso_t2        optional beginning and ending time (ISO format)",
                          version='3.0.0' )

    parser.add_option("-o", "--output", dest="outFile", default=outFile,
              help="Output NetCDF file name (default=%s)"\
                          %outFile )

    parser.add_option("-a", "--algorithm", dest="method", default=method,
              help="Interpolation algorithm, one of linear, nearest (default=%s)"\
                          %method)

    parser.add_option("-V", "--vars", dest="Vars", default=None,
              help="Variables to sample, comma delimited (default=All)")
    
    parser.add_option("-f", "--format", dest="format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_64BIT,NETCDF3_CLASSIC (default=%s)"%format )

    #parser.add_option("-I", "--isoTime",
    #                  action="store_true", dest="isoTime",
    #                  help="Include time in ISO format as well.")

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")
    
    (options, args) = parser.parse_args()
    
    if len(args) == 4 :
        stnFile, dataset, iso_t1, iso_t2 = args
        t1, t2 = (isoparser(iso_t1), isoparser(iso_t2))
    elif len(args) == 2 :
        stnFile, dataset = args
        t1, t2 = (None,None)
    else:
        parser.error("must have 2 or 4 arguments: stnFile inDataset [iso_t1 iso_t2]")

    if options.Vars is not None:
        options.Vars = options.Vars.split(',')

    if options.format not in ["NETCDF4","NETCDF4_CLASSIC","NETCDF3_64BIT","NETCDF3_CLASSIC"]:
        raise ValueError('Invalid format <%s>'%options.format)
        
    # Read coordinates from CSV file
    # ------------------------------
    df = pd.read_csv(stnFile, index_col=0)
    stations = df.index.values
    lons = df['lons'].values
    lats = df['lats'].values

    # Sample variables at station locations
    # -------------------------------------
    stn = STATION(stations,lons,lats,dataset,verbose=options.verbose)
    ds = stn.sample(Variables=options.Vars,method=method)
    if options.verbose:
        print(ds)
        print('- Writing',options.outFile)

        
    # Write out netcdf file
    # ---------------------
    ds.to_netcdf(options.outFile,format=options.format)

#......................................  Trajectory Sampler CLI ..........................................

def _getTrackTLE(tleFile,t1,t2,dt):
    """
    Get trajectory from TLE (.tle) file. It is assumed only 1 satellite per file.
    """
    from .tle import TLE
    time, lon, lat = TLE(tleFile).getSubpoint(t1,t2,dt)
    return (lon, lat, time)

def _getTrackICT(ictFile,dt_secs):
    """
    Get trajectory from ICART (.ict) file.
    """
    from .icartt import ICARTT
    m = ICARTT(ictFile)
    lon, lat, tyme = m.Nav['Longitude'], m.Nav['Latitude'], m.Nav['Time']
    mdt = (tyme[-1] - tyme[0]).total_seconds()/float(len(tyme)-1) # in seconds
    idt = int(dt_secs/mdt+0.5)
    return (lon[::idt], lat[::idt], tyme[::idt])

def _getTrackHSRL(hsrlFile,dt_secs=60):
    """
    Get trajectory from HSRL HDF-5 file.
    """
    from .hsrl import HSRL
    h = HSRL(hsrlFile,Nav_only=True)
    lon, lat, tyme = h.lon[:].ravel(), h.lat[:].ravel(), h.tyme[:].ravel()
    if dt_secs > 0:
        dt = tyme[1] - tyme[0] 
        idt = int(dt_secs/dt.total_seconds()+0.5)
        return (lon[::idt], lat[::idt], tyme[::idt])
    else:
        idt = 1
    return (lon[::idt], lat[::idt], tyme[::idt])

def _getTrackCSV(csvFile):
    """
    Get trajectory from a CSV with (lon,lat,time) coordinates.
    """
    df = pd.read_csv(csvFile, index_col=0)
    lon, lat, time = (df['lon'].values,df['lat'].values,pd.to_datetime(df.index).values)
    return (lon,lat,time)   

        
    return ( np.array(lon), np.array(lat), np.array(tyme) )
    
def _getTrackNPZ(npzFile):
    """
    Get trajectory from a NPZ with (lon,lat,time) coordinates.
    Notice that *time* is a datetime object.

    Note: These are simple NPZ usually generated during Neural
          Net or other type of python based utility. Not meant
          for general consumption, but could be since NPZ files
          are much more compact than CSV.

    """
    from .npz import NPZ
    n = NPZ(npzFile)
    if 'time' in n.__dict__:
        return ( n.lon, n.lat, n.time)
    elif 'tyme' in n.__dict__:
        return ( n.lon, n.lat, n.tyme)
    else:
        raise ValueError('NPZ file has neither *time* nor *tyme* attribute.')

def trj_sampler():
    
    """
    Parses command line and write files with resulting trajectory sampling results. 
    """

    from .waypoint import WAYPOINT
    from optparse        import OptionParser

    format = 'NETCDF4'
    rcFile  = 'trj_sampler.rc'
    outFile = 'trj_sampler.nc'
    dt_secs = 60
    method = 'linear'
    plane = 'DC8'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] trjFile inDataset [iso_t1 iso_t2]|[takeOff_isoLocalTime(s)]\n"+\
                                "where: \n"+
                                "   trjFile              Trajecotry specification (time,lon,lat). One of these\n"+\
                                "                        - csvFile        comman separated file\n"+\
                                "                        - wpFile         waypoint file; in this case t1,t2,dt are \n"+\
                                "                                         takeoff times\n"+\
                                "                        - tleFile        two line element file (1 sat per file)\n"+\
                                "                        - ictFile        ICARTT format file\n"+\
                                "                        - npzFile        Numpy NPZ file\n"+\
                                "   inDataset            GrADS-style ctl or a shell-style wildcard string\n"+\
                                "   iso_t1,iso_t2        optional beginning and ending time (ISO format)",
                          version='1.0.1' )

    parser.add_option("-a", "--algorithm", dest="method", default=method,
              help="Interpolation algorithm, one of linear, nearest (default=%s)"\
                          %method)

    parser.add_option("-o", "--output", dest="outFile", default=outFile,
              help="Output NetCDF file (default=%s)"\
                          %outFile )

    parser.add_option("-f", "--format", dest="format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC or NETCDF3_64BIT (default=%s)"%format )

    parser.add_option("-p", "--plane", dest="plane", default='DC8',
              help="aircraft: DC8, ER2, ... or 'snapshot' (default=%s)"%plane )

    parser.add_option("-V", "--vars", dest="Vars", default=None,
              help="Variables to sample, comma delimited (default=All)")
    
    parser.add_option("-t", "--trajectory", dest="traj", default=None,
                      help="Trajectory file format: one of tle, ict, csv, wp, npz (default=trjFile extension except for wp)" )

    parser.add_option("-d", "--dt_secs", dest="dt_secs", default=dt_secs,
              type='int',
              help="Timesetp in seconds for TLE sampling (default=%s)"%dt_secs )

    #parser.add_option("-I", "--isoTime",
    #                  action="store_true", dest="isoTime",
    #                  help="Include ISO format time in output file.")

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")

    (options, args) = parser.parse_args()
    
    if options.traj == 'WP':
        trjFile, dataset = args[0:2]
        TakeOff = args[2:]
    elif len(args) == 4:
        trjFile, dataset, iso_t1, iso_t2 = args
        t1, t2 = (isoparser(iso_t1), isoparser(iso_t2))
        dt = timedelta(seconds=options.dt_secs)
    elif len(args) == 2:
        trjFile, dataset = args
        t1, t2 = None, None
    else:
        parser.error("must have 2 or 4 arguments: tleFile|ictFile [iso_t1 iso_t2]")

    if options.traj is None:
        name, ext = os.path.splitext(trjFile)
        options.traj = ext[1:]
    options.traj = options.traj.upper()
        
    # Create consistent file name extension
    # -------------------------------------
    name, ext = os.path.splitext(options.outFile)
    if 'NETCDF4' in options.format:
        options.outFile = name + '.nc4'
    elif 'NETCDF3' in options.format:
        options.outFile = name + '.nc'
    else:
        raise ValueError('Invalid extension <%s>'%ext)

    # Create trajectory
    # -----------------
    if options.traj == 'TLE':
        if t1 is None:
            raise ValueError('time range (t1,t2) must be specified when doing TLE sampling.')
        lon, lat, time = _getTrackTLE(trjFile, t1, t2, dt)
    elif options.traj == 'ICT':
        lon, lat, time = _getTrackICT(trjFile,options.dt_secs)
    elif options.traj == 'CSV':
        lon, lat, time = _getTrackCSV(trjFile)
    elif options.traj == 'WP':
        pass # special handling
    elif options.traj == 'NPZ':
        lon, lat, time = _getTrackNPZ(trjFile)
    elif options.traj == 'HSRL' or options.traj == 'H5': # deprecated, undocumented for now
        lon, lat, time = _getTrackHSRL(trjFile,options.dt_secs)
    else:
        raise ValueError('cannot handle trajectory file format <%s>'%options.traj)


    # Waypoints (several takeoff times)
    # ---------------------------------
    if options.traj == 'WP':
        name, ext = os.path.splitext(options.outFile) # prepare to append to name
        outFile = name + '.@city_@aircraft_@takeoff' + ext # template for addition
        wp = WAYPOINT(trjFile, options.plane, verbose=options.verbose)
        for takeoff in TakeOff:
            outFile_ = outFile.replace('@city',wp.city).\
                               replace('@aircraft',wp.plane).\
                               replace('@takeoff',str(takeoff).replace(' ','T'))
            df = wp.getTraj(takeoff)
            time = df.index.values
            lon = df['lon'].values
            lat = df['lat'].values
            trj = TRAJECTORY(time,lon,lat,dataset,verbose=options.verbose)
            ds = trj.sample(Variables=options.Vars,method=method)

            if options.verbose:
                print('- Writing',outFile,'from',trjFile,'at takeoff',takeoff)

            ds.to_netcdf(outFile_,format=options.format,compute=True)

    # All else
    # --------
    else:
        
        trj = TRAJECTORY(time,lon,lat,dataset,verbose=options.verbose)
        ds = trj.sample(Variables=options.Vars,method=method)
        if options.verbose:
            #print(ds)
            print('- Writing',outFile,'from',trjFile,'(%s)'%options.traj)
        
        # Write out netcdf file
        # ---------------------
        ds.to_netcdf(options.outFile,format=options.format)

#...................................... Simple Minded Testing ..........................................

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


        
