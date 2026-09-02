"""

   Classes implementing station and trajectgory samplers.

"""

import os

import numpy  as np
import xarray as xr
import xesmf as xe
import pandas as pd

from datetime import datetime, timedelta
from dateutil.parser import parse as isoparser
from glob import glob

from . import xrctl as xc
import fsspec
import dask.distributed
import esmpy as ESMF
from scipy.spatial import ConvexHull

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

def build_dual_mesh(ds):
    """
    Builds a continuous triangulated mesh using the cell centers as the nodes.
    This allows ESMF to perform Bilinear interpolation from cubed-sphere centers.
    """
    # Grab the cell center coordinates (usually 'lon' and 'lat' with dims: nf, Ydim, Xdim)
    if 'time' in ds['lon'].dims:
        lon = ds['lon'].isel(time=0).values
        lat = ds['lat'].isel(time=0).values
    else:
        lon = ds['lon'].values
        lat = ds['lat'].values

    flat_lon = lon.ravel()
    flat_lat = lat.ravel()
    n_nodes = len(flat_lon)

    # Convert lat/lon to 3D Cartesian coordinates
    lon_rad = np.radians(flat_lon)
    lat_rad = np.radians(flat_lat)
    x = np.cos(lat_rad) * np.cos(lon_rad)
    y = np.cos(lat_rad) * np.sin(lon_rad)
    z = np.sin(lat_rad)
    
    # Use ConvexHull to automatically triangulate the sphere
    coords_3d = np.column_stack((x, y, z))
    hull = ConvexHull(coords_3d)
    
    # hull.simplices contains the 0-based indices for the 3 corners of every triangle
    triangles = hull.simplices

    # Create ESMF Mesh
    mesh = ESMF.Mesh(parametric_dim=2, spatial_dim=2, coord_sys=ESMF.CoordSys.SPH_DEG)
    
    # Add the cell centers as NODES
    node_ids = np.arange(1, n_nodes + 1, dtype=np.int32)
    node_coords = np.empty(2 * n_nodes, dtype=np.float64)
    node_coords[0::2] = flat_lon
    node_coords[1::2] = flat_lat
    node_owners = np.zeros(n_nodes, dtype=np.int32)
    
    mesh.add_nodes(n_nodes, node_ids, node_coords, node_owners)
    
    # Add the Triangles as ELEMENTS
    n_elems = len(triangles)
    elem_ids = np.arange(1, n_elems + 1, dtype=np.int32)
    elem_types = np.full(n_elems, ESMF.MeshElemType.TRI, dtype=np.int32)
    elem_conn = triangles.ravel().astype(np.int32) # ConvexHull is natively 0-based, perfect for esmpy!
    
    mesh.add_elements(n_elems, elem_ids, elem_types, elem_conn)
    
    return mesh

class STATION(object):

    def __init__(self, stations, lons, lats,
                 dataset, time_range=None, times=None, 
                 verbose=False,
                 parallel=True,chunks='auto',cs=False,
                 **kwargs):
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
        times: optional specific times to sample at, 
              otherwise output at model time resolution
        cs: sampling on the cubed sphere
        """

        self.verb = verbose

        # If dataset is an xarray dataset we are good to go
        # -------------------------------------------------
        if isinstance(dataset,xr.Dataset):
            self.ds = dataset # we are good to go...
        elif isinstance(dataset,dask.distributed.Future):
            self.ds = xr.open_dataset(dataset, engine="zarr",chunks=chunks, backend_kwargs={"consolidated": False})
        #  if dataset is a parquet referece file path
        elif dataset[-4:] == 'parq':
            fs = fsspec.filesystem("reference", fo=dataset, remote_protocol='file', lazy=True)
            self.ds = xr.open_dataset(fs.get_mapper(), engine="zarr",chunks=chunks, backend_kwargs={"consolidated": False})        
        # If dataset is a list of files...
        # OR GrADS-style ctl
        # OR a glob type of template
        # --------------------------------
        elif isinstance(dataset,(list,tuple,str)):
            self.ds = xc.open_mfdataset(dataset,time_range=time_range,parallel=parallel,chunks=chunks,cs=cs,**kwargs)

        else:
            raise SamplerError("Invalid dataset specification.")

        # Save coordinates
        # ----------------
        self.stations = xr.DataArray(stations, dims='station')
        self.lons = xr.DataArray(lons, dims='station',attrs=self.ds.coords['lon'].attrs)
        self.lats = xr.DataArray(lats, dims='station',attrs=self.ds.coords['lat'].attrs)

        if times is not None:
            self.times = xr.DataArray(times,dims='time')
        else:
            self.times = times

        # Use xESMF for regridding, pre-compute transforms here
        # -------------------------------------------------------------------
        if cs:
            # Build ESMF mesh
            self.src_mesh = build_dual_mesh(self.ds)

            n_stations = len(lons)

            # Create locstream
            self.dst_locstream = ESMF.LocStream(n_stations, coord_sys=ESMF.CoordSys.SPH_DEG)
            self.dst_locstream["ESMF:Lon"] = np.asarray(lons, dtype=np.float64)
            self.dst_locstream["ESMF:Lat"] = np.asarray(lats, dtype=np.float64)

            self.src_field = ESMF.Field(self.src_mesh, meshloc=ESMF.MeshLoc.NODE)
            self.dst_field = ESMF.Field(self.dst_locstream, name='dst')

            self.regridder = ESMF.Regrid(self.src_field,
                                         self.dst_field,
                                         regrid_method=ESMF.RegridMethod.BILINEAR,
                                         unmapped_action=ESMF.UnmappedAction.IGNORE)
        else:
            ds_loc = xr.Dataset({
                "lon": (["locations"], np.asarray(lons, dtype=np.float64)),
                "lat": (["locations"], np.asarray(lats, dtype=np.float64)),
            })

            self.regridder = xe.Regridder(self.ds, ds_loc, "bilinear", locstream_out=True)

        self.cs = cs

    #--
    def sample(self,Variables=None,method='linear'):
        """
        Sample variables and pre-determined obs locations

        """
        if Variables is None:
            if self.cs:
                # Only grab variables that actually have cubed-sphere spatial dimensions
                Variables = [
                    vn for vn in self.ds.data_vars
                    if all(d in self.ds[vn].dims for d in ('nf', 'Ydim', 'Xdim'))
                ]
            else:
                Variables = list(self.ds.data_vars)

        elif isinstance(Variables,str):
            Variables = [Variables,]

        sampled = dict()

        for vn in Variables:
            if self.verb: print('[ ] sampling ',vn)
            if self.cs:

                da_var = self.ds[vn]
                # Identify extra dimensions (anything not part of the spatial cubed-sphere grid)
                spatial_dims = ('nf', 'Ydim', 'Xdim')
                extra_dims = [d for d in da_var.dims if d not in spatial_dims]
                
                # Reorder the DataArray so extra dimensions come first, spatial dims at the end
                da_var = da_var.transpose(*extra_dims, *spatial_dims)
                
                # Calculate shapes
                extra_shape = tuple(da_var.sizes[d] for d in extra_dims)
                num_slices = int(np.prod(extra_shape)) if extra_shape else 1
                num_elements = int(np.prod([da_var.sizes[d] for d in spatial_dims]))
                num_stations = len(self.stations)
                
                # Reshape data into a 2D array: (total_extra_slices, total_spatial_elements)
                # If there are no extra dims, this just becomes (1, num_elements)
                flat_data_2d = da_var.values.reshape(num_slices, num_elements)
                
                # Pre-allocate array for the interpolated results
                interpolated_results = np.empty((num_slices, num_stations), dtype=da_var.dtype)
                
                # Loop over each slice, regrid, and save
                for i in range(num_slices):
                    # Assign the 1D spatial data slice to the source field
                    self.src_field.data[...] = flat_data_2d[i, :]
                    
                    # Perform the regridding
                    self.regridder(self.src_field, self.dst_field)
                    
                    # Store the result
                    interpolated_results[i, :] = self.dst_field.data
                
                # Reshape results back to (*extra_shape, num_stations)
                final_shape = extra_shape + (num_stations,)
                interpolated_results = interpolated_results.reshape(final_shape)
                
                # Package back into an xarray DataArray
                # Gather coordinates for the extra dimensions + the stations
                out_coords = {d: da_var.coords[d] for d in extra_dims if d in da_var.coords}
                out_coords['station'] = self.stations.values
                
                da_spatial = xr.DataArray(
                    interpolated_results,
                    dims=extra_dims + ['station'],
                    coords=out_coords,
                    name=vn,
                    attrs=da_var.attrs
                )

                sampled[vn] = da_spatial

            else:
                sampled[vn] = self.regridder(self.ds[vn])
                sampled[vn] = sampled[vn].assign_coords(time=('locations', self.stations.values)).rename({'locations': 'station'})

#            sampled[vn] = self.ds[vn].interp(time=self.times,lon=self.lons,lat=self.lats,method=method)

            if self.times is not None:
                sampled[vn] = sampled[vn].interp(time=self.times,method=method)

        return xr.Dataset(sampled)

#......................................................................................

class TRAJECTORY(object):

    def __init__(self, times, lons, lats, dataset, 
                  parallel=True,chunks='auto',
                  verbose=False,cs=False,
                  **kwargs):
        """
        Specifies dataset to be sampled at obs location.
        On input,

        times, lons, lats: trajectory coordinates
        dataset: the input dataset, it can be one of these
                 xr.Dataset: an xarray dataset
                 string : either a GrASDS-style control file
                          or a glob template (e.g., *.nc)
                 list,tuple: a list of file names
        parallel: bool, whether to open dataset in parallel and return
                  dask arrays.
        verbose: bool, what it says.
        cs: sample on cubed sphere
        """

        self.verb = verbose
        self.times = xr.DataArray(times,dims='time',coords={'time': times})
        time_range = times.min(), times.max()
        if isinstance(time_range[0],np.datetime64):
            time_range = pd.to_datetime(time_range)

        # If dataset is an xarray dataset we are good to go
        # -------------------------------------------------
        if isinstance(dataset,xr.Dataset):
            self.ds = dataset # we are good to go...
        #  if dataset is a parquet referece file
        elif dataset[-4:] == 'parq':
            fs = fsspec.filesystem("reference", fo=dataset, remote_protocol='file', lazy=True)
            self.ds = xr.open_dataset(fs.get_mapper(), engine="zarr", chunks=chunks, backend_kwargs={"consolidated": False})
        # If dataset is a list of files...
        # OR GrADS-style ctl 
        # OR a glob type of template
        # --------------------------------
        elif isinstance(dataset,(list,tuple,str)):
            self.ds = xc.open_mfdataset(dataset,time_range=time_range,parallel=parallel,
                                        chunks=chunks,cs=cs,**kwargs) # special handles GrADS-style ctl if found

        else:
            raise SamplerError("Invalid dataset specification.")


        # Save coordinates with proper attributes
        # ---------------------------------------
        self.lons = xr.DataArray(lons, dims='time',attrs=self.ds.coords['lon'].attrs)
        self.lats = xr.DataArray(lats, dims='time',attrs=self.ds.coords['lat'].attrs)

        if cs:
            # Build ESMF unstructured mesh
            self.src_mesh = build_dual_mesh(self.ds)

            nobs = len(lons)
            # Create locstream
            self.dst_locstream = ESMF.LocStream(nobs, coord_sys=ESMF.CoordSys.SPH_DEG)
            self.dst_locstream["ESMF:Lon"] = np.asarray(lons, dtype=np.float64)
            self.dst_locstream["ESMF:Lat"] = np.asarray(lats, dtype=np.float64)
            
            self.src_field = ESMF.Field(self.src_mesh, name='src', meshloc=ESMF.MeshLoc.NODE)
            self.dst_field = ESMF.Field(self.dst_locstream, name='dst')
            
            self.regridder = ESMF.Regrid(
                self.src_field, 
                self.dst_field,
                regrid_method=ESMF.RegridMethod.BILINEAR,
                unmapped_action=ESMF.UnmappedAction.IGNORE
            )

        self.cs = cs
            
    #--
    def sample(self,Variables=None,method='linear'):
        """
        Sample variables and pre-determined obs locations

        """
        if Variables is None:
            if self.cs:
                # Only include variables with standard horizontal dims
                Variables = [
                    vn for vn in self.ds.data_vars
                    if all(d in self.ds[vn].dims for d in ('nf', 'Ydim', 'Xdim'))
                ]
            else:
                Variables = list(self.ds.data_vars)


        elif isinstance(Variables,str):
            Variables = [Variables,]

        sampled = dict()

        for vn in Variables:
            if self.verb: print('[ ] sampling',vn)
            if self.cs:
                da_var = self.ds[vn]
                
                # Identify extra dims (like time, lev)
                spatial_dims = ('nf', 'Ydim', 'Xdim')
                extra_dims = [d for d in da_var.dims if d not in spatial_dims]
                
                # Reorder so extra dims are first
                da_var = da_var.transpose(*extra_dims, *spatial_dims)
                
                # Calculate shapes
                extra_shape = tuple(da_var.sizes[d] for d in extra_dims)
                num_slices = int(np.prod(extra_shape)) if extra_shape else 1
                num_elements = int(np.prod([da_var.sizes[d] for d in spatial_dims]))
                num_locations = len(self.lons)
                
                # Reshape data into 2D and pre-allocate results
                flat_data_2d = da_var.values.reshape(num_slices, num_elements)
                interpolated_results = np.empty((num_slices, num_locations), dtype=da_var.dtype)
                
                # Spatial Regridding Loop
                for i in range(num_slices):
                    self.src_field.data[...] = flat_data_2d[i, :]
                    self.regridder(self.src_field, self.dst_field)
                    interpolated_results[i, :] = self.dst_field.data
                    
                # Reshape spatial results back to multidimensional
                final_shape = extra_shape + (num_locations,)
                interpolated_results = interpolated_results.reshape(final_shape)
                
                # Package into xarray DataArray for temporal interpolation
                out_coords = {d: da_var.coords[d] for d in extra_dims if d in da_var.coords}
                out_coords['locations'] = np.arange(num_locations) # temporary coordinate
                
                ds_spatial = xr.DataArray(
                    interpolated_results,
                    dims=extra_dims + ['locations'],
                    coords=out_coords,
                    attrs=da_var.attrs
                )
                
                # Temporal interpolation (if 'time' is a dimension)
                if 'time' in extra_dims:
                    times_da = xr.DataArray(self.times.values, dims='locations')
                    da_traj = ds_spatial.interp(time=times_da, method='linear')
                else:
                    da_traj = ds_spatial
                
                # Rename locations -> time (matching the trajectory format)
                sampled[vn] = da_traj.assign_coords(time=('locations', self.times.values)).swap_dims({'locations': 'time'}).drop_vars('locations', errors='ignore')

            else:
                sampled[vn] = self.ds[vn].interp(time=self.times,lon=self.lons,lat=self.lats,method=method)
                sampled[vn] = sampled[vn].assign_coords({'time': self.times})

        return xr.Dataset(sampled)

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

def CLI_stnSampler():

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

    parser.add_option("--cs",
                      action="store_true", dest="cs",
                      help="Sampling on cubed sphere")

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
    stn = STATION(stations,lons,lats,dataset,verbose=options.verbose,cs=options.cs)
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

#....................................................................................
def addVertCoord(aer):
        """
        adds a GEOS vertical coordinate derived from DELP and AIRDENS
        to an existing Dataset
        Useful for plotting and comparing to observations

        use it this way:
        aer.pipe(addVertCoord)

        Z = mid-level altitude above surface in km
        DZ = level thickness in km
        """
        from .constants import MAPL_GRAV as GRAV

        # GEOS files can be inconsistent when it comes to case
        # ----------------------------------------------------
        try:
            dp = aer['DELP']
        except:
            dp = aer['delp']

        # Get layer thicnkness from DELP & AIRDENS
        # DELP: Pa = kg m-1 s-2
        # AIRDENS: kg m-3
        # GRAV: m s-2
        # -----------------------------------------
        rhodz = dp / GRAV
        dz = rhodz / aer['AIRDENS']       # column thickness in m

        # add up the thicknesses to get edge level altitudes
        nlev = dz.sizes['lev']
        ze = xr.concat([dz.isel(lev=slice(i,None)).sum(dim='lev',keep_attrs=True) for i in range(nlev)],dim='lev')
        # append surface level, altitude = 0
        surface = xr.DataArray(np.zeros((1,) +ze.shape[1:]),dims=ze.dims)
        ze = xr.concat([ze,surface],dim='lev')
        # transpose back to npts,nlev again
        dims = list(dz.dims)
        ze = ze.transpose(*dims)
        # convert from m to km
        ze = ze*1e-3

        # get mid-level altitudes
        z = (ze.isel(lev=slice(None,-1)) + ze.isel(lev=slice(1,None)))*0.5

        # Attributes
        # ----------
        A = dict (Z = {'long_name':'mid-level atltitude above surface', 'units':'km'},
                  DZ = {'long_name':'level altitude thickness', 'units':'km'}
                  )

        # Pack results into a DataArray
        # ---------------------------
        DA = dict(  Z = z.assign_attrs(A['Z']),
                    DZ = dz.assign_attrs(A['DZ'])
                 )

        # Add to Dataset
        # ---------------
        for var in DA:
            aer[var] = DA[var]


#................................................................................
def CLI_trjSampler():

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
                          version='3.0.0' )

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

    parser.add_option("--cs",
                      action="store_true", dest="cs",
                      help="Sampling on cubed sphere")

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

    if options.Vars is not None:
        options.Vars = options.Vars.split(',')

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
            trj = TRAJECTORY(time,lon,lat,dataset,verbose=options.verbose,cs=options.cs)
            ds = trj.sample(Variables=options.Vars,method=method)

            if options.verbose:
                print('- Writing',outFile,'from',trjFile,'at takeoff',takeoff)

            ds.to_netcdf(outFile_,format=options.format,compute=True)

    # All else
    # --------
    else:

        trj = TRAJECTORY(time,lon,lat,dataset,verbose=options.verbose,cs=options.cs)
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



