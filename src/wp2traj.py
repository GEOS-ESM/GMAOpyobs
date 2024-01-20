#!/usr/bin/env python
"""

Convert waypoint files to (time,lon,lat) trajectory CSV.

"""

import os
import pandas as pd
import numpy  as np
import pyproj as pj

from datetime        import datetime, timedelta
from dateutil.parser import parse as isoparser
from optparse        import OptionParser

# When running under pyobs we will use this
# from pyobs.aircrafts import platform

# For now, hardwire DC8.
platform = dict(
DC8 = {'Platform':'dc8','names':['dc8','DC8','DC-8','dc-8','DC 8','dc 8','Dc','dC'],
    'max_alt':13000.0,'base_speed':130.0,'speed_per_alt':0.0075,
    'mean_speed': 136,
    'max_speed':175.0,'max_speed_alt':6000.0,'descent_speed_decrease':15.0,
    'climb_vert_speed':15.0,'descent_vert_speed':-10.0,'alt_for_variable_vert_speed':0.0,
    'vert_speed_base':15.0,'vert_speed_per_alt':0.001,
    'rate_of_turn':None,'turn_bank_angle':15.0,
    'warning':False},
)

def write_traj(options,takeoff):
    """
    Write trajectory file for a givem take off.
    """
    
    # Default file name
    # -----------------
    if 'city_aircraft_takeoff' in options.outFile:
        outFile = options.outFile.replace('city',city).replace('aircraft',options.plane).replace('takeoff',takeoff)
    else:
        outFile = options.outFile
        
    # Generate time coordinates
    # -------------------------
    t0 = isoparser(takeoff) - timedelta(hours=int(utcOffset)) # UTC
    time = np.repeat(t0,N)

    if options.plane != 'snapshot':
        
        geod = pj.Geod(ellps='WGS84')
        _, _, dist = geod.inv(wp.lon[0:-1],wp.lat[0:-1],wp.lon[1:],wp.lat[1:])
        speed = platform[options.plane]['mean_speed'] # m/s
        dt = dist / speed
        time[1:] += np.array([timedelta(seconds=s) for s in dt.cumsum()])
        

    # Create trajectory DataFrame
    # ---------------------------
    traj = pd.DataFrame(dict(lon=wp.lon.values, lat=wp.lat.values), index=time)
    
    #if options.verbose:
    #    print(traj)

    # Write out results
    # -----------------
    if options.verbose:
        print('- Writing',outFile)
        
    if options.format == 'csv' or options.format == 'gzip':

        traj.index = traj.index.map(lambda x: datetime.strftime(x, '%Y-%m-%dT%H:%M:%SZ'))
        traj.to_csv(outFile,index_label='time')

    elif options.format == 'netcdf' or options.format == 'nc':

        traj.to_xarray().rename({'index':'time'}).to_netcdf(outFile)


if __name__ == "__main__":

    plane = 'DC8'
    outFile = 'city_aircraft_takeoff.csv'
    format = 'csv'

    #   Parse command line options
    #   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] wpFile takeOff_isoLocalTime(s)",
                          version='1.0.0' )

    parser.add_option("-o", "--output", dest="outFile", default=outFile,
              help="Output file (default=%s)"\
                          %outFile )

    parser.add_option("-p", "--platform", dest="plane", default=plane,
              help="Platform (default=%s). Specify '-p snapshot' for a snapshot at takeoff time "%plane)

    parser.add_option("-f", "--format", dest="format", default=format,
              help="Output file format: one of 'csv' or 'netcdf' (default=%s)"%format )


    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")

    (options, args) = parser.parse_args()
    
    if len(args) >= 2 :
        wpFile = args[0]
        TakeOff = args[1:]
    else:
        parser.error("must have 2 arguments: wpFile takeOff")

    # Parse CSV file
    # --------------
    f = open(wpFile,"r")
    f.readline()
    city, utcOffset = f.readline().replace('\n','').split(',')
    wp = pd.read_csv(f)
    N = wp.shape[0]


    # Create consistent file name extension
    # -------------------------------------
    name, ext = os.path.splitext(options.outFile)
    if ext.lower() == '.csv':
        options.format = 'csv'
    elif ext.lower() == '.nc':
        options.format = 'netcdf'
    elif ext.lower() == '.gz':
        options.format = 'gzip'
        
    if options.format == 'csv':
        options.outFile = name + '.csv'
    elif options.format == 'netcdf':
        options.outFile = name + '.nc'
    elif options.format == 'gzip':
        options.outFile = name + '.gz'
    else:
        raise ValueError('invalid extension <%s>'%ext)

    for takeoff in TakeOff:
        write_traj(options,takeoff)
    
        
