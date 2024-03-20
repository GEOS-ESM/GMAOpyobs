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

try:
    from .aircrafts import platform

except:
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

def _greatCircle(startlong, endlong, startlat, endlat,nsegs):
    """
    Generate nsegs line segments between bounding coordinates using the great circle distance.
    """

    # calculate distance between points
    g = pj.Geod(ellps='WGS84')
    (az12, az21, dist) = g.inv(startlong, startlat, endlong, endlat)

    # calculate line string along path with npts segments
    lonlats = g.npts(startlong, startlat, endlong, endlat,nsegs-1)

    # npts doesn't include start/end points, so prepend/append them
    lonlats.insert(0, (startlong, startlat))
    lonlats.append((endlong, endlat))

    N = len(lonlats)
    lon, lat = np.zeros(N), np.zeros(N)
    for n in range(N):
        lon[n], lat[n] = lonlats[n]

    return (lon,lat)
    
#..................................................................................................

class WAYPOINT(object):

    def __init__ (self, wpFile, plane, verbose=False, refine=1 ):
        """
        Loads CSV waypoint file. On input,

        plane: str, aircraft name as defined in module aircrafts.
               If plane='snapshot', a stationary trajectory at
               the takeoff time will be generated.

        refine, int, number of segments in between each waypoint
        
        """

        self.verbose = verbose

        # Parse CSV file
        # --------------
        f = open(wpFile,"r")
        f.readline()
        self.city, self.utcOffset = f.readline().replace('\n','').split(',')
        self.wp = pd.read_csv(f)
        self.N = self.wp.shape[0]
        self.plane = plane

        # Use great circle distance to refine waypoints
        # ---------------------------------------------
        if refine > 1:
            lon, lat = [], []
            for n in range(self.N-1):
                lon_, lat_ = _greatCircle(self.wp.lon[n],self.wp.lon[n+1],
                                          self.wp.lat[n],self.wp.lat[n+1],refine)
                lon.append(lon_[:-1])
                lat.append(lat_[:-1])
            lon.append(np.array([lon_[-1],]))
            lat.append(np.array([lat_[-1],]))
            lon = np.concatenate(lon)
            lat = np.concatenate(lat)
            #breakpoint()
            self.wp = pd.DataFrame({'lon':lon, 'lat':lat})
            self.N = self.wp.shape[0]
            
        if self.verbose:
            print('- '+self.city+' is %s hours later than UTC'%self.utcOffset)

    #--        
    def getTraj(self,takeoff):
        """
        Calculates trajectory for a given takeoff local time. On input,

        takeoff: str or time delta, local takeoff datetime in ISO format.

        refine: float, factor for refining waypoints. Refine=10 will refine the waypoints
                by adding 10 subintervals

        Returns DataFrame with trajecotry coordinates.
        
        """

        # Generate time coordinates
        # -------------------------
        if isinstance(takeoff,str):
            takeoff_ = isoparser(takeoff)
        elif isinstance(takeoff,datetime):
            takeoff_ = takeoff
        else:
            raise ValueError("takeoff must be str or datetime.")           
        t0 = takeoff_ - timedelta(hours=int(self.utcOffset)) # UTC
        time = np.repeat(t0,self.N)

        if self.plane != 'snapshot':

            geod = pj.Geod(ellps='WGS84')
            _, _, dist = geod.inv(self.wp.lon[0:-1],self.wp.lat[0:-1],self.wp.lon[1:],self.wp.lat[1:])
            speed = platform[self.plane]['mean_speed'] # m/s
            dt = dist / speed
            time[1:] += np.array([timedelta(seconds=s) for s in dt.cumsum()])


        # Create trajectory DataFrame
        # ---------------------------
        traj = pd.DataFrame(dict(lon=self.wp.lon.values, lat=self.wp.lat.values), index=time)

        return traj
    
    #---
    def writeTraj(self,takeoff,outFile='@city_@aircraft_@takeoff.csv',format=None):
        """
         Writes trajectory at local takeoff time. On input:

         outFile: str, output file name. Default is @city_@aircraft_@takeoff.csv, where
                  @city, @aircraft and @takeoff are replaced with actual values.
         format: str, either "csv" or "netcdf". Default: derived from outFile extension.
         
        """

        # Compute trajectory
        # ------------------
        traj = self.getTraj(takeoff)

        # Create consistent file name extension
        # -------------------------------------
        name, ext = os.path.splitext(outFile)
        if ext.lower() == '.csv':
            format = 'csv'
        elif ext.lower() == '.nc':
            format = 'netcdf'
        elif ext.lower() == '.gz':
            format = 'gzip'

        if format == 'csv':
            outFile = name + '.csv'
        elif format == 'netcdf':
            outFile = name + '.nc'
        elif format == 'gzip':
            outFile = name + '.gz'
        else:
            raise ValueError('invalid extension <%s>'%ext)

        # Default file name
        # -----------------
        if '@city_@aircraft_@takeoff' in outFile:
            outFile = outFile.replace('@city',self.city).\
                              replace('@aircraft',self.plane).\
                              replace('@takeoff',str(takeoff).replace(' ','T'))

        # Write out results
        # -----------------
        if self.verbose:
            print('- Writing',outFile)

        if format == 'csv' or format == 'gzip':

            traj.index = traj.index.map(lambda x: datetime.strftime(x, '%Y-%m-%dT%H:%M:%SZ'))
            traj.to_csv(outFile,index_label='time')

        elif self.format == 'netcdf' or self.format == 'nc':

            traj.to_xarray().rename({'index':'time'}).to_netcdf(outFile)


#..................................................................................................
def CLI_wp2traj():

    plane = 'DC8'
    outFile = '@city_@aircraft_@takeoff.csv'
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

    parser.add_option("-r", "--refine", dest="refine", default=1,
              help="Refine the waypoints with REFINE segments based on great circle distance (default=1)" )


    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")

    (options, args) = parser.parse_args()
    
    if len(args) >= 2 :
        wpFile = args[0]
        TakeOff = args[1:]
    else:
        parser.error("must have 2 arguments: wpFile takeOff")


    # Instantiate waypoint
    # --------------------
    wp = WAYPOINT(wpFile, options.plane, refine=int(options.refine),verbose=options.verbose)

    # Write out files
    # ---------------
    for takeoff in TakeOff:
        wp.writeTraj(takeoff,options.outFile,options.format)

if __name__ == "__main__":
    CLI_wp2traj()
