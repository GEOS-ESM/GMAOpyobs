#!/usr/bin/env python
#
# Creates a trj_sampler resource file for a given bounding box.
# Some named cross-sections are defined internally for consistency.
# Aircraft speed hardwired for now, but should eventually be turned into
# an option.
#
# Arlindo da Silva, May 2016.

import os
import sys

from datetime import datetime, timedelta
from numpy import linspace, pi, zeros, sin, cos, arccos, array
from dateutil.parser import parse as isoparser

a = 6376.e3 # earth's radius

Sections = dict(     jeju2seoul = [127.1,32.,127.1,38.],
                   koreanstrait = [126.5,32.5,132.,37.],
                    seoul2busan = [126.8,37.8,130.1,34.],
                        westsea = [124.5,32.,124.5,39.],
                 westeast2seoul = [124.,37.5,132.,37.5],
                )


def rdist(lon,lat,Chordal=False):
    """
    RDIST	Chordal/Great Circle distance.

    R = RDIST(LON,LAT) returns the pairwise chordal or great circle
        distance matrix for the earth surface locations defined by the
        longitude-latitude coordinate arrays LON,LAT

                   R(i,j) = dist(u(i),u(j)) 

        where u(i) has lon-lat coordinates LON(i),LAT(i), u(j) has
        lon-lat coordinates LON(j),LAT(j), and dist(.,.)  is chordal
        or great circle distance depending on the parameter *Chordal*.
		
       The lon-lat coordinates are assumed to be given in degrees; distances
       are returned in meters.

       """

    n = lon.size

#   Cartesian coords on unity sphere 
#   --------------------------------
    slon = sin((pi/180)*lon[:])
    clon = cos((pi/180)*lon[:])
    slat = sin((pi/180)*lat[:])
    clat = cos((pi/180)*lat[:])
    x, y, z = (clat*clon, clat*slon, slat) 

#   Compute distances on unit sphere
#   --------------------------------
    R = zeros((n,n))
    if Chordal:
        for i in range(n):
            xi, yi, zi = ( x[i], y[i], z[i] )
            r2 = (x-xi)**2 + (y-yi)**2 + (z-zi)**2
            r2[r2<0.0] = 0.0
            R[i] = sqrt(r2)
    else:
        for i in range(n):
            xi, yi, zi = ( x[i], y[i], z[i] )
            c = x*xi + y*yi + z*zi
            c[c>1.]  = 1.
            c[c<-1.] = -1.
            R[i] = arccos(c)

#   Multiply by Earth's radius
#   --------------------------
    R = a * R

    return R

#------------------------------------------------------------------------
if __name__ == "__main__":

    # Bounding box on input
    # ---------------------
    if len(sys.argv) < 4:
        print "Usage:     %s  out_filename iso_t0 lon1 lat1 lon2 lat2"
        print "           %s  out_filename iso_t0 section_name"
        print "Sections: ", Sections.keys()
        sys.exit(1)

    if len(sys.argv) == 7:
        filename, iso_t0, bbox_ = sys.argv[1], sys.argv[2], sys.argv[3:7]
        bbox = [ float(x) for x in bbox_ ] 
    elif len(sys.argv) == 4:
        filename, iso_t0, section = sys.argv[1:]
        bbox = Sections[section]
    else:
        raise ValueError, "not enough parameters on input"

    t0 = isoparser(iso_t0) # take-off time
    
    # DC8 speed --- make this an option at some point
    # -----------------------------------------------
    V_knots = (425. + 490.) / 2. # http://www.nasa.gov/centers/dryden/research/AirSci/DC-8/dc8_his.html
    V_mps = 0.514444 * V_knots # meter/sec
    V_mpm = 60.* V_mps         # meter/min
    dt_min = timedelta(seconds=60) # 1 minute
    
    # Length of segment
    # -----------------
    Lons = linspace(float(bbox[0]),float(bbox[2]),num=2, endpoint=True)
    Lats = linspace(float(bbox[1]),float(bbox[3]),num=2, endpoint=True)
    R = rdist(Lons,Lats) # 2x2 matrix of distances (re-used here)
    d = R[0,1]
    N = int(0.5 + d / V_mpm)

    # Lons, lats
    # ----------
    Lons = linspace(float(bbox[0]),float(bbox[2]),num=N, endpoint=True)
    Lats = linspace(float(bbox[1]),float(bbox[3]),num=N, endpoint=True)
    tyme = array([ t0 + i * dt_min for i in range(N) ])
    
    # Write out stn rc file
    # ---------------------
    f = open(filename,'w')
    f.write('lon,lat,time\n')
    for i in range(N):
        f.write('%f,%f,%s\n'%(Lons[i],Lats[i],tyme[i].isoformat()))
    f.close()

    
