#!/usr/bin/env python
#
# Creates a stn_sampler resource file for a given bounding box.
#

import os
import sys

from numpy import linspace

if __name__ == "__main__":

    # Bounding box on input
    if len(sys.argv) < 6:
        print "Usage:   %s  filename lon1 lat1 lon2 lat2 [N]"
        sys.exit(1)

    if len(sys.argv) == 7:
        filename, bbox, N = sys.argv[1], sys.argv[2:6], int(sys.argv[6])
    else:
        filename, bbox, N = sys.argv[1], sys.argv[2:6], 128

    # Lons, lats
    # ----------
    Lons = linspace(float(bbox[0]),float(bbox[2]),num=N, endpoint=True)
    Lats = linspace(float(bbox[1]),float(bbox[3]),num=N, endpoint=True)

    # Write out stn rc file
    # ---------------------
    f = open(filename,'w')
    f.write('name,lon,lat\n')
    for i in range(N):
        f.write('p%04d,%f,%f\n'%(i+1,Lons[i],Lats[i]))
    f.close()

    
