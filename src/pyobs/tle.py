"""
Class for parsing TLE and computing satellite subpoint given time.
"""

import os
import ephem as ep
import numpy as np

from datetime import datetime, timedelta

class TLE(object):

    def __init__(self,tle_filename):
        """
        Loads TLE file.
        """
        name = os.path.basename(tle_filename).split('.')[-1]
        tle = dict(name=name)
        for line in open(tle_filename).readlines():
            if line[0] == '1': tle[1] = line.replace('\n','')
            if line[0] == '2': tle[2] = line.replace('\n','')

        self.tle = tle
        self.ephem = ep.readtle(tle['name'],tle[1],tle[2])

    def getSubpoint(self,t1,t2,dt):
        """
        Returns 2-tuple with arrays of lons, lats (in degrees)
        for a time interval [t1,t2] with time step dt.

        t1, t2: datetime, time interval
        dt    : timedelta, timestep

        """
        sublon, sublat = [], []
        t = t1
        while t <= t2:
            self.ephem.compute(ep.Date(t.isoformat(sep=' ').replace('-','/')))
            sublon.append(np.rad2deg(self.ephem.sublong))
            sublat.append(np.rad2deg(self.ephem.sublat))
            t += dt
                       
        return (np.array(sublon), np.array(sublat))

#.......................................................................

if __name__ == "__main__":

    tle_fn = '/Users/adasilva/data/tle/terra_2008.tle'

    tle = TLE(tle_fn)

    t1 = datetime(2008,1,11,0,45,0)
    t2 = datetime(2008,1,11,6,30,0)
    dt = timedelta(minutes=1)

    sublon, sublat = tle.getSubpoint(t1,t2,dt)

    
                
