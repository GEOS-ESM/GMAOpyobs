"""
Class for parsing TLE and computing satellite subpoint given time.
"""

import os
import ephem as ep
import numpy as np

from datetime import datetime, timedelta

class TLE(object):

    """
    Uses PyEphem to compute satellite sublon/sublat from a TLE file.
    This code is not vectorized but seems good enough for government
    business.
    """
    
    def __init__(self,tleFilename):
        """
        Loads TLE file. Assumes 1 TLE per file.
        """
        name = os.path.basename(tleFilename).split('.')[-1]
        tle = dict(name=name)
        for line in open(tleFilename).readlines():
            if line[0] == '1': tle[1] = line.replace('\n','')
            if line[0] == '2': tle[2] = line.replace('\n','')

        self.tle = tle
        self.ephem = ep.readtle(tle['name'],tle[1],tle[2])

    def getSubpoint(self,t1,t2,dt):
        """
        Returns 3-tuple with arrays of times, lons, lats (in degrees)
        for a time interval [t1,t2] with timestep dt.

        t1, t2: datetime, time interval
        dt    : timedelta, timestep

        """
        times, sublon, sublat = [], [], []
        t = t1
        while t <= t2:
            self.ephem.compute(t)
            times.append(t)
            sublon.append(np.rad2deg(self.ephem.sublong))
            sublat.append(np.rad2deg(self.ephem.sublat))
            t += dt
                       
        return (np.array(times), np.array(sublon), np.array(sublat))

#.......................................................................

if __name__ == "__main__":

    tle_fn = '/Users/adasilva/data/tle/terra_2008.tle'

    tle = TLE(tle_fn)

    t1 = datetime(2008,1,11,0,45,0)
    t2 = datetime(2008,1,11,6,30,0)
    dt = timedelta(minutes=1)

    times, sublon, sublat = tle.getSubpoint(t1,t2,dt)

    
                
