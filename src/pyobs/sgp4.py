"""

  Interface to SGP4 orbital calculator.

"""

import sgp4_

from datetime import datetime, timedelta
from numpy import sin, cos, pi, arccos, dot, sqrt, array

def getTrack ( tleFile, t_beg, t_end, dt_secs=60):
    """
    Given a Two Line Element (TLE) file name, a time interval, returns
    tuple with (lon,lat) coordinates of satellite ground track.
    """
    dt = t_end-t_beg
    n = 1 + int(dt.total_seconds() / dt_secs)
    Dt = timedelta(seconds=dt_secs)
    
    nymd = [ _getNYMD(t_beg), _getNYMD(t_end) ]
    nhms = [ _getNHMS(t_beg), _getNHMS(t_end) ]

    tyme = array([ t_beg + i * Dt for i in range(n) ])

    lon, lat, rc = sgp4_.sgp4track(n, tleFile, nymd, nhms, dt_secs)

    if rc>1:
        raise ValueError('Error on return from sgp4Track: <%d>'%rc)

    return (lon,lat,tyme)

#--
def chDist (x1,y1,z1,x2,y2,z2):
    """
    Compute chordal distance given 2 set of points.
    """

    # Great-circle distance for each point
    d = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2

    # L2-norm
    return d.sum()/(3*d.size)

#--
def gcPeriod(lon,lat,n1,n2):
    """
    Try to find period in the range (n1, n2)
    """
    x, y, z = _getXYZ(lon,lat)
    N = list(range(n1,n2+1))
    D = []
    for n in N:
        i = list(range(0,n))
        j = list(range(n,n+n))
        d_ = gcDist(x[i],y[i],z[i],
                    x[j],y[j],z[j])
        D = D + [d_,]
        print(n, d_)

    return (N,D)

#--
def chPeriod(lon,lat,n1,n2):
    """
    Try to find period in the range (n1, n2)
    """
    x, y, z = _getXYZ(lon,lat)
    N = list(range(n1,n2+1))
    D = []
    for n in N:
        i = list(range(0,n))
        j = list(range(n,n+n))
        d_ = chDist(x[i],y[i],z[i],
                    x[j],y[j],z[j])
        D = D + [d_,]
        print(n, d_)

    return (N,D)

#--
def dayPeriod(lon,lat,n1,n2,day):
    """
    Try to find period in the range (n1, n2) days.
    """
    x, y, z = _getXYZ(lon,lat)
    N = list(range(n1,n2+1))
    D = []
    for n_ in N:
        n = n_ * day
        i = list(range(0,n))
        j = list(range(n,n+n))
        d_ = gcDist(x[i],y[i],z[i],
                    x[j],y[j],z[j])
        D = D + [d_,]
        print(n, d_)

    return (N,D)

#........................................................................................

def _getNYMD(t):
    return  10000*t.year + 100*t.month + t.day

def _getNHMS(t):
    return  10000*t.hour + 100*t.minute + t.second

def _getXYZ ( lon, lat ):
    """
    Returns cartesian coords on unit sphere.
    """
    d2r = pi / 180.
    rlon, rlat = ( d2r * lon, d2r * lat )
    x = cos(rlat) * cos(rlon)
    y = cos(rlat) * sin(rlon)
    z = sin(rlat)
    return (x,y,z)

#........................................................................................

if __name__ == "__main__":

    from datetime import datetime as T

    dt_secs = 30

#    lon, lat, tyme = getTrack ( 'TLE/quikscat_2005.tle', T(2000,1,1), T(2000,2,15), dt_secs=dt_secs )
    (lon, lat, tyme) = getTrack ( 'ISS-2016257.tle', T(2016,9,14,10,30), T(2016,9,14,11,30), dt_secs=dt_secs )

    for i in range(len(lon)):
        print("%f,%f,%02d:%02d"%(lon[i],lat[i],tyme[i].hour,tyme[i].minute))
     
def hold():

    tlocal = array([ tyme[i] + timedelta(seconds=int(240*lon[i])) for i in range(tyme.size) ])

    # Determine equatorial crossing times
    # -----------------------------------
    n = lat.size - 1
    z = lat[0:n] * lat[1:]
    I = z<=0
    tlocal = tlocal[0:n]
    teq = tlocal[I]
    meq = array([ t.minute for t in teq ])
    heq = array([ t.hour   for t in teq ]) + meq/60.
    
    day = 24 * 60 * 60 / dt_secs # minutes in a day

    n1, n2 = (int(15.8*day), int(16.2*day) )

    N, D = chPeriod(lon,lat,n1,n2)

    lon, lat, tyme = getTrack ( 'TLE/quikscat_2005.tle', T(2005,1,1), T(2005,2,15), dt_secs=60 )

    day = 24 * 60 # minutes in a day

    N, D = period(lon,lat,2,20)

    n1, n2 = (int(14 * day), int(17 * day))
    n1, n2 = (2*day, 5*day)

    N, D = period(lon,lat,n1,n2)

    lon1, lat1, t1 = getTrack ( 'TLE/quikscat_2005.tle', T(2005,1,1),  T(2005,1,14), dt_secs=60 )
    lon2, lat2, t2 = getTrack ( 'TLE/quikscat_2005.tle', T(2005,1,15), T(2005,1,28), dt_secs=60 )

