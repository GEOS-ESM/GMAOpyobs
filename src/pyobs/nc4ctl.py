"""
    Extends the netCDF4 class with agreegation and some CF capabilities.
 
"""

import os
from time     import time
from numpy    import linspace, ones, zeros, any, array, float32, tile
import numpy as np
from datetime import datetime, timedelta
from dateutil.parser import parse         as isoparser

from netCDF4 import Dataset as Dataset_

#from binObs_ import interpxy3d
from dateutil.relativedelta import relativedelta


def interpxy3d(gField,glon,glat,lon,lat,periodic):


    nobs       = len(lon)
    im,jm,km   = gField.shape
    oField     = np.ones([km,nobs])*-999

    if periodic == 1:
        m_dlon = float(im) / 360.     # assume global & periodic
        m_dlat = float(jm-1) / 180.
    else:
        m_dlon = float(im-1) / ( glon[im-1] - glon[0] )
        m_dlat = float(jm-1) / ( glat[jm-1] - glat[0] )

    # Loop over observations
    # ------------------------------
    for nob in np.arange(nobs):
        # Longitude
        # --------------
        o_lon = 1. + (lon[nob]-glon[0]) * m_dlon
        if ( periodic == 1 ):
            i   = min(im, int( o_lon ))
            alfa  = o_lon - i
            if (i == im):
                i1 = im-1
                i2 = 0
            else:
                i1 = i-1
                i2 = i
        else:
            i   = min(im-1, int( o_lon ))
            alfa  = o_lon - i
            i1 = i-1
            i2 = i

        # Latitude
        # -------------
        o_lat = 1. + (lat[nob] - glat[0]) * m_dlat
        j   = min( jm-1, int( o_lat ) )
        beta  = o_lat - j
        a11 = gField[i1,j-1,  :]   #!W-S
        a21 = gField[i2,j-1,  :]   #!E-S
        a12 = gField[i1,j,:]       #!W-N
        a22 = gField[i2,j,:]       #!E-N
        a00 = a11 + alfa * ( a21 - a11 )

        oField[:,nob] = a00 + beta * ( a12 + alfa * ( a22 - a12 ) - a00)

    return oField

class NC4ctlHandle(object):
    """
    A simple container class to collect NC4ctl arrays.
    """
    def __init__ (self, name):
        self.name = name
    
class NC4ctlError(Exception):
    """
    Defines NC4ctl general exception errors.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Dataset(Dataset_):

#---
    def coordVars(self):
        """
        Identify coordinate variables and type.  Avoid dependence on udunits for now.
        """
        coords = dict(lon=None,lat=None,time=None,lev=None, tunits=None, toff = None)
        Vars = list(self.variables.keys())
        for d in list(self.dimensions.keys()): 
            if d in Vars:          # potential coordinate variable
               try:
                   units = self.variables[d].units.lower()
                   if   'degrees_north' in units: coords['lat'] = d
                   elif 'degrees_south' in units: coords['lat'] = d
                   elif 'degrees_east'  in units: coords['lon'] = d
                   elif 'degrees_west'  in units: coords['lon'] = d
                   elif 'pa'            in units: coords['lev'] = d
                   elif 'mb'            in units: coords['lev'] = d
                   elif 'meter'         in units: coords['lev'] = d
                   elif 'level'         in units: coords['lev'] = d
                   elif 'sigma_level'   in units: coords['lev'] = d
                   elif 'seconds'       in units: coords['time'] = d
                   elif 'minutes'       in units: coords['time'] = d
                   elif 'hours'         in units: coords['time'] = d
                   elif 'days'          in units: coords['time'] = d
               except:
                   pass

        if coords['time'] != None:
            try:
                isoT = self.variables['time'].units.split('since ')[1]
                tunt = self.variables['time'].units.split(' since')[0].lower()
                if ':' not in isoT: isoT = isoT+':0:0'
                coords['toff'] = isoparser(isoT)
                if tunt == 'seconds': tunits = 1.
                elif tunt == 'minutes': tunits = 60.
                elif tunt == 'hours': tunits = 60.**2
                elif tunt == 'days': tunits = 24 * 60.**2
                coords['tunits'] = tunits
            except:
                self.variables['time'] = None
               
        return coords
#.........................................................................................................

class NC4ctl(object):
    """
    Generic class to implement an aggregation functionality similar to
    templates in a GrADS control file. For simplicity, this class is
    only implemented for reading files, not writing them.
    """

    def __init__ (self, ctlfile ):
        """
        Initialize an aggregated NC4ctl object.
        """
        # Parse CTL file
        # --------------
        CTL = open(ctlfile).readlines()
        dset, template, nt = (None,False, None)
        for line in CTL:
            tokens = line.replace('\r','').replace('\n','').split()
            keyw = tokens[0].upper()
            if keyw== 'DSET':
                dset = tokens[1]
            elif keyw == 'OPTIONS':
                if 'TEMPLATE' in line.upper():
                    template = True
            if keyw == 'TDEF':
                if len(tokens) == 5:
                    tdef, nt, linear, t0, dt = tokens
                elif len(tokens) == 6:
                    tdef, dim, nt, linear, t0, dt = tokens
                else:
                    raise ValueError('Invalid TDEF record: '+line)

        # Consistency check
        # -----------------
        if dset is None or nt is None:
            raise ValueError('<%s> does not seem to be a valid GrADS control file'%ctlfile)
        else:
            if '^' in dset:
                dirn = os.path.dirname(ctlfile)
                dset = dset.replace('^',dirn+'/')
        if template is False:
            raise ValueError('<%s> does not seem to be templated'%ctlfile)

        # Handle time attributes
        # ----------------------
        dt = dt.lower()
        if 'hr' in dt:
            secs = int(dt.replace('hr','')) * 60 * 60
        elif 'mn' in dt:
            secs = int(dt.replace('mn','')) * 60
        elif 'dy' in dt:
            secs = int(dt.replace('dy','')) * 24 * 60 * 60
        elif 'mo' in dt:
            mons = int(dt.replace('mo',''))            
        else:
            raise ValueError('invalid time step <%s>'%dt) 

        if 'mo' in dt:
            dt = relativedelta(months=+mons)
        else:
            dt = timedelta(seconds=secs)

        # Save this
        # ---------
        self.dset = dset
        self.lm = int(nt)
        self.dt = dt
        self.tbeg = _gat2dt(t0)
        self.tend = self.tbeg + (self.lm-1) * self.dt
        self.Files = dict()

        # Open first file in sequence
        # ---------------------------
        self.openTime(self.tbeg)

#---
    def openTime(self,tyme):
        """
        Open time for this time.
        """
        filename = _strTemplate(self.dset,tyme=tyme)
        try:
            self.nc = self.Files[filename]
            self.coords = self.nc.coordVars()
        except:
            self.Files[filename] = Dataset(filename)
            self.nc = self.Files[filename]
            self.coords = self.nc.coordVars()
        
#---
    def tbracket (self, t):
        """
        Given (t1,t2) find bracketing times on file.
        """
        if t<self.tbeg:
            raise ValueError('%s before %s'%(str(t),str(self.tbeg)))
        elif t>self.tend:
            raise ValueError('%s after %s'%(str(t),str(self.tend)))

        # dt = t - self.tbeg
        # i = int(dt.total_seconds() / self.dt.total_seconds())
        # t1 = self.tbeg + i * self.dt
        # t2 = t1 + self.dt

        t1 = self.tbeg
        while t1 <= t:
          t1 = t1 + self.dt

        t1 = t1 - self.dt
        t2 = t1 + self.dt

        return (t1,t2)
    
#---
    def trange(self, t1, t2):
        """
        Return file times between t1 & t2.
        """
        ta, tb = self.tbracket(t1)
        tc, td = self.tbracket(t2)
        if tc==t2:
            td = tc
        T = []
        t = ta
        while t<=td:
            T.append(t)
            t += self.dt

        return T
                
    #---
    def interpXY_LatLon(self,vname,lon,lat,tyme,
                        squeeze=True,Transpose=False,
                        kbeg = 1, kount=None,
                        algorithm='linear',**kwds):

        """
        Interpolate grided Lat/Lon fields. Not as general as it should be,
        therefore, not made the default interpXY method.

        On input, kbeg assumes fortran indexing, starting at 1.
        """

        if tyme is not None:
            self.openTime(tyme)

        V = self.nc.variables[vname]
        Lon = self.nc.variables[self.coords['lon']][:]
        Lat = self.nc.variables[self.coords['lat']][:]

        # Tighten domain to minimize I/O
        # ------------------------------
        lon_min, lon_max = lon.min(), lon.max()
        lat_min, lat_max = lat.min(), lat.max()

        if (lon_min<=Lon.min()) or (lon_max>=Lon.max()):
            periodic = 1
        else:
            periodic = 0

        dLon, dLat = Lon[1]-Lon[0], Lat[1]-Lat[0] # assume constant

        I = (Lon>=(lon_min-dLon)) & (Lon<=(lon_max+dLon))
        J = (Lat>=(lat_min-dLat)) & (Lat<=(lat_max+dLat))

        X, Y = Lon[I], Lat[J]

        rank = len(V.shape)
        n = len(lon)
        if rank == 3:
            t, y, x = V.shape
            if t>1:
                raise ValueError('InterpXY_LatLon assumes 1 time per file for now')
            Z = V[0:1,J,I] # time here standing in for Z 
            v = interpxy3d(Z.T,X,Y,lon,lat,periodic)
            v = v[0,:]
        elif rank == 4:
            t, z, y, x = V.shape
            if kount is None:
                k1, k2 = kbeg-1, z # fortran indexing on input
                kount = z
            else:
                k1, k2 = kbeg-1, kbeg-1+kount # fortran indexing on input
            if k1<0 or k2>z:
                raise ValueError('invalid kbeg=%d, kount=%d'%(kbeg,kount))
            if t>1:
                raise ValueError('InterpXY_LatLon assumes 1 time per file for now, ask a guru to generalize me.')
            Z = V[0,k1:k2,J,I] 
            v = interpxy3d(Z.T,X,Y,lon,lat,periodic) # cache optimized
        else:
            raise ValueError('invalid rank = %d'%rank)

        if Transpose: 
            v = v.T

        if squeeze:
            v = v.squeeze()

        return v
            
#---
    def interpXY_PlaceHolder(self,vname,lon,lat,tyme=None,
                             squeeze=True,Transpose=False,
                             algorithm='linear',**kwds):
        """
        Reads a variable at a given time/date, or the first time on
        file if nymd/nhms is not specified, and interpolates it to the
        observation locations given by the list of (lon,lat) on
        input. By default it returns all vertical levels of the
        variable, but a range can be specified with the first vertical
        index (kbeg) and the number of vertical levels to be read
        (kount). Keyword arguments **kwds are passed to the NC4ctl
        interp() method.

        Example:
    
        v = self.interpXY('delp',lon,lat)

        By default, when the variable being interpolated is 2-D, the
        output array will have shape (nobs,); when the variable being
        interpolated is 3-D, the output array will have shape
        (km,nobs), where nobs = len(lon), and km is the number of
        vertical levels being requested. If Transpose=True, the output
        array is transposed.
        
        """

        raise ValueError("user must supply own interp() method for now.")


#---

    interpXY = interpXY_PlaceHolder # By default user must choose own

#---

    def sample(self,vname,lon,lat,time,
               squeeze=True,Transpose=False,
               algorithm='linear',Verbose=False,**kwopts):
        """
        Interpolates *vname* to (lon,lat,time) trajectory;
        time must be in ascending order. Keyword arguments **kwds are
        passed to the interpXY method().

        By default, output arrays have shape (km,nobs), where *km* is the
        requested number of verical levels. If Transpose=True is specified,
        then the output arrays will be transposed to (nobs,km).
    
        """

        # Inputs must be 1D arrays
        # ------------------------
        if len(lon.shape)!=1 or len(lat.shape)!=1 or len(time.shape)!=1:
            raise ValueError("lons, lats, time must be 1D arrays")
        
        # Find times bracketing the input time array
        # ------------------------------------------
        Times = self.trange(time.min(),time.max())
        # dt, dt_secs = (self.dt, self.dt.total_seconds())
        dt = self.dt

        # Loop over time, producing XY interpolation at each time
        # -------------------------------------------------------
        V, I = [], []
        for now in Times:
            if Verbose: print(" [] XY Interpolating <%s> at "%vname,now)
            i = (time>=now-dt) & (time<=now+dt)            
            if any(i):
                v = self.interpXY(vname, lon[i], lat[i],tyme=now,
                  Transpose=True, # shape will be (nobs,km)
                  squeeze=False,algorithm=algorithm,**kwopts)                
            else:
                v = None
            V.append(v)
            I.append(i)
            
        # Now perform the time interpolation
        # ----------------------------------
        if Verbose: print(" () T  Interpolating <%s>"%vname)
        N = len(lon)
        if len(v.shape)==2:
            km = v.shape[1]
        else:
            km = 0
        if km>0:
            shp = [N,km]
        else:
            shp = [N,]

        if len(Times) > 1:
            v  = zeros(shp,dtype=float32)
            v1, v2 = v.copy(), v.copy() # scratch space
            n = 0
            for now in Times[:-1]:
                v1[I[n]], v2[I[n+1]] = V[n], V[n+1]
                j = (time>=now) & (time<=now+dt)
                dt_secs = ((now+dt)-now).total_seconds()
                if any(j): 
                    a = array([r.total_seconds()/dt_secs for r in time[j]-now],dtype=float32) 
                    if len(shp)==2: # has vertical levels
                        a = tile(a,(shp[1],1)).T # replicate array
                    v[j] = (1-a) * v1[j] + a * v2[j]
                n += 1
        else:
            v = V[0]

        if Transpose == False: v = v.T # back to NC4ctl's (km,nobs)
        if squeeze == True:    v = v.squeeze()

        return v

#...........................................................................

__Months__ = ['JAN','FEB','MAR','APR','MAY','JUN',
              'JUL','AUG','SEP','OCT','NOV','DEC']

def _strTemplate(templ,expid=None,nymd=None,nhms=None,
                    yy=None,mm=None,dd=None,h=None,m=None,s=None,
                    tyme=None):
    """
    Expands GrADS template in string *templ*. On input,

       expid ---  experiment id, expands %s
       yy    ---  year, expands %y4 and %y2
       mm    ---  month, expands %m2 or %m3
       dd    ---  day, expands %d2
       h     ---  hour, expands %h2
       m     ---  minute, expands %n2
       s     ---  minute, expands %S2 (notice capital "S")

       nymd  ---  same as yy*10000 + mm*100 + dd
       nhms  ---  same as h *10000 + h*100  + s

       tyme ---  python datetime

    Unlike GrADS, notice that seconds are expanded using the %S2 token. 
    Input date/time can be either strings or integers.

    Examples:

    >>> templ = "%s.aer_f.eta.%m3%y2.%y4%m2%d2_%h2:%n2:%S2z.nc"
    >>> print strTemplate(templ,expid="e0054A",yy=2008,mm=6,dd=30,h=1,m=30,s=47)
    e0054A.aer_f.eta.jun08.20080630_01:30:47z.nc
    >>> print strTemplate(templ,expid="e0054A",nymd=20080630,nhms=13000)
    e0054A.aer_f.eta.jun08.20080630_01:30:00z.nc

    NOTE: This function exists in MAPL/config.py; it is copied here for
          dependency management.
          
    """

    MMM = ( 'jan', 'feb', 'mar', 'apr', 'may', 'jun', 
            'jul', 'aug', 'sep', 'oct', 'nov', 'dec' ) 
    
    str_ = templ[:]

    if tyme is not None:
        yy = tyme.year
        mm = tyme.month
        dd = tyme.day
        h  = tyme.hour
        m  = tyme.minute
        s  = tyme.second

    if nymd is not None:
        nymd = int(nymd)
        yy = nymd/10000
        mm = (nymd - yy*10000)/100
        dd = nymd - (10000*yy + 100*mm )

    if nhms is not None:
        nhms = int(nhms)
        h = nhms/10000
        m = (nhms - h * 10000)/100
        s = nhms - (10000*h + 100*m)

    if expid is not None: 
        str_ = str_.replace('%s',expid)
    if yy is not None: 
        y2 = yy%100
        str_ = str_.replace('%y4',str(yy))
        str_ = str_.replace('%y2',"%02d"%y2)
    if mm is not None: 
        mm = int(mm)
        mmm = MMM[mm-1]
        str_ = str_.replace('%m2',"%02d"%mm)
        str_ = str_.replace('%m3',mmm)
    if dd is not None: 
        str_ = str_.replace('%d2',"%02d"%int(dd))
    if h  is not None: 
        str_ = str_.replace('%h2',"%02d"%int(h))
    if m  is not None: 
        str_ = str_.replace('%n2',"%02d"%int(m))
    if s  is not None: 
        str_ = str_.replace('%S2',"%02d"%int(s))

    return str_

#...........................................................................

def _splitDate(nymd):
    """
    Split nymd into year, month, date tuple.
    """
    nymd = int(nymd)
    yy = nymd/10000
    mm = (nymd - yy*10000)/100
    dd = nymd - (10000*yy + 100*mm )

    return (yy,mm,dd)

def _nymd2dt(nymd,nhms=0):
    """
    Return python time given (nymd,nhms).
    """
    yy, mm, dd = _splitDate(nymd)
    h, m, s = _splitDate(nhms)
    return datetime(yy,mm,dd,h,m,s)
    
def _hms2td(nhms):
    """
    Return python time delta given nhms.
    """
    h, m, s = _splitDate(nhms)
    return timedelta(seconds=(h*3600+m*60+s))
    
def _gat2dt(gat):
    """
    Convert grads time to datetime.
    """
    try:
        time, date = gat.upper().split('Z')
    except:
        time = '0'
        date = gat.upper()
    if time.count(':') > 0:
        h, m = time.split(":")
    else:
        h = time
        m = '0'
    mmm = date[-7:-4]
    dd, yy = date.split(mmm)
    mm = __Months__.index(mmm) + 1
    dt = datetime(int(yy),int(mm),int(dd),int(h),int(m))
    return dt

#...........................................................................

class NC4ctl_(NC4ctl):
    interpXY = NC4ctl.interpXY_LatLon

if __name__ == "__main__":

    from time import time as now
    from gfio import GFIOctl

    kbeg = 36
    kount = 72 - kbeg + 1

    filename = '/home/adasilva/iesa/MERRAero/inst3d_prog_v.ddf'

    # interp using NC4ctl
    # -------------------
    print("------ NC4ctl ----------")
    f = NC4ctl_(filename)
    lons = linspace(-45.,-30.,100)
    lats = linspace(30.,50.,100)
    dt = f.dt/10
    t0 = datetime(2008,6,29,12)
    times = array([ t0 + i * dt for i in range(100)])
    slp = f.sample('SLP',lons,lats,times,Verbose=True)
    t = f.sample('T',lons,lats,times,kbeg=kbeg,kount=kount,Verbose=True)

    # Redo with GFIOctl for testing
    # -----------------------------
    print("\n------ GFIOctl ----------")
    g = GFIOctl(filename)
    SLP = g.sample('SLP',lons,lats,times,Verbose=True)
    T = g.sample('T',lons,lats,times,Verbose=True,kbeg=kbeg,kount=kount)
