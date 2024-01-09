"""
   Extends Xarray's open_mfdaset by recognizing GrADS-style control (ctl) files.

   Note: CHSUB not supported. 

"""

import os

import xarray as xr
import numpy  as np

from datetime import datetime, timedelta
from dateutil.parser import parse         as isoparser
from dateutil.relativedelta import relativedelta
    
class XRctlError(Exception):
    """
    Defines NC4ctl general exception errors.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

#...........................................................................

def open_mfdataset(paths,*args, time_range=None,**kwargs):
    """
    Intercepts call to xarray open_mfdataset() and if *paths*
    is a GrADS-style ctl file, parses it generating a list of
    files that are then passed down to xr.open_mfdataset(). 
    """
    paths_ = paths
    if isinstance(paths,str):
        if paths_.split('.')[-1] in ('ctl','xdf', 'ddf'):  # GrADS style control file
                paths_ = parse_ctl(paths,time_range)
    return xr.open_mfdataset(paths_,*args,**kwargs)

#...........................................................................


def parse_ctl(ctlfile, time_range=None):
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
                raise XRctlError('Invalid TDEF record: '+line)

    # Consistency check
    # -----------------
    if dset is None or nt is None:
        raise XRctlError('<%s> does not seem to be a valid GrADS control file'%ctlfile)
    else:
        if '^' in dset:
            dirn = os.path.dirname(ctlfile)
            dset = dset.replace('^',dirn+'/')
    if template is False:
        raise XRctlError('<%s> does not seem to be templated'%ctlfile)

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
        raise XRctlError('invalid time step <%s>'%dt) 

    if 'mo' in dt:
        dt = relativedelta(months=+mons)
    else:
        dt = timedelta(seconds=secs)

    lm = int(nt)
    
    if time_range is None:
        tbeg = _gat2dt(t0)
        tend = tbeg + (lm-1) * dt
    else:
        tbeg, tend = time_range
        
    # Create file list
    # ----------------
    Files = []
    t = tbeg
    while t<=tend:
        Files.append(_strTemplate(dset,time=t))
        t += dt

    return np.unique(Files)



__Months__ = ['JAN','FEB','MAR','APR','MAY','JUN',
              'JUL','AUG','SEP','OCT','NOV','DEC']


def _strTemplate(templ,expid=None,nymd=None,nhms=None,
                    yy=None,mm=None,dd=None,h=None,m=None,s=None,
                    time=None):
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

       time ---  python datetime

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

    if time is not None:
        yy = time.year
        mm = time.month
        dd = time.day
        h  = time.hour
        m  = time.minute
        s  = time.second

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

if __name__ == "__main__":

    ctlfile = '/Users/adasilva/data/merra2/ctl/tavg1_2d_aer_Nx.ctl'

    tbeg, tend = datetime(2023,4,7,0,30), datetime(2023,4,15,23,30)
    
    Files_all = parse_ctl(ctlfile)

    Files = parse_ctl(ctlfile, time_range=(tbeg,tend) )

    #print('Full Month\n', Files_all)
    #print('Partial Month\n', Files)

    #ds1 = open_mfdataset(Files,parallel=True)

    ds2 = open_mfdataset(ctlfile,time_range=(tbeg,tend),parallel=True)

 

    

