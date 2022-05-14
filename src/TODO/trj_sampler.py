#!/usr/bin/env python

"""
    Utility to generate a orbital trajectory and sample model
    output along it.

    Arlindo da Silva, February 2014.

"""

import os

from numpy import zeros, arange, array, float32

from optparse        import OptionParser
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser
from csv             import DictReader

from MAPL           import Config, eta
from MAPL.constants import *
from pyobs.sgp4     import getTrack as getTrackTLE
from pyobs          import ICARTT, NPZ, HSRL, ORACLES

class TleVar(object):
    """                                                                                  
    Generic container for Variables
    """
    def __init__(self,name):
        self.name = name


#---
def Open(filename):
    """
    Uses GFIO or GFIOctl to open either a NetCDF-4 or a control file.
    Very heuristic.
    """
    from gfio import GFIOurl, GFIOctl
    name, ext = os.path.splitext(filename)

    # Open the GFIO dataset
    # ---------------------
    f = None
    if 'HTTP://' == name[:7].upper():
        f = GFIOurl(filename)
        f.lower = True # force variable names to be lower case when sampling.
    elif ext.upper() in ('.NC4','.NC','.HDF','.H5'):
        f = GFIOurl(filename)
        f.lower = False
    else:
        f = GFIOctl(filename)
        f.lower = False

    # Create variable dictionary
    # --------------------------
    Vars = dict()
    if len(f.vtitle)<len(f.vname):
        f.vtitle = f.vname[:]      # In case vtitle is not filled (hack)
    for i in range(len(f.vname)):
        if f.lower:
            v = f.vname[i].upper()
        else:
            v = f.vname[i]
        var = TleVar(v)
        var.title = f.vtitle[i]
        var.km = f.kmvar[i]
        if var.km>0:
            var.levunits = f.levunits[:]
            var.levs = f.levs[:]    
        try:
            var.units = f.vunits[i]
        except:
            var.units = 'unknown'  # opendap currently lacks units
        Vars[v] = var

    f.Vars = Vars

    return f

def getTrackICT(ictFile,dt_secs):
    """
    Get trajectory from ICART (.ict) file.
    """
    m = ICARTT(ictFile)
    lon, lat, tyme = m.Nav['Longitude'], m.Nav['Latitude'], m.Nav['Time']
    mdt = (tyme[-1] - tyme[0]).total_seconds()/float(len(tyme)-1) # in seconds
    idt = int(dt_secs/mdt+0.5)
    return (lon[::idt], lat[::idt], tyme[::idt])

def getTrackHSRL(hsrlFile,dt_secs=60):
    """
    Get trajectory from HSRL HDF-5 file.
    """
    h = HSRL(hsrlFile,Nav_only=True)
    lon, lat, tyme = h.lon[:].ravel(), h.lat[:].ravel(), h.tyme[:].ravel()
    if dt_secs > 0:
        dt = tyme[1] - tyme[0] 
        idt = int(dt_secs/dt.total_seconds()+0.5)
        return (lon[::idt], lat[::idt], tyme[::idt])
    else:
        idt = 1
    return (lon[::idt], lat[::idt], tyme[::idt])

def getTrackORACLES(mergeFile,dt_secs=60):
    """
    Get trajectory from ORACLES merge file (netCDF).
    """
    m = ORACLES(mergeFile)
    lon, lat, tyme = m.lon, m.lat, m.tyme
    mdt = (tyme[-1] - tyme[0]).total_seconds()/float(len(tyme)-1) # in seconds
    idt = int(dt_secs/mdt+0.5)
    return (lon[::idt], lat[::idt], tyme[::idt])

def getTrackCSV(csvFile):
    """
    Get trajectory from a CSV with (lon,lat,time) coordinates.


                  2014-02-05T12:30:45
    """
    CSV = DictReader(open(csvFile))
    lon, lat, tyme = [], [], []
    for row in CSV:
        lon  += [float(row['lon']),]
        lat  += [float(row['lat']),]
        tyme  += [isoparser(row['time']),]
        
    return ( array(lon), array(lat), array(tyme) )
    
def getTrackNPZ(npzFile):
    """
    Get trajectory from a NPZ with (lon,lat,time) coordinates.
    Notice that *time* is a datetime object.

    Note: These are simple NPZ usually generated during Neural
          Net or other type of python based utility. Not meant
          for general consumption, but could be since NPZ files
          are much more compact than CSV.

    """
    n = NPZ(npzFile)
    if 'time' in n.__dict__:
        return ( n.lon, n.lat, n.time)
    elif 'tyme' in n.__dict__:
        return ( n.lon, n.lat, n.tyme)
    else:
        raise ValueError, 'NPZ file has neither *time* nor *tyme* attribut.e'

def getVars(rcFile):
    """
    Parse reource file, create variable dictionary with relevant
    metadata.
    """

    cf = Config(rcFile)
    Vars = dict()
    AllVars = dict()
    levUnits = 'none'
    levs = []
    for V in cf.keys():
        path = cf(V)
        f = Open(path)
        varList = []
        if '*' in V:
            VARS = f.Vars.keys()
        else:
            VARS = V.split(',') 
        for v in VARS:
            v = v.strip()
            var = f.Vars[v]
            if AllVars.__contains__(v): # unique variable names for output
                print " >< Skipping duplicate variable <%s> in %s"%(v,path)
                continue
            elif v.upper() == "TAITIME":
                continue # annoying HDFEOS crap
            else:
                AllVars[v] = True
            if var.km>0:
                levUnits = var.levunits
                levs = var.levs  
            varList += [var,]
        Vars[path] = varList
        
    return (Vars, levs, levUnits)

#---
def writeNC ( lons, lats, tyme, Vars, levs, levUnits, trjFile, options,
              title='GEOS-5 Trajectory Sampler',
              doAkBk=False, zlib=False):
    """
    Write a NetCDF file with sampled GEOS-5 variables along the satellite track
    described by (lon,lat,tyme).
    """
    from netCDF4 import Dataset

    km = len(levs)
    
    # Open NC file
    # ------------
    nc = Dataset(options.outFile,'w',format=options.format)

    # Set global attributes
    # ---------------------
    nc.title = title
    nc.institution = 'NASA/Goddard Space Flight Center'
    nc.source = 'Global Model and Assimilation Office'
    nc.history = 'Created from GEOS-5 standard collections by trj_sampler.py'
    nc.references = 'n/a'
    nc.comment = 'This file contains GEOS-5 related parameters along a satellite or aircraft track.'
    nc.contact = 'Arlindo da Silva <arlindo.dasilva@nasa.gov>'
    nc.Conventions = 'CF'
    nc.trjFile = trjFile
 
    # Create dimensions
    # -----------------
    nt = nc.createDimension('time',len(tyme))
    ls = nc.createDimension('ls',19)
    if km>0:
        nz = nc.createDimension('lev',km)
        if doAkBk:
            ne = nc.createDimension('ne',km+1)
    x = nc.createDimension('x',1)
    y = nc.createDimension('y',1)

    # Coordinate variables
    # --------------------
    time = nc.createVariable('time','i4',('time',),zlib=zlib)
    time.long_name = 'Time'
    t0 = tyme[0]
    time.units = 'seconds since %s'%t0.isoformat(' ')
    time[:] = array([int((t-t0).total_seconds()+0.5) for t in tyme])
    if km > 0: # pressure level not supported yet
        lev = nc.createVariable('lev','f4',('lev',),zlib=zlib)
        lev.long_name = 'Vertical Level'
        lev.units = levUnits.strip()
        lev.positive = 'down'
        lev.axis = 'z'
        lev[:] = levs[:]

        if doAkBk:
            ae, be = eta.getEdge(km) # Coefficients for Hybrid coordinates
            ak = nc.createVariable('ak','f4',('ne',),zlib=zlib)
            ak.long_name = 'Eta coordinate coefficient ak (p = ak + bk * ps)'
            ak.units = 'Pa'
            ak = ae[:]
            bk = nc.createVariable('bk','f4',('ne',),zlib=zlib)
            bk.long_name = 'Eta coordinate coefficient bk (p = ak + bk * ps)'
            bk.units = '1'
            bk = be[:]
    
    # Add fake dimensions for GrADS compatibility
    # -------------------------------------------
    x = nc.createVariable('x','f4',('x',),zlib=zlib)
    x.long_name = 'Fake Longitude for GrADS Compatibility'
    x.units = 'degrees_east'
    x[:] = zeros(1)
    y = nc.createVariable('y','f4',('y',),zlib=zlib)
    y.long_name = 'Fake Latitude for GrADS Compatibility'
    y.units = 'degrees_north'
    y[:] = zeros(1)
    
    # Trajectory coordinates
    # ----------------------
    lon = nc.createVariable('trjLon','f4',('time',),zlib=zlib)
    lon.long_name = 'Trajectory Longitude'
    lon.units = 'degrees_east'
    lon[:] = lons[:]
    lat = nc.createVariable('trjLat','f4',('time',),zlib=zlib)
    lat.long_name = 'Trajectory Latitude'
    lat.units = 'degrees_north'
    lat[:] = lats[:]
    
    # Time in ISO format if so desired
    # ---------------------------------
    if options.isoTime:
        isotime = nc.createVariable('isotime','S1',('time','ls'),zlib=zlib)
        isotime.long_name = 'Time (ISO Format)'
        isotmp = zeros((len(lons),19),dtype='S1')
        for i in range(len(lons)):
            isotmp[i][:] = list(tyme[i].isoformat())
        isotime[:] = isotmp[:]
      
    # Loop over datasets, sample and write each variable
    # --------------------------------------------------
    for path in Vars:
        if options.verbose:
            print " <> opening "+path
        g = Open(path) 
        for var in Vars[path]:
            if var.km == 0:
                dim = ('time',)
            else:
                dim = ('time','lev')
            this = nc.createVariable(var.name,'f4',dim,zlib=zlib)
            this.standard_name = var.title
            this.long_name = var.title.replace('_',' ')
            this.missing_value = float32(MAPL_UNDEF)
            this.units = var.units
            if g.lower:
                name = var.name.lower() # GDS always uses lower case
            else:
                name = var.name
            if options.verbose:
                print " [] %s interpolating <%s>"%\
                        (options.algo.capitalize(),name.upper())        
            Z = g.sample(name,lons,lats,tyme,algorithm=options.algo,
                         Transpose=True,squeeze=True)
            Z[abs(Z)>MAPL_UNDEF/1000.] = MAPL_UNDEF # detect undef contaminated interp
            this[:] = Z
            
    # Close the file
    # --------------
    nc.close()

    if options.verbose:
        print " <> wrote %s file %s"%(options.format,options.outFile)
    
#---
def writeXLS ( lons, lats, tyme, Vars, trjFile, options,
              title='GEOS-5 Trajectory Sampler'):
    """
    Write a Excel Spreadsheet file with sampled GEOS-5 variables along the satellite track
    described by (lons,lats,tyme).
    """
    from xlwt import Workbook 

    km = len(levs)
    
    # Open XLS file
    # -------------
    book = Workbook()
    meta = book.add_sheet('Metadata')
    
    # Set global attributes on its own sheet
    # --------------------------------------
    meta.write(0,0,'Title:')
    meta.write(0,1,title)

    meta.write(1,0,'Institution:')
    meta.write(1,1,'NASA/Goddard Space Flight Center')

    meta.write(2,0,'Source:')
    meta.write(2,1,'Global Model and Assimilation Office')
    
    meta.write(3,0,'History:')
    meta.write(3,1,'Created from GEOS-5 standard collections by trj_sampler.py')
    
    meta.write(4,0,'References:')  
    meta.write(4,1,'n/a')  
    
    meta.write(5,0,'Comment:') 
    meta.write(5,1,'This file contains GEOS-5 parameters sampled along a satellite or aircraft track.')

    meta.write(6,0,'Contact:')
    meta.write(6,1,'Arlindo da Silva <arlindo.dasilva@nasa.gov>')

    meta.write(7,0,'Trajectory File:')
    meta.write(7,1,trjFile)
 
    # Time in ISO format
    # ------------------
    isoTime = array([t.isoformat() for t in tyme])

    # Data sheet
    # ----------
    sheet = book.add_sheet('Data')

    # Header: coordinates
    # -------------------
    sheet.write(1,0,'Time')
    sheet.write(1,1,'Longitude')
    sheet.write(2,1,'degrees')
    sheet.write(1,2,'Latitude')
    sheet.write(2,2,'Degrees')
    
    # Coordinate variables
    # --------------------
    for n in range(len(isoTime)):
        sheet.write(n+3,0,isoTime[n])
        sheet.write(n+3,1,lons[n])
        sheet.write(n+3,2,lats[n])

    # Loop over datasets, sample and write each variable
    # --------------------------------------------------
    j = 3
    for path in Vars:
        
        if options.verbose:
            print " <> opening "+path
        g = Open(path) 
        for var in Vars[path]:
            if var.km > 0:
                print 'Warning: ignoring <%s>, only single-level variables supported for now'%var.name
                continue # no profiles for now
            if g.lower:
                name = var.name.lower() # GDS always uses lower case
            else:
                name = var.name
                
            # Variable header
            # ---------------    
            sheet.write(0,j,var.name.upper())
            sheet.write(1,j,var.title.replace('_',' ').replace('ensemble',''))
            sheet.write(2,j,var.units)

            # Interpolate
            # -----------
            if options.verbose:
                print " [] Interpolating <%s>"%name.upper()
            Z = g.sample(name,lons,lats,tyme,Transpose=True,squeeze=True)
            Z[abs(Z)>MAPL_UNDEF/1000.] = MAPL_UNDEF # detect undef contaminated interp

            # Write to sheet
            # --------------
            Z = Z.astype('float')
            for n in range(len(isoTime)):
                sheet.write(n+3,j,Z[n])

            j += 1
                
    # Close the file
    # --------------
    book.save(options.outFile)

    if options.verbose:
        print " <> wrote %s file %s"%(options.format,options.outFile)
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    
    format = 'NETCDF3_CLASSIC'
    rcFile  = 'trj_sampler.rc'
    outFile = 'trj_sampler.nc'
    dt_secs = 60
    algo = 'linear'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] tleFile|ictFile|csvFile|npzFile [iso_t1 iso_t2]",
                          version='1.0.1' )

    parser.add_option("-a", "--algorithm", dest="algo", default=algo,
              help="Interpolation algorithm, one of linear, nearest (default=%s)"\
                          %algo)

    parser.add_option("-r", "--rcFile", dest="rcFile", default=rcFile,
              help="Resource file defining parameters to sample (default=%s)"\
                          %rcFile )

    parser.add_option("-o", "--output", dest="outFile", default=outFile,
              help="Output NetCDF file (default=%s)"\
                          %outFile )

    parser.add_option("-f", "--format", dest="format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC, NETCDF3_64BIT or EXCEL (default=%s)"%format )

    parser.add_option("-t", "--trajectory", dest="traj", default=None,
                      help="Trajectory file format: one of tle, ict, csv, npz, hsrl, oracles (default=trjFile extension)" )

    parser.add_option("-d", "--dt_secs", dest="dt_secs", default=dt_secs,
              type='int',
              help="Timesetp in seconds for TLE sampling (default=%s)"%dt_secs )

    parser.add_option("-I", "--isoTime",
                      action="store_true", dest="isoTime",
                      help="Include ISO format time in output file.")

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")

    (options, args) = parser.parse_args()
    
    if len(args) == 3:
        trjFile, iso_t1, iso_t2 = args
        t1, t2 = (isoparser(iso_t1), isoparser(iso_t2))
    elif len(args) == 1:
        trjFile = args[0]
        t1, t2 = None, None
    else:
        parser.error("must have 1 or 3 arguments: tleFile|ictFile [iso_t1 iso_t2]")

    if options.traj is None:
        name, ext = os.path.splitext(trjFile)
        options.traj = ext[1:]
    options.traj = options.traj.upper()
        
    # Create consistent file name extension
    # -------------------------------------
    name, ext = os.path.splitext(options.outFile)
    if ext.upper() == '.XLS':
        options.format = 'EXCEL'
    if 'NETCDF4' in options.format:
        options.outFile = name + '.nc4'
    elif 'NETCDF3' in options.format:
        options.outFile = name + '.nc'
    elif 'EXCEL' in options.format:
        options.outFile = name + '.xls'
    else:
        raise ValueError, 'invalid extension <%s>'%ext
  
    # Get Variables and Metadata
    # --------------------------
    Vars, levs, levUnits = getVars(options.rcFile)

    # Create trajectory
    # -----------------
    if options.traj == 'TLE':
        if t1 is None:
            raise ValueError, 'time range (t1,t2) must be specified when doing TLE sampling.'
        lon, lat, tyme = getTrackTLE(trjFile, t1, t2, options.dt_secs)
    elif options.traj == 'ICT':
        lon, lat, tyme = getTrackICT(trjFile,options.dt_secs)
    elif options.traj == 'CSV':
        lon, lat, tyme = getTrackCSV(trjFile)
    elif options.traj == 'NPZ':
        lon, lat, tyme = getTrackNPZ(trjFile)
    elif options.traj == 'HSRL' or options.traj == 'H5':
        lon, lat, tyme = getTrackHSRL(trjFile,options.dt_secs)
    elif options.traj == 'ORACLES':
        lon, lat, tyme = getTrackORACLES(trjFile,options.dt_secs)
    else:
        raise ValueError, 'cannot handle trajectory file format <%s>'%options.traj

    # Make sure longitudes in [-180,180]
    # ----------------------------------
    if lon.max()>180.:
        lon[lon>180] = lon[lon>180] - 360.

    # Write output file
    # -----------------
    if options.format == 'EXCEL':
        writeXLS(lon,lat,tyme,Vars,trjFile,options)
    else:
        writeNC(lon,lat,tyme,Vars,levs,levUnits,trjFile,options)
    
