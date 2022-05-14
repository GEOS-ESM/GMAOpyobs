#!/usr/bin/env python

"""
    Utility to sample GEOS-5  files at fixed station locations.

"""

import os

from numpy import zeros, ones, arange, array

from optparse        import OptionParser
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser
from csv             import DictReader

from MAPL           import Config
from MAPL.constants import *

class StnVar(object):
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
        var = StnVar(v)
        var.title = f.vtitle[i]
        var.km = f.kmvar[i]
        try:
            var.units = f.vunits[i]
        except:
            var.units = 'unknown'  # opendap currently lacks units
        Vars[v] = var

    f.Vars = Vars

    return f

def getStations(csvFile):
    """
    Parse CSV file.
    """
    CSV = DictReader(open(csvFile))
    name, lon, lat = [], [], []
    for row in CSV:
        name += [row['name'],]
        lon  += [float(row['lon']),]
        lat  += [float(row['lat']),]
    return ( array(name), array(lon), array(lat) )

def getTyme(f,t1,t2):
    tbeg, tend = f.tbracket(t1)[0], f.tbracket(t2)[1]
    if tend>t2:
        tend -= f.dt
    t = tbeg
    tyme = [t,]
    while t < tend:
        t += f.dt
        tyme += [t,]
    tyme = array(tyme)
    return tyme

#---
def stnSample(f,V,stnLon,stnLat,tyme,options,squeeze=True):
    """
    Sample file at station locations.
    """
    if options.verbose:
        print " [] Interpolating <%s>"%V.name
    ns, nt, nz = (len(stnLon), len(tyme), V.km)
    if f.lower:
        name = V.name.lower() # GDS always uses lower case
    else:
        name = V.name
    if nz>0:
        Z = zeros((ns,nt,nz))
    else:
        Z = zeros((ns,nt))
    n = 0
    for t in tyme:
        try:
            z = f.interp(name,stnLon,stnLat,tyme=t,algorithm=options.algo,
                         Transpose=True,squeeze=squeeze)
        except:
            print "    - Interpolation failed for <%s> on %s"%(V.name,str(t))
            if nz>0:
                z = MAPL_UNDEF * ones((ns,nz))
            else:
                z = MAPL_UNDEF * ones(ns)

        if nz>0:
            Z[:,n,:] = z
        else:
            Z[:,n] = z
        n += 1
        Z[abs(Z)>MAPL_UNDEF/1000.] = MAPL_UNDEF # detect undef contaminated interp
    return Z

#---
def writeXLS ( stnName, stnLon, stnLat, tyme, f, options,
               title='GEOS-5 Station Sampler'):
    """
    Write a Excel Spreadsheet file with sampled GEOS-5 variables at the
    station locations described by (lon,lat).
    """
    from xlwt import Workbook 

    ns_, nt_, nz_ = ( len(stnLon), len(tyme), f.km )

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
    meta.write(3,1,'Created from GEOS-5 standard collections by stn_sampler.py')
    
    meta.write(4,0,'References:')  
    meta.write(4,1,'n/a')  
    
    meta.write(5,0,'Comment:') 
    meta.write(5,1,'This file contains GEOS-5 parameters sampled at station locations.')

    meta.write(6,0,'Contact:')
    meta.write(6,1,'Arlindo da Silva <arlindo.dasilva@nasa.gov>')

    # Time in ISO format
    # ------------------
    isoTime = array([t.isoformat() for t in tyme])

    Vars = options.Vars or f.vname

    # Each station on a separate sheet
    # --------------------------------
    Sheets = dict()
    for stn in stnName:
        sheet = book.add_sheet(stn)
        sheet.write(1,0,'Time')
        j = 1
        for v in Vars:
            var = f.Vars[v.upper()]
            if var.km > 0:
                print 'Warning: ignoring <%s>, only single-level variables supported for now'
                continue 
            sheet.write(0,j,var.name.upper())
            sheet.write(1,j,var.title.replace('_',' ').replace('ensemble',''))
            sheet.write(2,j,var.units)
            j += 1
        Sheets[stn] = sheet
        for n in range(len(isoTime)):
            sheet.write(n+3,0,isoTime[n])
            
    # Loop over variables on GFIO file.
    # --------------------------------
    k = 0
    for v in Vars:
        var = f.Vars[v.upper()]
        if var.km == 0:
            dim = ('station','time',)
            shp = ( ns_, nt_)
        else:
            continue # only 2D for now

        # Interpolate
        # -----------
        if options.dryrun:
            this = zeros(shp)
            if options.verbose:
                print "[] Zero-filling <%s>"%var.name
        else:
            this = stnSample(f,var,stnLon,stnLat,tyme,options)

        # Loop over stations
        # ------------------
        for s in range(len(stnName)):
            sheet = Sheets[stnName[s]]
            for n in range(len(isoTime)):
                sheet.write(n+3,k+1,this[s,n])

        k += 1
        
    # Close the workbook
    # ------------------
    book.save(options.outFile)

    if options.verbose:
        print " <> wrote %s file %s"%(options.format,options.outFile)
    
#---
def writeNC ( stnName, stnLon, stnLat, tyme, f, options,
              zlib=False, title='GEOS-5 Station Sampler'):
    """
    Write a NetCDF file with sampled GEOS-5 variables at the station locations
    described by (lon,lat).
    """
    from netCDF4 import Dataset

    ns_, nt_, nz_ = ( len(stnLon), len(tyme), f.km )

    # Open NC file
    # ------------
    nc = Dataset(options.outFile,'w',format=options.format)

    # Set global attributes
    # ---------------------
    nc.Title = title
    nc.Institution = 'NASA/Goddard Space Flight Center'
    nc.Source = 'Global Model and Assimilation Office'
    nc.History = 'Created from GEOS-5 standard collections by stn_sampler.py'
    nc.References = 'n/a'
    nc.Comment = 'This file contains GEOS-5 parameters sampled at station locations.'
    nc.Contact = 'Arlindo da Silva <arlindo.dasilva@nasa.gov>'
    nc.Conventions = 'CF'

    # Because of GrADS limitations on attribute size, break Station list into
    # multiple attributes. This is Hackish, but necessary
    # ------------------------------------------------------------------------
    for i in range(ns_):
        setattr(nc,'Station_%03d'%(i+1),stnName[i])

    # Create dimensions
    # -----------------
    ns = nc.createDimension('station', ns_ )
    nz = nc.createDimension('lev', nz_ )
    nt = nc.createDimension('time', nt_ )
    ls = nc.createDimension('ls',19)
    x = nc.createDimension('x',1)
    y = nc.createDimension('y',1)
     
    # Station names
    # -------------
    stnName_ = nc.createVariable('stnName','S1',('station','ls'),zlib=zlib)
    stnName_.long_name = 'Station Names'
    stnName_.axis = 'e'
    stntmp = zeros((ns_,19),dtype='S1')
    for i in range(ns_):
        stntmp[i][:] = list('%-19s'%stnName[i])
    stnName_[:] = stntmp[:]

    # Coordinate variables
    # --------------------
    lon = nc.createVariable('stnLon','f4',('station',),zlib=zlib)
    lon.long_name = 'Longitude'
    lon.units = 'degrees_east'
    lon[:] = stnLon[:]
    
    lat = nc.createVariable('stnLat','f4',('station',),zlib=zlib)
    lat.long_name = 'Latitude'
    lat.units = 'degrees_north'
    lat[:] = stnLat[:]

    # Add fake dimensions for GrADS compatibility
    # ------------------------------------------
    x = nc.createVariable('x','f4',('x',),zlib=zlib)
    x.long_name = 'Fake Longitude for GrADS Compatibility'
    x.units = 'degrees_east'
    x[:] = zeros(1)
    y = nc.createVariable('y','f4',('y',),zlib=zlib)
    y.long_name = 'Fake Latitude for GrADS Compatibility'
    y.units = 'degrees_north'
    y[:] = zeros(1)
    e = nc.createVariable('station','i4',('station',),zlib=zlib)
    e.long_name = 'Station Ensemble Dimension'
    e.axis = 'e'
    e.grads_dim = 'e'
    e[:] = range(ns_)
    
    if f.km > 0:
        lev = nc.createVariable('lev','f4',('lev',),zlib=zlib)
        lev.long_name = 'Vertical Level'
        lev.units = f.levunits.rstrip()
        lev.positive = 'down'
        lev.axis = 'z'
        lev[:] = f.levs[:]

    time = nc.createVariable('time','i4',('time',),zlib=zlib)
    time.long_name = 'Time'
    t0 = tyme[0]
    time.units = 'seconds since %s'%t0.isoformat(' ')
    time[:] = array([(t-t0).total_seconds() for t in tyme])

    # Time in ISO format if so desired
    # ---------------------------------
    if options.isoTime:
        isotime = nc.createVariable('isotime','S1',('time','ls'),zlib=zlib)
        isotime.long_name = 'Time (ISO Format)'
        isotmp = zeros((nt_,19),dtype='S1')
        for i in range(nt_):
            isotmp[i][:] = list(tyme[i].isoformat())
        isotime[:] = isotmp[:]
      
    # Loop over variables on GFIO file.
    # --------------------------------
    Vars = options.Vars or f.vname
    for v in Vars:
        #var = f.Vars[v.upper()]
        var = f.Vars[v]
        if var.km == 0:
            dim = ('station','time',)
            shp = ( ns_, nt_)
        else:
            dim = ('station','time','lev')
            shp = ( ns_, nt_, nz_)
        this = nc.createVariable(var.name,'f4',dim,zlib=zlib)
        this.standard_name = var.title
        this.long_name = var.title.replace('_',' ')
        this.missing_value = MAPL_UNDEF
        this.units = var.units
        if options.dryrun:
            this_ = zeros(shp)
            if options.verbose:
                print "[] Zero-filling <%s>"%var.name
        else:
            this_ = stnSample(f,var,stnLon,stnLat,tyme,options)
        this[:] = this_[:]
            
    # Close the file
    # --------------
    nc.close()

    if options.verbose:
        print " <> wrote %s file %s"%(options.format,options.outFile)
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    
    format = 'NETCDF3_CLASSIC'
    outFile = 'stn_sampler.nc'
    algo = 'linear'
    
#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] stnFile g5File [iso_t1 iso_t2]",
                          version='1.0.0' )

    parser.add_option("-o", "--output", dest="outFile", default=outFile,
              help="Output NetCDF file (default=%s)"\
                          %outFile )

    parser.add_option("-a", "--algorithm", dest="algo", default=algo,
              help="Interpolation algorithm, one of linear, nearest (default=%s)"\
                          %algo)

    parser.add_option("-V", "--vars", dest="Vars", default=None,
              help="Variables to sample (default=All)")
    
    parser.add_option("-f", "--format", dest="format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC, NETCDF3_64BIT or EXCEL (default=%s)"%format )

    parser.add_option("-I", "--isoTime",
                      action="store_true", dest="isoTime",
                      help="Include time in ISO format as well.")

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")

    parser.add_option("-n", "--dryrun",
                      action="store_true", dest="dryrun",
                      help="Dry-run mode: fill variables with zeros.")


    (options, args) = parser.parse_args()
    
    if len(args) == 4 :
        stnFile, g5File, iso_t1, iso_t2 = args
        t1, t2 = (isoparser(iso_t1), isoparser(iso_t2))
    elif len(args) == 2 :
        stnFile, g5File = args
        t1, t2 = (None,None)
    else:
        parser.error("must have 2 or 4 arguments: stnFile g5File [iso_t1 iso_t2]")

    if options.Vars is not None:
        options.Vars = options.Vars.split(',')

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
        
    # Open the file and gather metadata
    # ---------------------------------
    f = Open(g5File)

    # If to time specified, then use whole file
    # -----------------------------------------
    if t1 is None:
        t1, t2 = f.tbeg, f.tend
    
    # Station Locations
    # -----------------
    stnName, stnLon, stnLat = getStations(stnFile)

    # Time range
    # ----------
    tyme = getTyme(f,t1,t2)
    
    # Write output file
    # -----------------
    if options.format == 'EXCEL':
        writeXLS (stnName,stnLon,stnLat,tyme,f,options)
    else:
        writeNC(stnName,stnLon,stnLat,tyme,f,options)

