#!/usr/bin/env python

"""
    Utility to compute optical properties along a sampled file.

    Ed Nowottnick, January, 2015.

    Modified to return lidar intensive properties
    P. Castellanos June 2017

"""

import os
import sys
import MieObs_
from netCDF4 import Dataset
from mieobs import getAOPext, getAOPint, getAOPscalar, getEdgeVars
import numpy as np
from math import pi, sin, cos, asin, acos

from optparse        import OptionParser
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser
from csv             import DictReader

from MAPL           import Config, eta
from MAPL.constants import *
from pyobs          import ICARTT
from pyobs.sgp4     import getTrack 

# Generic Lists of Varnames and Units
VNAMES_DU = ['DU001','DU002','DU003','DU004','DU005']
VNAMES_SS = ['SS001','SS002','SS003','SS004','SS005']
VNAMES_BC = ['BCPHOBIC','BCPHILIC']
VNAMES_OC = ['OCPHOBIC','OCPHILIC']
VNAMES_SU = ['SO4']
VNAMES_NI = ['NO3AN1','NO3AN2','NO3AN3']
VNAMES_BRC = ['BRCPHOBIC','BRCPHILIC']

MieVarsNames = ['ext','scatext','backscat','aback_sfc','aback_toa','depol','ext2back','tau','ssa','g']
MieVarsUnits = ['km-1','km-1','km-1 sr-1','sr-1','sr-1','unitless','sr','unitless','unitless','unitless']
MieVarsLongNames = ['total aerosol extinction','scattering extinction','total aerosol backscatter',
                    'attenuated aerosol backscatter at the surface','attenuated aerosol backscatter at TOA',
                    'depolarization ratio','extinction to backscatter ratio','aerosol optical depth',
                    'single scattering albedo','assymetry parameter']

IntVarsUnits = {'area':'m2 m-3', 
                'vol': 'm3 m-3',
                'refi':'unitless',
                'refr':'unitless',
                'reff':'m'}

IntVarsLongNames = {'area':'aerosol cross sectional area',
                    'vol' :'aerosol volume',
                    'refi':'imaginary refractive index',
                    'refr':'real refractive index',
                    'reff':'aerosol effective raidus'}                

class MieCalc(object):
    pass
    """                                                                                  
    Generic container for Variables
    """

#---
def getVars(inFile):
    """
#    Parse input file, create variable dictionary
#    """

    Vars       = MieCalc()
    file       = Dataset(inFile)
    names      = file.variables.keys()
    MIENAMES   = names
    for n, name in enumerate(MIENAMES):
        var = file.variables[name]
        if (name == 'trjLon') or (name == 'stnLon'):
            name = 'LONGITUDE'
        if (name == 'trjLat') or (name == 'stnLat'):
            name = 'LATITUDE'
        name = name.upper()
        size = len(var.shape)
        if size == 3:
            setattr(Vars,name,var[:,:,:])
        if size == 2:
            setattr(Vars,name,var[:,:])
        if size == 1:
            setattr(Vars,name,var[:])
    return Vars        

#---
def computeMie(Vars, channel, varnames, rcFile, options):
    """
#    Computes optical quantities and combines into a dictionary
#   """

    #STN Sampled?
    if options.station:
        NAMES = varnames + ['PS','DELP','RH','AIRDENS']
        nstn = len(Vars.STATION)
        nobs = len(Vars.TIME)

        for v in range(nstn):
            VarsIn = MieCalc()
            for n, name in enumerate(NAMES):
                Var = getattr(Vars,name)
                size = len(Var.shape)
                if (size == 2):
                    #1D Variables, ex. PS
                    setattr(VarsIn,name,Var[v,:])
                if size == 3:
                    #2D Variables, ex. DU001
                    setattr(VarsIn,name,Var[v,:,:])

  
            if (v==0):
                pe, ze, te = getEdgeVars(VarsIn)
                tau,ssa,g = getAOPscalar(VarsIn,channel,vnames=varnames,vtypes=varnames,Verbose=options.verbose,rcfile=rcFile)
                ext,sca,backscat,aback_sfc,aback_toa,depol = getAOPext(VarsIn,channel,I=None,vnames=varnames,vtypes=varnames,Verbose=options.verbose,rcfile=rcFile)
                ext2back = np.ones(backscat.shape)*MAPL_UNDEF
                I = backscat > 0
                ext2back[I] = ext[I]/backscat[I]
                MieVars = {"ext":[ext],"scatext":[sca],"backscat":[backscat],"aback_sfc":[aback_sfc],"aback_toa":[aback_toa],"depol":[depol],"ext2back":[ext2back],"tau":[tau],"ssa":[ssa],"g":[g]}

                if options.intensive:
                    vol, area, refr, refi, reff = getAOPint(VarsIn,channel,I=None,vnames=varnames,vtypes=varnames,Verbose=options.verbose,rcfile=rcFile)
                    MieVars['vol']  = [vol]
                    MieVars['area'] = [area]
                    MieVars['refr'] = [refr]
                    MieVars['refi'] = [refi]
                    MieVars['reff'] = [reff]
                
                MieVars['pe'] = [pe]
                MieVars['ze'] = [ze]
                MieVars['rh'] = [np.transpose(VarsIn.RH)]         


            else:
                pe, ze, te = getEdgeVars(VarsIn)
                tau,ssa,g = getAOPscalar(VarsIn,channel,vnames=varnames,vtypes=varnames,Verbose=options.verbose,rcfile=rcFile)
                ext,sca,backscat,aback_sfc,aback_toa,depol = getAOPext(VarsIn,channel,I=None,vnames=varnames,vtypes=varnames,Verbose=options.verbose,rcfile=rcFile)
                ext2back = np.ones(backscat.shape)*MAPL_UNDEF
                I = backscat > 0
                ext2back[I] = ext[I]/backscat[I]
                MieVars['ext'].append(ext)
                MieVars['scatext'].append(sca)
                MieVars['backscat'].append(backscat)
                MieVars['aback_sfc'].append(aback_sfc)
                MieVars['aback_toa'].append(aback_toa)
                MieVars['depol'].append(depol) 
                MieVars['ext2back'].append(ext2back)
                MieVars['tau'].append(tau)
                MieVars['ssa'].append(ssa)
                MieVars['g'].append(g)

                MieVars['pe'].append(pe)
                MieVars['ze'].append(ze)
                MieVars['rh'].append(np.transpose(VarsIn.RH))

                if options.intensive:
                    vol, area, refr, refi, reff = getAOPint(VarsIn,channel,I=None,vnames=varnames,vtypes=varnames,Verbose=options.verbose,rcfile=rcFile)
                    MieVars['vol'].append(vol)
                    MieVars['area'].append(area)
                    MieVars['refr'].append(refr)
                    MieVars['refi'].append(refi)
                    MieVars['reff'].append(reff)                    

    #TRJ Sampled?
    else:
        pe, ze, te = getEdgeVars(Vars)
        tau,ssa,g = getAOPscalar(Vars,channel,vnames=varnames,vtypes=varnames,Verbose=options.verbose,rcfile=rcFile)
        ext,sca,backscat,aback_sfc,aback_toa,depol = getAOPext(Vars,channel,I=None,vnames=varnames,vtypes=varnames,Verbose=options.verbose,rcfile=rcFile)
        ext2back = np.ones(backscat.shape)*MAPL_UNDEF
        I = backscat > 0
        ext2back[I] = ext[I]/backscat[I]
        MieVars = {"ext":[ext],"scatext":[sca],"backscat":[backscat],"aback_sfc":[aback_sfc],"aback_toa":[aback_toa],"depol":[depol],"ext2back":[ext2back],"tau":[tau],"ssa":[ssa],"g":[g]}       

        MieVars['pe'] = [pe]
        MieVars['ze'] = [ze]
        MieVars['rh'] = [np.transpose(Vars.RH)]
        if options.intensive:
            vol, area, refr, refi, reff  = getAOPint(Vars,channel,I=None,vnames=varnames,vtypes=varnames,Verbose=options.verbose,rcfile=rcFile)
            MieVars['vol']  = [vol]
            MieVars['area'] = [area]
            MieVars['refr'] = [refr]
            MieVars['refi'] = [refi]
            MieVars['reff'] = [reff]

    return MieVars

#---
def writeNC ( stations, lons, lats, tyme, isotimeIn, MieVars, MieVarsNames, MieVarsLongNames,
              MieVarsUnits, inFile, outFile, options, zlib=False):
    """
    Write a NetCDF file with sampled GEOS-5 variables along the satellite track
    described by (lon,lat,tyme).
    """
    km = 72

    # Open NC file
    # ------------
    nc = Dataset(outFile,'w',format=options.format)

    # Set global attributes
    # ---------------------
    nc.title = 'GEOS-5 Sampled Aerosol Optical Properties File'
    nc.institution = 'NASA/Goddard Space Flight Center'
    nc.source = 'Global Model and Assimilation Office'
    nc.history = 'Created from sampled GEOS-5 collections'
    nc.references = 'n/a'
    nc.comment = 'This file contains sampled GEOS-5 aerosol optical properties.'
    fixrh = float(options.fixrh)
    if(fixrh >= 0.):
       nc.comment = nc.comment+' Calculations carried out at fixed RH = '+str(fixrh)
    nc.contact = 'Ed Nowottnick <edward.p.nowottnick@nasa.gov>'
    nc.Conventions = 'CF'
    nc.inFile = inFile
 
    # Create dimensions
    # -----------------
    nt = nc.createDimension('time',len(tyme))
    if options.station:
        ns = nc.createDimension('station',len(stations))
    ls = nc.createDimension('ls',19)
    if km>0:
        nz = nc.createDimension('lev',km)
        nze = nc.createDimension('leve',km+1)
    x = nc.createDimension('x',1)
    y = nc.createDimension('y',1)

    if options.station:
        # Station names
        # -------------
        stnName_ = nc.createVariable('stnName','S1',('station','ls'),zlib=zlib)
        stnName_.long_name = 'Station Names'
        stnName_.axis = 'e'
        stnName_[:] = stations[:]   

    # Coordinate variables
    # --------------------
    time = nc.createVariable('time','i4',('time',),zlib=zlib)
    time.long_name = 'Time'
    t0 = tyme[0]
    isot0 = isotimeIn[0]
    date0 = ''.join(isot0[:10])
    time0 = ''.join(isot0[-8:])
    time.units = 'seconds since '+date0+' '+time0
    time[:] = tyme
    if km > 0: # pressure level not supported yet
        lev = nc.createVariable('lev','f4',('lev',),zlib=zlib)
        lev.long_name = 'Vertical Level'
        lev.units = 'none'
        lev.positive = 'down'
        lev.axis = 'z'
        lev[:] = range(1,km+1)

        lev = nc.createVariable('leve','f4',('leve',),zlib=zlib)
        lev.long_name = 'Vertical Level Edge'
        lev.units = 'none'
        lev.positive = 'down'
        lev.axis = 'z'
        lev[:] = range(1,km+2)

    # Add fake dimensions for GrADS compatibility
    # -------------------------------------------
    x = nc.createVariable('x','f4',('x',),zlib=zlib)
    x.long_name = 'Fake Longitude for GrADS Compatibility'
    x.units = 'degrees_east'
    x[:] = np.zeros(1)
    y = nc.createVariable('y','f4',('y',),zlib=zlib)
    y.long_name = 'Fake Latitude for GrADS Compatibility'
    y.units = 'degrees_north'
    y[:] = np.zeros(1)
    if options.station:
        e = nc.createVariable('station','i4',('station',),zlib=zlib)
        e.long_name = 'Station Ensemble Dimension'
        e.axis = 'e'
        e.grads_dim = 'e'
        e[:] = range(len(stations))
    
    # Lat/Lon Coordinates
    # ----------------------
    if options.station:
        lon = nc.createVariable('longitude','f4',('station',),zlib=zlib)
        lon.long_name = 'Longitude'
        lon.units = 'degrees_east'
        lon[:] = lons[:]
        lat = nc.createVariable('latitude','f4',('station',),zlib=zlib)
        lat.long_name = 'Latitude'
        lat.units = 'degrees_north'
        lat[:] = lats[:]
    else:
        lon = nc.createVariable('longitude','f4',('time',),zlib=zlib)
        lon.long_name = 'Longitude'
        lon.units = 'degrees_east'
        lon[:] = lons[:]
        lat = nc.createVariable('latitude','f4',('time',),zlib=zlib)
        lat.long_name = 'Latitude'
        lat.units = 'degrees_north'
        lat[:] = lats[:]        
    
    # Time in ISO format if so desired
    # ---------------------------------
    isotime = nc.createVariable('isotime','S1',('time','ls'),zlib=zlib)
    isotime.long_name = 'Time (ISO Format)'
    isotime[:] = isotimeIn[:]

    # Write each variable
    # --------------------------------------------------
    for n, name in enumerate(MieVarsNames):
        
        var = np.squeeze(MieVars[name])
        size = len(var.shape)
        if options.station:
            # make sure profiles are always 3 dimensions
            # even if only 1 station or time interval
            if size == 1:
                var.shape = (len(stations),km,len(tyme))
            if size == 2:
                if len(tyme) == 1:
                    var.shape = (len(stations),km,len(tyme))
            size = len(var.shape)
        if size == 3:
            dim = ('station','time','lev')
        if size == 2:
            dim = ('time','lev')
        if size == 1:
            dim = ('time')
        this = nc.createVariable(name,'f4',dim,zlib=zlib)
        this.standard_name = name
        this.long_name = MieVarsLongNames[n]
        this.units = MieVarsUnits[n]
        this.missing_value = np.float32(MAPL_UNDEF)
        if options.station:
            this[:] = np.transpose(var,(0,2,1))
        else:
            this[:] = np.transpose(var)

    uu = ['Pa','m','none']
    long_name = ['Pressure','Altitude','Relative humidity']
    for n,name in enumerate(['pe','ze','rh']):
        var = np.squeeze(MieVars[name])
        if options.station:
            # make sure profiles are always 3 dimensions
            # even if only 1 station or time interval
            if name == 'rh':
                var.shape = (len(stations),km,len(tyme))
            else:
                var.shape = (len(stations),km+1,len(tyme))

            if name == 'rh':
                dim = ('station','time','lev')
            else:
                dim = ('station','time','leve')
        else:
            if name == 'rh':
                dim = ('time','lev')
            else:
                dim = ('time','leve')       

        this = nc.createVariable(name,'f4',dim,zlib=zlib)
        this.standard_name = name
        this.long_name = long_name[n]
        this.units = uu[n]
        this.missing_value = np.float32(MAPL_UNDEF)
        if options.station:
            this[:] = np.transpose(var,(0,2,1))
        else:
            this[:] = np.transpose(var)

    if options.intensive:
        for name in IntVarsUnits:
            var = np.squeeze(MieVars[name])
            size = len(var.shape)
            if options.station:
                if size == 2:
                    if len(tyme) == 1:
                        var.shape = (len(stations),1,km)
                    else:
                        var.shape = (1,len(tyme),km)
                    size = len(var.shape)
            if size == 3:
                dim = ('station','time','lev')
            if size == 2:
                dim = ('time','lev')
            if size == 1:
                dim = ('time')
            this = nc.createVariable(name,'f4',dim,zlib=zlib)
            this.standard_name = name
            this.long_name     = IntVarsLongNames[name]
            this.units = IntVarsUnits[name]
            this.missing_value = np.float32(MAPL_UNDEF)
            if options.station:
                this[:] = np.transpose(var,(0,2,1))
            else:
                this[:] = np.transpose(var)            

    # Close the file
    # --------------
    nc.close()

    if options.verbose:
        print " <> wrote %s file %s"%(options.format,options.outFile)
    
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    
    format = 'NETCDF3_CLASSIC'
    inFile  = 'trj_sampler.nc'
    outFile = 'ext_sampler.nc'
    channel = (532)
    rcFile = 'Aod3d_532nm.rc'
    chmFile = 'Chem_MieRegistry.rc'
    intensive = False
    fixrh   = (-1.)

#   Parse command line options
#   --------------------------
    parser = OptionParser()

    parser.add_option("-i", "--input", dest="inFile", default=inFile,
              help="Sampled input file (default=%s)"%inFile)

    parser.add_option("-o", "--output", dest="outFile", default=outFile,
              help="Output file containing optical properties (default=%s)"%outFile)

    parser.add_option("-r", "--rc", dest="rcFile", default=rcFile,
              help="Resource file pointing to optical tables (default=%s)"%rcFile)

    parser.add_option("-C", "--chm", dest="chmFile", default=chmFile,
              help="Chem Registry Resource file (default=%s)"%chmFile)

    parser.add_option("-f", "--format", dest="format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC, NETCDF3_64BIT or EXCEL (default=%s)"%format )

    parser.add_option("-c", "--channel", dest="channel", default=channel,
              help="Channel for Mie calculation (default={})".format(channel))

    parser.add_option("--rh", dest="fixrh", default=fixrh,
              help="If specified use provided RH (0-1) in calculations")

    parser.add_option("--vnames", dest="VNAMES", default=None,
              help="Comma sepearted list of species to include in calculation (default=read from chem registry file)")    

    parser.add_option("-I", "--intensive",default=intensive,
                      action="store_true", dest="intensive",
                      help="return intensive variables")

    parser.add_option("-v", "--verbose",default=False,
                      action="store_true", dest="verbose",
                      help="Verbose mode")

    parser.add_option("--du",
                      action="store_true", dest="dust",
                      help="Output separate file with DUST only properties")

    parser.add_option("--ss",
                      action="store_true", dest="seasalt",
                      help="Output separate file with SEASALT only properties")

    parser.add_option("--su",
                      action="store_true", dest="sulfate",
                      help="Output separate file with SULFATE only properties")

    parser.add_option("--bc",
                      action="store_true", dest="bcarbon",
                      help="Output separate file with BLACK CARBON only properties")

    parser.add_option("--oc",
                      action="store_true", dest="ocarbon",
                      help="Output separate file with ORGANIC CARBON only properties")

    parser.add_option("--ni",
                      action="store_true", dest="nitrate",
                      help="Output separate file with NITRATE only properties")

    parser.add_option("--brc",
                      action="store_true", dest="brcarbon",
                      help="Output separate file with BROWN CARBON only properties")

    parser.add_option("--stn",
                      action="store_true", dest="station",
                      help="Input File is from stn_sampler.py")

    (options, args) = parser.parse_args()

         
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
    
    # Get Variables
    # --------------------------
    Vars = getVars(options.inFile)

    # Check FIXRH option
    # --------------------------
    fixrh = float(options.fixrh)
    if(fixrh >= 0.):
        if(fixrh > 1.):
           print "Your --RH > 1, must be between 0 - 1; exit and fix"
           sys.exit()
        Vars.RH = Vars.RH*0.0+fixrh

    # Variables to be included in calculation
    if options.VNAMES is None:
        VNAMES = []
        cf = Config(options.chmFile)
        try:
            du = cf('doing_DU')
            if du.upper() == 'YES':
                VNAMES += VNAMES_DU
        except:
            pass
        try:
            ss = cf('doing_SS')
            if ss.upper() == 'YES':
                VNAMES += VNAMES_SS
        except:
            pass
        try: 
            su = cf('doing_SU')
            if su.upper() == 'YES':
                VNAMES += VNAMES_SU
        except:
            pass
        try:
            oc = cf('doing_OC')
            if oc.upper() == 'YES':
                VNAMES += VNAMES_OC
        except:
            pass
        try:
            bc = cf('doing_BC')
            if bc.upper() == 'YES':
                VNAMES += VNAMES_BC
        except:
            pass
        try:
            brc = cf('doing_BRC')
            if brc.upper() == 'YES':
                VNAMES += VNAMES_BRC
        except:
            pass
        try:
            ni = cf('doing_NI')
            if ni.upper() == 'YES':
                VNAMES += VNAMES_NI
        except:
            pass

    if type(options.VNAMES) is str:
        VNAMES = options.VNAMES.split(',')

    # Run Mie Calculator and Write Output Files
    # --------------------------
    if options.station:
        StnNames = Vars.STNNAME
    else:
        StnNames = ''

    channelIn = float(options.channel)
    MieVars = computeMie(Vars,channelIn,VNAMES,options.rcFile,options)
    writeNC(StnNames,Vars.LONGITUDE,Vars.LATITUDE,Vars.TIME,Vars.ISOTIME,
            MieVars,MieVarsNames,MieVarsLongNames,MieVarsUnits,options.inFile,options.outFile,options)

    name, ext = os.path.splitext(options.outFile)
    if options.dust:
        outFile = name+'.du'+ext
        MieVars = computeMie(Vars,channelIn,VNAMES_DU,options.rcFile,options)
        writeNC(StnNames,Vars.LONGITUDE,Vars.LATITUDE,Vars.TIME,Vars.ISOTIME,
                MieVars,MieVarsNames,MieVarsLongNames,MieVarsUnits,options.inFile,outFile,options)

    if options.seasalt:
        outFile = name+'.ss'+ext
        MieVars = computeMie(Vars,channelIn,VNAMES_SS,options.rcFile,options)
        writeNC(StnNames,Vars.LONGITUDE,Vars.LATITUDE,Vars.TIME,Vars.ISOTIME,
                MieVars,MieVarsNames,MieVarsLongNames,MieVarsUnits,options.inFile,outFile,options)

    if options.sulfate:
        outFile = name+'.su'+ext
        MieVars = computeMie(Vars,channelIn,VNAMES_SU,options.rcFile,options)
        writeNC(StnNames,Vars.LONGITUDE,Vars.LATITUDE,Vars.TIME,Vars.ISOTIME,
                MieVars,MieVarsNames,MieVarsLongNames,MieVarsUnits,options.inFile,outFile,options)

    if options.bcarbon:
        outFile = name+'.bc'+ext
        MieVars = computeMie(Vars,channelIn,VNAMES_BC,options.rcFile,options)
        writeNC(StnNames,Vars.LONGITUDE,Vars.LATITUDE,Vars.TIME,Vars.ISOTIME,
                MieVars,MieVarsNames,MieVarsLongNames,MieVarsUnits,options.inFile,outFile,options)

    if options.ocarbon:
        outFile = name+'.oc'+ext
        MieVars = computeMie(Vars,channelIn,VNAMES_OC,options.rcFile,options)
        writeNC(StnNames,Vars.LONGITUDE,Vars.LATITUDE,Vars.TIME,Vars.ISOTIME,
                MieVars,MieVarsNames,MieVarsLongNames,MieVarsUnits,options.inFile,outFile,options)
   
    if options.nitrate:
        outFile = name+'.ni'+ext
        MieVars = computeMie(Vars,channelIn,VNAMES_NI,options.rcFile,options)
        writeNC(StnNames,Vars.LONGITUDE,Vars.LATITUDE,Vars.TIME,Vars.ISOTIME,
                MieVars,MieVarsNames,MieVarsLongNames,MieVarsUnits,options.inFile,outFile,options)

    if options.brcarbon:
        outFile = name+'.brc'+ext
        MieVars = computeMie(Vars,channelIn,VNAMES_BRC,options.rcFile,options)
        writeNC(StnNames,Vars.LONGITUDE,Vars.LATITUDE,Vars.TIME,Vars.ISOTIME,
                MieVars,MieVarsNames,MieVarsLongNames,MieVarsUnits,options.inFile,outFile,options)
