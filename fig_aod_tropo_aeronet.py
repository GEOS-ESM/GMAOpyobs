#!/usr/bin/env python
"""
Do co-location of AERONET and MODIS 
"""

import os,sys
from   pyobs.mxd04     import MxD04_L2, MISSING, granules, BEST 
#from   mxd04_nnr       import  SDS, ALIAS, CHANNELS, SCHANNELS, TranslateInput, TranslateTarget
from   glob            import glob
from   pyobs.aeronet   import AERONET_L2
from   optparse        import OptionParser   # Command-line args
from   datetime        import datetime, timedelta
from   dateutil.parser import parse as isoparse
import numpy           as     np
from   haversine        import haversine
from   netCDF4         import Dataset, stringtoarr
from   pyobs.bits  import BITS

Ident = dict( modo = ('MOD04','ocean'),
              modl = ('MOD04','land'),
              modd = ('MOD04','deep'),
              mydo = ('MYD04','ocean'),
              mydl = ('MYD04','land'),
              mydd = ('MYD04','deep')
            )

META = ( 'ScatteringAngle','GlintAngle',
         'SolarAzimuth', 'SolarZenith',
         'SensorAzimuth','SensorZenith',
         'cloud','qa_flag'  )


for var in META:
    TranslateInput[var] = (var,)

META += ('aod','reflectance','ISO_DateTime')

SDSmod = {'OCEAN'  : META,
            'LAND' : META + ('sfc_reflectance',),
            'DEEP' : META + ('sfc_reflectance',)}


SDSanet = ('AERONET_Site','lat','lon',
           'AOT_1640' , 
           'AOT_1020', 
           'AOT_870', 
           'AOT_675', 
           'AOT_667',
           'AOT_555', 
           'AOT_551', 
           'AOT_532', 
           'AOT_531', 
           'AOT_500',
           'AOT_490'  ,
           'AOT_443' , 
           'AOT_440', 
           'AOT_412', 
           'AOT_380', 
           'AOT_340',
           'AOT_550')

SDSmerra = ('fdu','fss','fcc','fsu')


# -----------------------------------------------------------------------------
class MODIS(MxD04_L2):
    """ Does QA/QC selections for different MODIS collections """
    def __init__(self,l2_path,prod,algo,year,julday,
                 cloud_thresh=0.70,
                 glint_thresh=40.0,
                 scat_thresh=170.0,
                 coll='006',verbose=0):

#        Files = glob('{}/{}/{}/{}/{}/*hdf'.format(l2_path,coll,prod,year,str(julday).zfill(3)))        
        Files = sorted(glob('{}/{}/{}/{}/{}/*hdf'.format(l2_path,coll,prod,year,str(julday).zfill(3))))
        if algo != "DEEP":
            MxD04_L2.__init__(self,Files,algo,
                              only_good=True,
                              SDS=SDS,
                              Verb=verbose)   
        else:
            MxD04_L2.__init__(self,Files,algo,
                              only_good=True,
                              SDS=SDS,                            
                              alias=ALIAS,
                              Verb=verbose)   


        if self.nobs < 1:
            return # no obs, nothing to do

        # Reorganize Reflectance Arrays
        # -----------------------------
        self.rChannels = CHANNELS[algo]
        if algo in SCHANNELS:
            self.sChannels = SCHANNELS[algo]

        if algo == "OCEAN":
            self.reflectance = self.reflectance[:,0:7]  #not using 412, 443, and 745 for now
        if algo == "LAND":
            self.reflectance = self.reflectance[:,0:-1]  #not using 745 for now


        # 3-Ch Algorithm only used when Dark Target data is unavailable
        # --------------------------------------------------------------

        if algo == "DEEP":
            # Get DARK TARGET qa_flag
            self.qa_flag_lnd = BITS(self.Quality_Assurance_Land[:,0])[1:4]            
            lndGood = self.qa_flag_lnd == BEST
            lndGood = lndGood & (self.cloud_lnd < cloud_thresh)
            rChannels = CHANNELS["LAND"]
            sChannels = SCHANNELS["LAND"]
            for i,c in enumerate(rChannels):
                lndGood = lndGood & (self.reflectance_lnd[:,i]>0)

            for i,c in enumerate(sChannels):
                lndGood = lndGood & (self.sfc_reflectance_lnd[:,i]>0)

            self.iGood = (self.qa_flag == BEST) & ~lndGood

#            # Keep only "good" observations
#            # -----------------------------
#            m = self.iGood
#            for sds in self.SDS:
#                rank = len(self.__dict__[sds].shape)
#                if rank == 1:
#                    self.__dict__[sds] = self.__dict__[sds][m]
#                elif rank == 2:
#                    self.__dict__[sds] = self.__dict__[sds][m,:]
#                else:
#                    raise IndexError, 'invalid rank=%d'%rank
#
#            # Reset aliases
#            for sds in self.SDS:
#                if sds in self.ALIAS:
#                    self.__dict__[self.ALIAS[sds]] = self.__dict__[sds] 
#
#
#            self.qa_flag = self.qa_flag[m]
#            self.aod     = self.aod[m,:]
#            self.Time    = self.Time[m]
#            self.iGood   = self.iGood[m] 
#            self.nobs    = self.Longitude.shape[0]         
#
#            if self.nobs < 1:
#                return # no obs, nothing to do             


        # Q/C
        # ---        
        self.iGood = self.cloud<cloud_thresh 
        for i,c in enumerate(self.rChannels):
            self.iGood = self.iGood & (self.reflectance[:,i]>0)

        if algo in SCHANNELS:
            for i,c in enumerate(self.sChannels):
                self.iGood = self.iGood & (self.sfc_reflectance[:,i]>0)

        if algo == "OCEAN":
            self.iGood = self.iGood & (self.GlintAngle > glint_thresh)

        if algo != "OCEAN":
            self.iGood = self.iGood & (self.ScatteringAngle < scat_thresh)

        if any(self.iGood) == False:
            print("WARNING: Strange, no good obs left to work with")
            print("Setting nobs to zero")
            self.nobs=0
            return

# -----------------------------------------------------------------------------

class AERONET(AERONET_L2):
    """ Gets AERONET data """

    def __init__(self,l2_path,date,verbose=0):

        MM = str(date.month).zfill(2)
        nymd = str(date.date()).replace('-','')
        Files = glob('{}/Y{}/M{}/aeronet_v20.{}.txt'.format(l2_path,date.year,MM,nymd))
        
        AERONET_L2.__init__(self,Files,Verbose=verbose)

        if self.nobs > 0:
            iGood = self.AOT_550 >= 0

            self.reduce(iGood)
# -----------------------------------------------------------------------------

def coarseColocate(anet,mod):
    # For each unique AERONET site
    # Loop through and look for possible distance & time matches
    # Coarse initial search

    usites = np.unique(anet.AERONET_Site)
    dtMax = timedelta(minutes=60)
    dlMax = 0.5
    matches = []

    for site in usites:
        isite = anet.AERONET_Site == site
        lon,lat = anet.lon[isite][0], anet.lat[isite][0]

        distfound = (np.abs(lon - mod.Longitude) < dlMax) & (np.abs(lat - mod.Latitude) < dlMax)
        distfound = distfound & mod.iGood

        times = anet.tyme[isite]
        timefound = np.zeros(times.shape)
        for i,t in enumerate(times):
            tmin = t - dtMax
            tmax = t + dtMax

            timefound[i] = np.any((mod.Time[distfound] <= tmax) & (mod.Time[distfound] >= tmin))

        if np.any(timefound):
            matches.append(np.arange(mod.nobs)[distfound])
        else:
            matches.append([])

    return matches
# -----------------------------------------------------------------------------


def fineColocate(anet,mod,Cmatches):
    # For each aeronet site that has coarse matches
    # loops thrgh and looks for distance and time matches

    usites = np.unique(anet.AERONET_Site)
    dtMax = timedelta(minutes=30)
    dlMax = 27.5 #km
    modMatches = []
    anetMatches = []

    
    for site,mm in zip(usites,Cmatches):
        if any(mm):
            isite = anet.AERONET_Site == site

            # Get MODIS obs within a dlMax circle of AERONET site
            alat,alon = anet.lat[isite][0],anet.lon[isite][0]
            mlat,mlon = mod.lat[mm], mod.lon[mm]
            distance = haversine(alat,alon,mlat,mlon)

            found = distance <= dlMax
            modMatches.append(mm[found])

            if not any(found):
                anetMatches.append([])
            else:
                # Get AERONET obs within +/- dtMax of median MODIS overpass
                atyme = anet.tyme[isite]
                mdT   = (mod.Time[mm[found]].max() - mod.Time[mm[found]].min()).total_seconds()
                mtyme = mod.Time[mm[found]].min() + timedelta(seconds=mdT*0.5)

                tmin = mtyme - dtMax
                tmax = mtyme + dtMax
                found = (atyme <= tmax) & (atyme >= tmin)
                anetMatches.append(np.arange(anet.nobs)[isite][found])

        else:
            modMatches.append([])
            anetMatches.append([])


    return modMatches, anetMatches    

# -----------------------------------------------------------------------------


def createGiant(mod,anet,options,algo):
    if not os.path.exists(options.outpath):
        os.makedirs(options.outpath)
    outFile = options.outpath + '/{}_{}_giant.nc4'.format(options.ident,options.coll)


    ncOut  = Dataset(outFile,'w',format='NETCDF4_CLASSIC',clobber=True)
    ncOut.comment  = 'MODIS-AERONET colocation, follows MAPSS protocal'

    obs = ncOut.createDimension('obs',None)
    length = ncOut.createDimension('length',100)
    channels = ncOut.createDimension('channels',len(mod.channels))
    rChannels = ncOut.createDimension('rChannels',len(mod.rChannels))
    sChannels = ncOut.createDimension('sChannels',len(mod.sChannels))

    VARS = SDSmod[algo.upper()]
    for varname in VARS:
        if varname == 'ISO_DateTime':
            varobj = ncOut.createVariable(varname,'c',('obs','length',))
        elif varname == 'aod':
            varobj = ncOut.createVariable(varname,'f4',('obs','channels',))
        elif varname == 'reflectance':
            varobj = ncOut.createVariable(varname,'f4',('obs','rChannels',))
        elif varname == 'sfc_reflectance':
            varobj = ncOut.createVariable(varname,'f4',('obs','sChannels',))
        else:
            varobj = ncOut.createVariable(varname,'f4',('obs',))

    VARS = SDSanet
    for varname in VARS:
        if varname == 'AERONET_Site':
            varobj = ncOut.createVariable(varname,'c',('obs','length',))        
        else:
            varobj = ncOut.createVariable(varname,'f4',('obs',))

    VARS = SDSmerra
    for varname in VARS:
        varobj = ncOut.createVariable(varname,'f4',('obs',))    

    varobj = ncOut.createVariable('channels','f4',('channels',))        
    varobj[:] = mod.channels

    varobj = ncOut.createVariable('rChannels','f4',('rChannels',))
    varobj[:] = mod.rChannels

    varobj = ncOut.createVariable('sChannels','f4',('sChannels'))
    varobj[:] = mod.sChannels

    ncOut.close()



# -----------------------------------------------------------------------------

if __name__ == '__main__':
    
    # Defaults
    aer_x        = '/nobackup/NNR/Misc/tavg1_2d_aer_Nx'
    outpath  = 'giants/'

    l2_path  = '/nobackup/MODIS/Level2/'
    ident    = 'mydl'
    coll     = '006'

    anet_path = '/nobackup/AERONET/Level2/'

    isotime = '2015-06-01'
    DT      = 3  
    dType   = 'days'
    verbose = False
    append  = False


#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("-x", "--aer_x", dest="aer_x", default=aer_x,
                      help="MERRA2 aer_Nx file (default=%s)"\
                           %aer_x )  

    parser.add_option("-o", "--outpath", dest="outpath", default=outpath,
                      help="Outpath (default=%s)"\
                           %outpath )    

    parser.add_option("-l", "--l2_path", dest="l2_path", default=l2_path,
                      help="MODIS Level2 Path (default=%s)"\
                           %l2_path )

    parser.add_option("-a", "--anet_path", dest="anet_path", default=anet_path,
                      help="AERONET Level2 Path (default=%s)"\
                           %anet_path )    

    parser.add_option("-i", "--ident", dest="ident", default=ident,
                      help="Algorithm code (default=%s)"\
                           %ident )

    parser.add_option("-I", "--isotime", dest="isotime", default=isotime,
                      help="isotime (default=%s)"\
                           %isotime )  

    parser.add_option("-D", "--DT", dest="DT", default=DT,
                      help="delta time (default=%s)"\
                           %DT ) 

    parser.add_option("-t", "--dType", dest="dType", default=dType,
                      help="DT type (default=%s)"\
                           %dType )     

    parser.add_option("-c", "--coll", dest="coll", default=coll,
                      help="Collection (default=%s)"\
                           %coll )    

    parser.add_option("-v", "--verbose", dest="verbose", default=verbose,action="store_true",
                      help="verbose (default=%s)"\
                           %verbose )       

    parser.add_option("-A", "--append", dest="append", default=append,action="store_true",
                      help="append mode - don't overwrite (default=%s)"\
                           %append )                                                          

    (options, args) = parser.parse_args()
    
    prod, algo = Ident[options.ident]

    nymd     = isoparse(options.isotime)
    if options.dType == 'days':
        enymd =  nymd + timedelta(days=int(options.DT))
    if options.dType == 'months':
        enymd = nymd + timedelta(months=options.DT)

    while nymd < enymd:
        print('++++Working on ',nymd)
        julday   = nymd - datetime(nymd.year,1,1) + timedelta(days=1)
        mod = MODIS(options.l2_path,prod,algo.upper(),nymd.year,julday.days,
                      coll=options.coll,
                      cloud_thresh=0.7,
                      verbose=options.verbose)
        

        anet = AERONET(options.anet_path,nymd,verbose=options.verbose)

        if (anet.nobs >0) & (mod.nobs >0):
            Cmatches = coarseColocate(anet,mod)

            modMatches, anetMatches  = fineColocate(anet,mod,Cmatches)

            if any(map(any,anetMatches)):
                print('Has Data',nymd)

                # Create giant file if you need to
                if not options.append:
                    createGiant(mod,anet,options,algo)
                    options.append =True
                                     
                else:
                    outFile = options.outpath + '/{}_{}_giant.nc4'.format(options.ident,options.coll)
                    if not os.path.exists(outFile):
                        createGiant(mod,anet,options,algo)

                # Create class of colocated data and append to giant file
                co = COLOCATE(mod,anet,modMatches,anetMatches,options,algo)
                co.writeGiant(options)                



        nymd += timedelta(days=1)


    
