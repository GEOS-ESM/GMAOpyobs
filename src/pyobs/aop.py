"""
    Implements calculation of aerosol optical properties based on gridded GOCART
mixing ratio files (aer_Nv) and GEOSmie optics tables.

"""

__version__ = '1.0.0'

import numpy  as np
import xarray as xr
import yaml

from . import mietable  as mt
from . import xrctl     as xc

from .constants import MAPL_GRAV as GRAV

# Default YAML file mapping GOCART tracers in aer_Nv and the optics files
# -----------------------------------------------------------------------
G2G_MieMap = """
#
# GEOS Aerosol Mie table Definition for each of species.
# The order of the tracers correspond to the bins in the optics netcdf files.
#

DU:
  monoFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/optics_DU.v15_3.nc4
  bandFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/opticsBands_DU.v15_3.RRTMG.nc4
  tracers:
    - DU001
    - DU002
    - DU003
    - DU004
    - DU005

SS:
  monoFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/optics_SS.v3_3.nc4
  bandFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/opticsBands_SS.v3_3.RRTMG.nc4
  tracers:
    - SS001
    - SS002
    - SS003
    - SS004
    - SS005

OC:
  monoFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/optics_OC.v1_3.nc4
  bandFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/opticsBands_OC.v1_3.RRTMG.nc4 
  tracers:
    - OCPHOBIC
    - OCPHILIC

BC:
  monoFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/optics_BC.v1_3.nc4
  bandFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/opticsBands_BC.v1_3.RRTMG.nc4
  tracers:
    - BCPHOBIC
    - BCPHILIC

BR:
  monoFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/optics_BRC.v1_5.nc4
  bandFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/opticsBands_BRC.v1_5.RRTMG.nc4
  tracers:
    - BRPHOBIC
    - BRPHILIC

SU:
  monoFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/optics_SU.v1_3.nc4
  bandFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/opticsBands_SU.v1_3.RRTMG.nc4
  tracers:
    - SO4
    
NI:
  monoFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/optics_NI.v2_5.nc4
  bandFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/opticsBands_NI.v2_5.RRTMG.nc4
  tracers:
    - NO3AN1
    - NO3AN2
    - NO3AN3

"""

# Not all parameters in the MieTables are supported here because of complex mixing rules
# (these will be implemented as needed.) Use the MIETABLE class directly for single
# species intensive properties.
#

class AOPError(Exception):
    """
    Defines general exception errors.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class G2GAOP(object):

    def __init__ (self,aerFiles,config=None,mieRootDir=None,band=False,verbose=False):
        """
        Lazy loads GOCART mixing ratio *=(aer_NV) files and corresponding Mie tables.

        aerFiles:  str, list, or Dataset with aerosol tracers
        config:    str or YAML file handle with Mie Table file names and mapping of GOCART
                   variables to specific bins in Mie Tables. If None, uses internal default.
        mieRootDir: str, prepend string to mieTable file names.
        band:       bool, by default monochromatic tables are loaded.
                    If band is True, tables for radiation bands will be loaded instead.

        """

        if config is None:
            config = G2G_MieMap
            
        self.verbose = verbose
        if isinstance(aerFiles,xr.Dataset):
            self.aer = aerFiles
        else:
            self.aer = xc.open_mfdataset(aerFiles)

        # Load YAML Config File
        # ---------------------
        self.mieTable = yaml.safe_load(config)

        if mieRootDir is None:
            edir = ''
        else:
            edir = mieRootDir + '/'

        # Band or monochromatic files
        # ---------------------------
        for s in self.mieTable:
            m = self.mieTable[s]
            if band:
                m['mie'] = mt.MIETABLE(edir+m['bandFile'])
            else:
                m['mie'] = mt.MIETABLE(edir+m['monoFile'])

        # Check consistency of Mie tables accross species
        # -----------------------------------------------
        dims =  self.mieTable['DU']['mie'].getDims() 
        self.vector = True
        self.p, self.m = (0,0)
        for s in self.mieTable:
           dims_ = self.mieTable[s]['mie'].getDims() # dimensions of Mie Tables, a dict
           if dims_['p'] is None:
               self.vector = False  # if some species is missing pmom, cannot do vector RT
               print('Warning: PMOM is missing for '+s)
               self.p, self.m = None, None
               break
           if self.vector and dims_['p'] != dims['p']:
               self.vector = False # variable size phase matrix not yet implemented
               self.p, self.m = None, None
               print('Warning: cannot handle variable size phase matrix for PMOM')
               break
           self.p = max(self.p,dims_['p']) # max number of entries in phase matrix
           self.m = max(self.m,dims_['m']) # max number of moments in phase matrix
        
    def getAOPrt(self,Species=None,wavelength=None,vector=False):
        """
        Returns an xarray Dataset with (aot,ssa,g) if vector is
        False, otherwise (aot,ssa,g,pmon) if vector is True.

        Species:  None, str, or list. If None, all species on file,
                  otherwise subset of species.
                  
        Wavelength: float, wavelength in nm.

        vector:     bool, whether to return full phase matrix or
                    asymmetry parameter.

        """

        # Tables must have be consistent across species
        # ---------------------------------------------
        if vector:
            if not self.vector:
                print('Warning: will not calculate PMOM because of inconsistent Mie Tables.')
                vector = False            
        
        # All species on file or a subset
        # -------------------------------
        if Species is None:
            Species = list(self.mieTable.keys())
        if isinstance(Species,str):
            Species = [Species,]

        a = self.aer    # aerosol mixing ratio tracers

        # GEOS files can be inconsistent when it comes to case
        # ----------------------------------------------------
        try:
            dp = a['DELP'] 
        except:
            dp = a['delp']

        # Handy arrays for extensive properties
        # -------------------------------------
        rhodz = dp / GRAV
        dz = rhodz / a['AIRDENS']       # column thickness
        rh = a['RH']
        
        # Relevant dimensions
        # -------------------
        space = rh.shape
        aot, sca, g = np.zeros(space), np.zeros(space), np.zeros(space)
        if vector:
            ns = np.prod(space)
            p, m = self.p, self.m
            pmom = np.zeros((ns,p,m)) # flatten space dimensions for convenience
        
        for s in Species:   # loop over species

            if self.verbose:
                print('[] working on',s)
            
            Tracers = self.mieTable[s]['tracers']
            mie = self.mieTable[s]['mie']

            bin = 1
            for q in Tracers:

                if self.verbose:
                    print('   -',q)

                q_mass = rhodz * a[q]
                aot_   = mie.getAOP('aot',bin, rh, q_mass, wavelength).values
                sca_   = mie.getAOP('sca',bin, rh, q_mass, wavelength).values
                g_     = mie.getAOP('g',bin, rh, q_mass, wavelength).values

                aot += aot_
                sca += sca_    
                g   += sca_ * g_ 

                if vector:
                    
                    pmom_ = mie.getAOP('pmom', bin, rh, q_mass=q_mass,
                                        wavelength=wavelength)
                    p_, m_ = pmom_.shape[-2:]
                    pmom_ = pmom_.values.reshape((ns,p_,m_))
                    
                    pmom[:,:,:m_] += pmom_[:,:,:] # If species have fewer moments, pad wih zeros
                else:

                    g   += sca_ * g_ 

                    
                bin += 1
                
        # Final normalization of SSA and g
        # --------------------------------
        ssa = sca / aot   
        if vector:
             pmom = pmom / sca.reshape((ns,1,1))
             pmom = pmom.reshape(space+(p,m))
        else:
             g = g / sca      


        A = dict (AOT = {'long_name':'Aerosol Optical Thickness', 'units':'1'},
                  SSA = {'long_name':'Aerosol Single Scattering Albedo', 'units':'1'},
                  G = {'long_name':'Aerosol Asymmetry Parameter', 'units':'1'},
                  PMOM = {'long_name':'Aerosol Phase Matrix (non-zero elements)', 'units':'1'}
                  )

        # Pack results into a Dataset
        # ---------------------------
        DA = dict( AOT = xr.DataArray(aot,dims=rh.dims,coords=rh.coords,attrs=A['AOT']),
                   SSA = xr.DataArray(ssa,dims=rh.dims,coords=rh.coords,attrs=A['SSA']),
                 )

        DA['DELP'] = dp
        DA['AIRDENS'] = a['AIRDENS']

        if vector:
            coords = dict(rh.coords).copy()
            coords['p'] = mie.ds.coords['p']
            dims = space + ('p', 'm')
            DA['PMOM'] = xr.DataArray(pmom, dims=rh.dims+('p','m'),coords=coords)
        else:
            DA['G'] = xr.DataArray(g,dims=rh.dims,coords=rh.coords)
            
        return xr.Dataset(DA)
     
    def getAOPext(self,Species=None,wavelength=None):
        """
        Returns an xarray Dataset with the following variables:

        EXT:     aerosol extinction profile
        SCA:     aerosol scattering profile
        BSC:     aerosol backscatter profile
        DEPOL:   aerosol depolarization ratio

        On inout,
        
        Species:    None, str, or list. If None, all species on file,
                    otherwise subset of species.
                  
        Wavelength: float, wavelength in nm.

        TO DO: total attenuated backscatter, including molecular component

        """

        # All species on file or a subset
        # -------------------------------
        if Species is None:
            Species = list(self.mieTable.keys())
        if isinstance(Species,str):
            Species = [Species,]

        a = self.aer    # aerosol mixing ratio tracers

        # GEOS files can be inconsistent when it comes to case
        # ----------------------------------------------------
        try:
            dp = a['DELP'] 
        except:
            dp = a['delp']

        rh = a['RH']

        # Relevant dimensions
        # -------------------
        space = rh.shape

        ext, sca, bsc, depol1, depol2 = (np.zeros(space), np.zeros(space),
                                         np.zeros(space), np.zeros(space),
                                         np.zeros(space))

        for s in Species:   # species

            if self.verbose:
                print('[] working on',s)
            
            Tracers = self.mieTable[s]['tracers']
            mie = self.mieTable[s]['mie']

            bin = 1
            for q in Tracers:

                if self.verbose:
                    print('   -',q)

                
                q_conc = (a['AIRDENS'] * a[q]).values
                ext_ = mie.getAOP('bext', bin, rh, wavelength=wavelength).values
                sca_ = mie.getAOP('bsca', bin, rh, wavelength=wavelength).values
                bsc_ = mie.getAOP('bbck', bin, rh, wavelength=wavelength).values
                p11_ = mie.getAOP('p11',  bin, rh, wavelength=wavelength).values
                p22_ = mie.getAOP('p22',  bin, rh, wavelength=wavelength).values

                ext_ = ext_ * q_conc
                sca_ = sca_ * q_conc
                bsc_ = bsc_ * q_conc

                ext += ext_
                sca += sca_
                bsc += bsc_
                depol1 += (p11_-p22_) * sca_
                depol2 += (p11_+p22_) * sca_

                bin += 1
                
        # Final normalization
        # -------------------
        ext *= 1000. # m-1 to km-1
        sca *= 1000. # m-1 to km-1
        bsc *= 1000. # m-1 to km-1
        depol = depol1 / depol2

        # Attributes
        # ----------
        A = dict (EXT = {'long_name':'Aerosol Extinction Coefficient', 'units':'km-1'},
                  SCA = {'long_name':'Aerosol Scattering Coefficient', 'units':'km-1'},
                  BSC = {'long_name':'Aerosol Backscatter Coefficient', 'units':'km-1'},
                  DEPOL = {'long_name':'Depolarization Ratio', 'units':'1'}
                  )
        
        # Pack results into a Dataset
        # ---------------------------
        DA = dict(  EXT = xr.DataArray(ext.astype('float32'),dims=rh.dims,coords=rh.coords,attrs=A['EXT']),
                    SCA = xr.DataArray(sca.astype('float32'),dims=rh.dims,coords=rh.coords,attrs=A['SCA']),
                    BSC = xr.DataArray(bsc.astype('float32'),dims=rh.dims,coords=rh.coords,attrs=A['BSC']),
                    DEPOL = xr.DataArray(depol.astype('float32'),dims=rh.dims,coords=rh.coords,attrs=A['DEPOL'])
                 )

        DA['DELP'] = dp
        DA['AIRDENS'] = a['AIRDENS']
        
        return xr.Dataset(DA)
     

    def getAOPintensive(self,Species=None,wavelength=None):
        """
        Returns an xarray Dataset with intensive properties.

        Species:  None, str, or list. If None, all species on file,
                  otherwise subset of emissions.
                  
        Wavelength: float, wavelength in nm. 

        """
        raise AOPError("not implemented yet")
        
    
    def getAOP(self,what,Species=None,wavelength=None):
        """
        Returns an xarray Dataset with the aerosol aerosol optical
        property requested.

        what:     str, list with AOPs to calculate
        Species:  None, str, or list. If None, all species on file,
                  otherwise subset of emissions.
        Wavelength: float, wavelength in nm. 
        
        """

        raise AOPError("not implemented yet")


#....................................................................................

def CLI_aop():
    """
    Parses command line and write files with Aerosol Optical Properties.
    """

    import sys
    import os
    
    from optparse        import OptionParser

    format = 'NETCDF4'
    config = None   # use internal config.
    outYAML = 'aop.yaml'
    outFile = 'aop_%{w}nm.nc4'
    aop = 'ext'
    rootDir = './'
    wavelengths='550'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] aerDataset [iso_t1 iso_t2]\n"+\
                                "   aerDataset          GrADS-style ctl or a shell-style wildcard string\n"+\
                                "                       with aerosol mixing ratios, either gridded or sampled"+\
                                "   iso_t1,iso_t2       optional beginning and ending time (ISO format)",
                          version=__version__ )

    parser.add_option("-a", "--aop", dest="aop", default='ext',
              help="AOP collection, one of 'rt' or'ext' (default=%s)"%aop)
    
    parser.add_option("-c", "--config", dest="config", default=None,
              help="optional configuration YAML file (default='buit-in')")

    parser.add_option("-d", "--dump",
                      action="store_true", dest="dump",
                      help="Dumps internal YAML configuration to stdout and stops.")

    parser.add_option("-f", "--format", dest="format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC or NETCDF3_64BIT (default=%s)"%format )

    parser.add_option("-o", "--output", dest="outFile", default=outFile,
              help="Output NetCDF file name; use %%{w} as a placeholder for wavelength (default=%s)"\
                          %outFile )

    parser.add_option("-r", "--root", dest="rootDir", default=rootDir,
              help="Root directory for MieTables (default=%s)"\
                          %rootDir )

    parser.add_option("-V", "--vector",
                      action="store_true", dest="vector",
                      help="Vector mode.")

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")

    parser.add_option("-w", "--wavelengths", dest="wavelengths", default=wavelengths,
              help="Comma separated wavelengths (default=%s)"\
                          %wavelengths )

    (options, args) = parser.parse_args()
    
    if options.dump:
        print(G2G_MieMap)
        sys.exit(0)

    if len(args) == 1:
        aerDataset = args[0]
        t1, t2 = None, None
    elif len(args) == 3:
        aerDataset, t1, t2 = args
        t1, t2 = None, None
    else:
        parser.error("must have 1 or 3 arguments: aerDataset [iso_t1 iso_t2]")

        
    # Create consistent file name extension
    # -------------------------------------
    name, ext = os.path.splitext(options.outFile)
    if 'NETCDF4' in options.format:
        options.outFile = name + '.nc4'
    elif 'NETCDF3' in options.format:
        options.outFile = name + '.nc'
    else:
        raise ValueError('Invalid extension <%s>'%ext)

    if options.config is not None:
        config = open(options.config,'r')
    else:
        config = None

    # Compute AOPs
    # ------------
    aer = xc.open_mfdataset(aerDataset,parallel=True) 
    g = G2GAOP(aer,config=config,mieRootDir=options.rootDir,verbose=options.verbose)
    for w_ in options.wavelengths.split(','):
        w = float(w_)
        if options.aop == 'ext':
            ds = g.getAOPext(wavelength=w)
        elif options.aop == 'rt':
            ds = g.getAOPrt(wavelength=w,vector=options.vector)
        else:
            print(options.aop)
            raise AOPError('Unknow AOP option '+options.aop)

        filename = options.outFile.replace('%{w}',w_)
        if options.verbose:
            print('Writing',filename)
        ds.to_netcdf(filename)  # TO DO: Chunking and compression
    
#....................................................................................
def Test_g2g_aop():
    """
    Simple tests.
    """

    # yaml.dump(rc,open('test.yml','w'))

    data = '/Users/adasilva/data/'
    aer_Nv = '/Users/adasilva/data/sampled/aer_Nv/CAMP2Ex-GEOS-MODISonly-aer-Nv-P3B_Model_*.nc'
    
    aer = xr.open_mfdataset(aer_Nv) # still having trouble with parallel

    g = G2GAOP(aer,mieRootDir=data,verbose=True)
    rts = None # g.getAOPrt(wavelength=550,vector=False)
    rtv = g.getAOPrt(wavelength=550,vector=True)
    ext = None # g.getAOPext(wavelength=550)

    return (g, rts, rtv, ext)

if __name__ == "__main__":

    g, rts, rtv, ext = Test_g2g_aop



    

    

    
