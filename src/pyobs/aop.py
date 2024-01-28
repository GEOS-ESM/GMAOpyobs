"""
    Implements calculation of aerosol optical properties based on gridded GOCART
mixing ratio files (aer_Nv) and GEOSmie optics tables.

"""

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
  monoFile: ExtData/chemistry/AerosolOptics/v0.0.0/x/optics_DU.v15_3.nc
  bandFile: ExtData/chemistry/AerosolOptics/v0.0.0/x/opticsBands_DU.v15_3.RRTMG.nc
  tracers:
    - DU001
    - DU002
    - DU003
    - DU004
    - DU005

SS:
  monoFile: ExtData/chemistry/AerosolOptics/v0.0.0/x/optics_SS.v3_3.nc
  bandFile: ExtData/chemistry/AerosolOptics/v0.0.0/x/opticsBands_SS.v3_3.RRTMG.nc
  tracers:
    - SS001
    - SS002
    - SS003
    - SS004
    - SS005

OC:
  monoFile: ExtData/chemistry/AerosolOptics/v0.0.0/x/optics_OC.v1_3.nc
  bandFile: ExtData/chemistry/AerosolOptics/v0.0.0/x/opticsBands_OC.v1_3.RRTMG.nc 
  tracers:
    - OCPHOBIC
    - OCPHILIC

BC:
  monoFile: ExtData/chemistry/AerosolOptics/v0.0.0/x/optics_BC.v1_3.nc
  bandFile: ExtData/chemistry/AerosolOptics/v0.0.0/x/opticsBands_BC.v1_3.RRTMG.nc
  tracers:
    - BCPHOBIC
    - BCPHILIC

BR:
  monoFile: ExtData/chemistry/AerosolOptics/v0.0.0/x/optics_BRC.v1_5.nc
  bandFile: ExtData/chemistry/AerosolOptics/v0.0.0/x/opticsBands_BRC.v1_5.RRTMG.nc
  tracers:
    - BRPHOBIC
    - BRPHILIC

SU:
  monoFile: ExtData/chemistry/AerosolOptics/v0.0.0/x/optics_SU.v1_3.nc
  bandFile: ExtData/chemistry/AerosolOptics/v0.0.0/x/opticsBands_SU.v1_3.RRTMG.nc
  tracers:
    - SO4
    
NI:
  monoFile: ExtData/chemistry/AerosolOptics/v0.0.0/x/optics_NI.v2_5.nc
  bandFile: ExtData/chemistry/AerosolOptics/v0.0.0/x/opticsBands_NI.v2_5.RRTMG.nc
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

    def __init__ (self,aerFiles,varMap=None,mieRootDir=None,band=False,verbose=False):
        """
        Lazy loads GOCART missing ration files and corresponding Mie tables.

        aerFiles:  str, list, or Dataset with aerosol tracers
        varMap:    str or YAML file mapping GOCART variables to specific
                   bins in Mie Tables. If None, uses internal default.
        mieRootDir: str, prepend string to mieTable file names.
        band:       bool, by default monochromatic tables are loaded.
                    If band is True, tables for radiation will be loaded instead.

        """

        if varMap is None:
            varMap = G2G_MieMap
            
        self.verbose = verbose
        self.ds = xc.open_mfdataset(aerFiles)
        self.mt = yaml.safe_load(varMap)

        if mieRootDir is None:
            edir = ''
        else:
            edir = mieRootDir + '/'
        
        for s in self.mt:
            m = self.mt[s]
            if band:
                m['mie'] = mt.MIETABLE(edir+m['bandFile'])
            else:
                m['mie'] = mt.MIETABLE(edir+m['monoFile'])

    
        (p, m, b, r, c) = self.mt['DU']['mie'].getDims() # dimensions of Mie Tables
        self.mtdims =  (p, m, b, r, c)
        self.vector = True
        for s in self.mt:
           mtdims_ = self.mt[s]['mie'].getDims() # dimensions of Mie Tables
           if self.mtdims != mtdims_:
               print('Warning: inconsistent Mie Table dimensions, will not be able to calculate PMOM')
               print('- DU (p,m,b,r,c):',self.mtdims)
               print('- '+s+' (p,m,b,r,c):',mtdims_)
               self.vector = False
               break

        
    def getAOPrt(self,Species=None,wavelength=None,vector=False):
        """
        Returns an xarray Dataset with (aot,ssa,g) if vector is
        False, otherwise (aot,ssa,g,pmon) if vector is True.

        Species:  None, str, or list. If None, all species on file,
                  otherwise subset of species.
                  
        Wavelength: float, wavelength in nm. 

        """


        if vector and not self.vector:
            print('Warning: will not calculate PMOM because of inconsistent Mie Tables.')
            vector = False
        
        # All species on file or a subset
        # -------------------------------
        if Species is None:
            Species = list(self.mt.keys())
        if isinstance(Species,str):
            Species = [Species,]

        a = self.ds    # aerosol mixing ratio tracers

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
        
        aot, ssa, g = np.zeros(space), np.zeros(space), np.zeros(space)
        if vector:
            pmom = np.zeros(space+(p,m))
        for s in Species:   # species

            if self.verbose:
                print('[] working on',s)
            
            Tracers = self.mt[s]['tracers']
            mie = self.mt[s]['mie']

            bin = 1
            for q in Tracers:

                if self.verbose:
                    print('   -',q)

                
                q_mass = rhodz * a[q]
                aot_, ssa_, g_ = mie.getAOPscalar(bin, rh, q_mass, wavelength)

                aot_ = aot_.values.squeeze()
                ssa_ = ssa_.values.squeeze()
                g_ = g_.values.squeeze()
                
                aot += aot_
                scat_ = ssa_ * aot_  # scattering AOT
                ssa += scat_         # accumulate scattering AOT here
                g += g_ * scat_

                if vector:
                    pmom_ = mie.getAOP('pmom', bin, rh,
                                       q_mass=q_mass, wavelength=wavelength)
                    pmom_ = pmom_.values.squeeze()
                    print('pmom',pmom.shape)
                    print('pmom_',pmom_.shape)
                    pmom += pmom_ * scat_.reshape(space+(1,1))  
                    
                bin += 1
                
        # Final normalization of SSA and g
        # --------------------------------
        g = g / ssa       # ssa here as TOTAL scattering AOT
        if vector:
             pmom = pmom / ssa.reshape(space+(1,1))
        ssa = ssa / aot   # ssa is now single scattering albedo.

        # Pack results into a Dataset
        # ---------------------------
        DA = dict( AOT = xr.DataArray(aot,dims=rh.dims,coords=rh.coords),
                    SSA = xr.DataArray(ssa,dims=rh.dims,coords=rh.coords),
                    G = xr.DataArray(g,dims=rh.dims,coords=rh.coords)
                 )

        DA['DELP'] = dp
        DA['AIRDENS'] = a['AIRDENS']
        
        if vector:
            dims = space + ('p', 'm')
            DA['pmom'] = DataArray(pmom, dims=dims)
         
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

def CLI_g2g_aop():
    """
    Command line interface
    """

    # yaml.dump(rc,open('test.yml','w'))

    data = '/Users/adasilva/data/'
    g = G2GAOP(data+'/sampled/*.nc',mieRootDir=data,verbose=True)
    ds = g.getAOPrt(wavelength=550,vector=True)

    return (g, ds)

    

    

    
