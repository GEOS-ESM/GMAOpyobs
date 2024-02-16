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

    def __init__ (self,aerFiles,varMap=None,mieRootDir=None,band=False,verbose=False):
        """
        Lazy loads GOCART mixing ratio *=(aer_NV) files and corresponding Mie tables.

        aerFiles:  str, list, or Dataset with aerosol tracers
        varMap:    str or YAML file handle mapping GOCART variables to specific
                   bins in Mie Tables. If None, uses internal default.
        mieRootDir: str, prepend string to mieTable file names.
        band:       bool, by default monochromatic tables are loaded.
                    If band is True, tables for radiation will be loaded instead.

        """

        if varMap is None:
            varMap = G2G_MieMap
            
        self.verbose = verbose
        if isinstance(aerFiles,xr.Dataset):
            self.ds = aerFiles
        else:
            self.ds = xc.open_mfdataset(aerFiles)
        self.mt = yaml.safe_load(varMap)

        if mieRootDir is None:
            edir = ''
        else:
            edir = mieRootDir + '/'

        # Band or monochromatic files
        # ---------------------------
        for s in self.mt:
            m = self.mt[s]
            if band:
                m['mie'] = mt.MIETABLE(edir+m['bandFile'])
            else:
                m['mie'] = mt.MIETABLE(edir+m['monoFile'])

        # Check consistency of Mie tables accross species
        # -----------------------------------------------
        dims =  self.mt['DU']['mie'].getDims() 
        self.vector = True
        self.p, self.m = (0,0)
        for s in self.mt:
           dims_ = self.mt[s]['mie'].getDims() # dimensions of Mie Tables, a dict
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
           self.p = max(self.p,dims_['p'])
           self.m = max(self.m,dims_['m'])
        
    def getAOPrt(self,Species=None,wavelength=None,vector=False):
        """
        Returns an xarray Dataset with (aot,ssa,g) if vector is
        False, otherwise (aot,ssa,g,pmon) if vector is True.

        Species:  None, str, or list. If None, all species on file,
                  otherwise subset of species.
                  
        Wavelength: float, wavelength in nm. 

        """

        # This is not working yet
        # -----------------------
        if vector:
            if not self.vector:
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

        # Bound RH
        # --------
        rh = a['RH'].values
        rh[rh<0] = 0.0
        rh[rh>0.99] = 0.99
        rh = xr.DataArray(rh,dims=dp.dims, coords=dp.coords)

        # Relevant dimensions
        # -------------------
        space = rh.shape
        aot, ssa, g = np.zeros(space), np.zeros(space), np.zeros(space)
        if vector:
            ns = np.prod(space)
            p, m = self.p, self.m
            pmom = np.zeros(space+(p,m)).reshape((ns,p,m))
        
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
                    #breakpoint()
                    pmom_ = mie.getAOP('pmom', bin, rh,
                                       q_mass=q_mass, wavelength=wavelength)
                    #breakpoint()
                    
                    pmom_ = pmom_.values.squeeze() * scat_.reshape(space+(1,1))
                    p_, m_ = pmom_.shape[-2:]
                    
                    pmom[:,:,:m_] += pmom_[:,:,:] # If species have fewer moments, pad wih zeros
                    
                bin += 1
                
        # Final normalization of SSA and g
        # --------------------------------
        g = g / ssa       # ssa here as TOTAL scattering AOT
        if vector:
             pmom = pmom / ssa.reshape(space+(1,1))
             pmom = pmom.reshape(space+(p,m))
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
            DA['pmom'] = DataArray(pmom, dims=rh.dims+('p','m') )
         
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

        # Bound RH
        # --------
        rh = a['RH']
        rh[rh<0] = 0.0
        rh[rh>0.99] = 0.99
        rh = xr.DataArray(rh,dims=dp.dims, coords=dp.coords)

        # Relevant dimensions
        # -------------------
        space = rh.shape

        ext, sca, bsc, depol1, depol2 = (np.zeros(space), np.zeros(space),
                                         np.zeros(space), np.zeros(space),
                                         np.zeros(space))

        for s in Species:   # species

            if self.verbose:
                print('[] working on',s)
            
            Tracers = self.mt[s]['tracers']
            mie = self.mt[s]['mie']

            bin = 1
            for q in Tracers:

                if self.verbose:
                    print('   -',q)

                
                q_conc = (a['AIRDENS'] * a[q]).values
                ext_ = mie.getAOP('bext', bin, rh, wavelength=wavelength)
                sca_ = mie.getAOP('bsca', bin, rh, wavelength=wavelength)
                bsc_ = mie.getAOP('bbck', bin, rh, wavelength=wavelength)
                p11_ = mie.getAOP('p11',  bin, rh, wavelength=wavelength)
                p22_ = mie.getAOP('p22',  bin, rh, wavelength=wavelength)
                
                ext_ = ext_.values.squeeze() * q_conc
                sca_ = sca_.values.squeeze() * q_conc
                bsc_ = bsc_.values.squeeze() * q_conc
                p11_ = p11_.values.squeeze() * q_conc
                p22_ = p22_.values.squeeze() * q_conc

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

        # Pack results into a Dataset
        # ---------------------------
        DA = dict(  EXT = xr.DataArray(ext,dims=rh.dims,coords=rh.coords),
                    SCA = xr.DataArray(sca,dims=rh.dims,coords=rh.coords),
                    BSC = xr.DataArray(bsc,dims=rh.dims,coords=rh.coords),
                    DEPOL = xr.DataArray(depol,dims=rh.dims,coords=rh.coords)
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

def CLI_g2g_aop():
    """
    Command line interface
    """

    # yaml.dump(rc,open('test.yml','w'))

    data = '/Users/adasilva/data/'
    aer_Nv = '/Users/adasilva/data/sampled/CAMP2Ex-GEOS-MODISonly-aer-Nv-P3B_Model_*.nc'
    
    aer = xr.open_mfdataset(aer_Nv) # still having trouble with parallel

    g = G2GAOP(aer,mieRootDir=data,verbose=True)
    ds = g.getAOPrt(wavelength=550,vector=True)

    #g = G2GAOP(data+'/sampled/*aer-Nv*.nc',mieRootDir=data,verbose=True)
    #ds = g.getAOPext(wavelength=550)

    return (g, ds)

if __name__ == "__main__":

    g, ds = CLI_g2g_aop



    

    

    
