"""
    Implements calculation of aerosol optical properties based on gridded GOCART
mixing ratio files (aer_Nv) and GEOSmie optics tables.

"""

__version__ = '1.1.0'

import numpy  as np
import xarray as xr
import yaml

from . import mietable  as mt
from . import xrctl     as xc

from .constants import MAPL_GRAV as GRAV
from .constants import MAPL_AVOGAD
from .constants import MAPL_RUNIV
# Default YAML file mapping GOCART tracers in aer_Nv and the optics files
# -----------------------------------------------------------------------
G2G_MieMap = """
#
# GEOS Aerosol Mie table Definition for each of species.
# The order of the tracers and rhod correspond to the bins in the optics netcdf files.
#
#  rhod: particle density in kg m-3
#  shapefactor: factor that accounts for aerodynamic resistance of non-spherical particles
#               this is used to calculate the aerodynamic radius for PM calculations when the aerodymic flag is turned on
#               see the following reference for further documentation
#               GMAO Office Note No. 22 (Version 1.1): 
#               Collow, A., V. Buchard, M. Chin, P. Colarco, A. Darmenov, and A. da Silva, 2023. 
#               Supplemental Documentation for GEOS Aerosol Products
#  pmconversion: additional factor for unaccounted aerosol species. was implemented to allow for sulfate to represent missing ammonium in MERRA-2.
#              pmconversion = 1.3756 for SU for MERRA-2, otherwise = 1

DU:
  monoFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/optics_DU.v15_3.nc4
  bandFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/opticsBands_DU.v15_3.RRTMG.nc4
  tracers:
    - DU001
    - DU002
    - DU003
    - DU004
    - DU005
  shapefactor: 1.4
  rhod:
    - 2500
    - 2650
    - 2650
    - 2650
    - 2650
  pmconversion: 1


SS:
  monoFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/optics_SS.v3_3.nc4
  bandFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/opticsBands_SS.v3_3.RRTMG.nc4
  tracers:
    - SS001
    - SS002
    - SS003
    - SS004
    - SS005
  shapefactor: 1
  rhod:
    - 2200
    - 2200
    - 2200
    - 2200
    - 2200
  pmconversion: 1


OC:
  monoFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/optics_OC.v1_3.nc4
  bandFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/opticsBands_OC.v1_3.RRTMG.nc4
  tracers:
    - OCPHOBIC
    - OCPHILIC
  shapefactor: 1
  rhod:
    - 1800
    - 1800
  pmconversion: 1

BC:
  monoFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/optics_BC.v1_3.nc4
  bandFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/opticsBands_BC.v1_3.RRTMG.nc4
  tracers:
    - BCPHOBIC
    - BCPHILIC
  shapefactor: 1
  rhod:
    - 1800
    - 1800
  pmconversion: 1

BR:
  monoFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/optics_BRC.v1_5.nc4
  bandFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/opticsBands_BRC.v1_5.RRTMG.nc4
  tracers:
    - BRCPHOBIC
    - BRCPHILIC
  shapefactor: 1
  rhod:
    - 1800
    - 1800
  pmconversion: 1

SU:
  monoFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/optics_SU.v1_3.nc4
  bandFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/opticsBands_SU.v1_3.RRTMG.nc4
  tracers:
    - SO4
  shapefactor: 1
  rhod:
    - 1700
  pmconversion: 1 

NI:
  monoFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/optics_NI.v2_5.nc4
  bandFile: ExtData/chemistry/AerosolOptics/v1.0.0/x/opticsBands_NI.v2_5.RRTMG.nc4
  tracers:
    - NO3AN1
    - NO3AN2
    - NO3AN3
  shapefactor: 1
  rhod:
    - 1725
    - 2200
    - 2650
  pmconversion: 1

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
        elif type(config) is str:
            # get a file handle
            config = open(config)

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
        # start with last species read
        dims =  self.mieTable[s]['mie'].getDims()
        self.vector = True
        self.p, self.m = (0,0)
        # loop through species and compare dims
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
        
    def getAOPrt(self,Species=None,wavelength=None,vector=False,fixrh=None,m=None):

        """
        Returns an xarray Dataset with (aot,ssa,g) if vector is
        False, otherwise (aot,ssa,g,pmon) if vector is True.

        Species:  None, str, or list. If None, all species on file,
                  otherwise subset of species.

        Wavelength: float, wavelength in nm.

        vector:     bool, whether to return full phase matrix or
                    asymmetry parameter.

        m: number of pmom moments to read in. If None, reads all on file.
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

        # pre-load RH, AIRDENS, and DELP so you don't hit dask
        # repeatedly looping through  AOP calculations
        # -------------------------------------------------------
        a['DELP'].load()
        a['AIRDENS'].load()
        a['RH'].load()

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


        # Check FIXRH option
        # --------------------------
        if fixrh is not None:
            fixrh = float(fixrh)
            if (fixrh<0.0) or (fixrh>1.0):
                raise ValueError("Your fixrh is {}, it must be between 0 - 1".format(fixrh))
            else:
                rh[:] = fixrh

        # Relevant dimensions
        # -------------------
        space = rh.shape
        aot, sca, g = np.zeros(space), np.zeros(space), np.zeros(space)
        if vector:
            ns = np.prod(space)
            p = self.p
            if m is None:
                m = self.m

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
                                        wavelength=wavelength,m=m)
                    p_, m_ = pmom_.shape[-2:]
                    pmom_ = pmom_.values.reshape((ns,p_,m_)) * sca_.reshape((ns,1,1))

                    pmom[:,:,:m_] += pmom_[:,:,:] # If species have fewer moments, pad wih zeros
                else:

                    g   += sca_ * g_


                bin += 1

        # Final normalization of SSA and g
        # protect against divide by zero
        # this can happen if you ask for the AOP of an individual species
        # and its' concentration in a layer is zero        
        # --------------------------------
        ssa = np.empty(space)
        ssa[:] = np.nan
        I = np.where(aot != 0.0)
        ssa[I] = sca[I] / aot[I]

        if vector:
             I = np.where(sca.reshape(ns) != 0.0)[0]
             pmom[I,:,:] = pmom[I,:,:] / sca.reshape((ns,1,1))[I,:,:]
             I = np.where(sca.reshape(ns) == 0.0)[0]
             pmom[I,:,:] = np.nan
             pmom = pmom.reshape(space+(p,m))
        else:
             I = np.where(sca != 0.0)
             g[I] = g[I] / sca[I]
             I = np.where(sca == 0.0)
             g[I] = np.nan


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

     
    def getAOPext(self,Species=None,wavelength=None,fixrh=None,doaback=True):

        """
        Returns an xarray Dataset with the following variables:

        EXT:     aerosol extinction profile
        SCA:     aerosol scattering profile
        BSC:     aerosol backscatter profile
        DEPOL:   aerosol depolarization ratio
        ABACKTOA: total attenuated backscatter from the TOA
        ABACKSFC: total attenuated backscatter from the surface
        On inout,

        Species:    None, str, or list. If None, all species on file,
                    otherwise subset of species.

        Wavelength: float, wavelength in nm.

        fixrh: value between 0 and 1 representing realtive humidity

 	doaback: flag to turn on/off attenuated backscatter calculation

        """

        # All species on file or a subset
        # -------------------------------
        if Species is None:
            Species = list(self.mieTable.keys())
        if isinstance(Species,str):
            Species = [Species,]

        a = self.aer    # aerosol mixing ratio tracers

        # pre-load RH, AIRDENS, T and DELP so you don't hit dask
        # repeatedly looping through AOP calculations
        # -------------------------------------------------------
        #go back and see if this can be loaded in one line

        try:
            dp = a['DELP'].load()
        except:
            dp = a['delp'].load()
        airdens = a['AIRDENS'].load()
        rh = a['RH'].load()
        
        # Check FIXRH option
        # --------------------------
        if fixrh is not None:
            fixrh = float(fixrh)
            if (fixrh<0.0) or (fixrh>1.0):
                raise ValueError("Your fixrh is {}, it must be between 0 - 1".format(fixrh))
            else:
                rh[:] = fixrh

        # Relevant dimensions
        # -------------------
        space = rh.shape
        ext, sca, bsc, depol1, depol2=  (np.zeros(space), np.zeros(space),
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
                pback11_ = mie.getAOP('pback11',  bin, rh, wavelength=wavelength).values
                pback22_ = mie.getAOP('pback22',  bin, rh, wavelength=wavelength).values

                ext_ = ext_ * q_conc
                sca_ = sca_ * q_conc
                bsc_ = bsc_ * q_conc

                ext += ext_
                sca += sca_
                bsc += bsc_
                depol1 += (pback11_-pback22_) * sca_
                depol2 += (pback11_+pback22_) * sca_

                bin += 1

        if doaback:
            # Compute Molecular Scattering and Total Attenuated Backscatter Coefficient
            # following the methodology begining on page 147 of
            # http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19960051003.pdf
            # -----------------------------------------
            T = a['T'].load()
            delz  = dp / (GRAV * airdens)
            abackTOA, abackSFC = self.calcABACK(wavelength,T,rh,dp,airdens,delz,ext,bsc)

            abackTOA *= 1000. # km-1 sr-1
            abackSFC *= 1000. # km-1 sr-1


        # Final normalization
        # -------------------
        ext *= 1000. # m-1 to km-1
        sca *= 1000. # m-1 to km-1
        bsc *= 1000. # m-1 to km-1

        # protect against divide by zero
        # this can happen if you ask for the AOP of an individual species
        # and its' concentration in a layer is zero
        # -----------------------------------------
        depol = np.empty(space)
        depol[:] = np.nan
        I = np.where(depol2 != 0.0)
        depol[I] = depol1[I] / depol2[I]

        # Attributes
        # ----------
        A = dict (EXT = {'long_name':'Aerosol Extinction Coefficient', 'units':'km-1'},
                  SCA = {'long_name':'Aerosol Scattering Coefficient', 'units':'km-1'},
                  BSC = {'long_name':'Aerosol Backscatter Coefficient', 'units':'km-1'},
                  DEPOL = {'long_name':'Depolarization Ratio', 'units':'1'}
                  )
        if doaback:
            A['ABACKTOA'] = {'long_name':'Total Attenuated Backscatter Coefficient from TOA','units':'km-1 sr-1'}
            A['ABACKSFC'] = {'long_name':'Total Attenuated Backscatter Coefficient from Surface','units':'km-1 sr-1'}

        # Pack results into a Dataset
        # ---------------------------
        DA = dict(  EXT = xr.DataArray(ext.astype('float32'),dims=rh.dims,coords=rh.coords,attrs=A['EXT']),
                    SCA = xr.DataArray(sca.astype('float32'),dims=rh.dims,coords=rh.coords,attrs=A['SCA']),
                    BSC = xr.DataArray(bsc.astype('float32'),dims=rh.dims,coords=rh.coords,attrs=A['BSC']),
                    DEPOL = xr.DataArray(depol.astype('float32'),dims=rh.dims,coords=rh.coords,attrs=A['DEPOL'])
                 )

        if doaback:
            DA['ABACKTOA'] = xr.DataArray(abackTOA.astype('float32'),dims=rh.dims,coords=rh.coords,attrs=A['ABACKTOA'])
            DA['ABACKSFC'] = xr.DataArray(abackSFC.astype('float32'),dims=rh.dims,coords=rh.coords,attrs=A['ABACKSFC'])

        DA['DELP'] = dp
        DA['AIRDENS'] = a['AIRDENS']

        return xr.Dataset(DA)
        
    def calcABACK(self,wavelength,T,rh,dp,airdens,delz,ext,bsc):
        # Compute Molecular Scattering and Total Attenuated Backscatter Coefficient
        # following the methodology begining on page 147 of
        # http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19960051003.pdf
        # -----------------------------------------
        space = rh.shape
        ndims = len(space)
        km = space[1]

        abackTOA, abackSFC, pressure = (np.zeros(space), np.zeros(space), np.zeros(space))
  
        # get the pressure vertical profile
        if ndims==2:
            pe=np.zeros((space[0],km+1))
            pe[:,0]=1 #assume TOA has a pressure of 1 Pa
            for k in range(0,km):
                pe[:,k+1]=pe[:,k]+dp[:,k]
            for k in range(km):
                pressure[:,k]=(pe[:,k]+pe[:,k+1])/2
        elif ndims==4:
            pe=np.zeros((space[0],km+1,space[2],space[3]))
            pe[:,0,:,:]=1 #assume TOA has a pressure of 1 Pa
            for k in range(0,km):
                pe[:,k+1,:,:]=pe[:,k,:,:]+dp[:,k,:,:]
            for k in range(km):
                pressure[:,k,:,:]=(pe[:,k,:,:]+pe[:,k+1,:,:])/2

        # calcualte molecular backscatter
        backscat_mol = (5.45e-32/(MAPL_RUNIV/MAPL_AVOGAD)) * (wavelength/550)**-4  * pressure / T #molecular backscatter coefficient in m-1

        # calculate backscatter
        tau_mol_layer = backscat_mol * 8 * np.pi /3 * delz
        tau_aer_layer = ext * delz
        if ndims==2:
            ###TOA
            abackTOA[:,0]=(bsc[:,0]+ backscat_mol[:,0]) * np.exp(-tau_aer_layer[:,0]) * np.exp(-tau_mol_layer[:,0])
            for k in range(1,km):
                tau_aer=0
                tau_mol=0
                for kk in range(0,k):
                    tau_aer += tau_aer_layer[:,kk]
                    tau_mol += tau_mol_layer[:,kk]
                tau_aer += 0.5 * tau_aer_layer[:,k]
                tau_mol += 0.5 * tau_mol_layer[:,k]
                abackTOA[:,k] = (bsc[:,k] + backscat_mol[:,k]) * np.exp(-2*tau_aer) * np.exp(-2*tau_mol)

            ###Surface
            abackSFC[:,0]=(bsc[:,km-1]+ backscat_mol[:,km-1]) * np.exp(-tau_aer_layer[:,km-1]) * np.exp(-tau_mol_layer[:,km-1])
            for k in range(km-2,-1,-1):
                tau_aer=0
                tau_mol=0
                for kk in range(km-1,k-1,-1):
                    tau_aer += tau_aer_layer[:,kk]
                    tau_mol += tau_mol_layer[:,kk]
                tau_aer += 0.5 *  tau_aer_layer[:,k]
                tau_mol += 0.5 *  tau_mol_layer[:,k]
                abackSFC[:,k] = (bsc[:,k] + backscat_mol[:,k]) * np.exp(-2*tau_aer) * np.exp(-2*tau_mol)

        if ndims==4:
            ###TOA
            abackTOA[:,0,:,:]=(bsc[:,0,:,:]+ backscat_mol[:,0,:,:]) * np.exp(-tau_aer_layer[:,0,:,:]) * np.exp(-tau_mol_layer[:,0,:,:])
            for k in range(1,km):
                tau_aer=0
                tau_mol=0
                for kk in range(0,k):
                    tau_aer += tau_aer_layer[:,kk,:,:]
                    tau_mol += tau_mol_layer[:,kk,:,:]
                tau_aer += 0.5 *  tau_aer_layer[:,k,:,:]
                tau_mol += 0.5 *  tau_mol_layer[:,k,:,:]
                abackTOA[:,k,:,:] = (bsc[:,k,:,:] + backscat_mol[:,k,:,:]) * np.exp(-2*tau_aer) * np.exp(-2*tau_mol)

            ###Surface
            abackSFC[:,0,:,:]=(bsc[:,km-1,:,:]+ backscat_mol[:,km-1,:,:]) * np.exp(-tau_aer_layer[:,km-1,:,:]) * np.exp(-tau_mol_layer[:,km-1,:,:])
            for k in range(km-2,-1,-1):
                tau_aer=0
                tau_mol=0
                for kk in range(km-1,k-1,-1):
                    tau_aer += tau_aer_layer[:,kk,:,:]
                    tau_mol += tau_mol_layer[:,kk,:,:]
                tau_aer += 0.5 *  tau_aer_layer[:,k,:,:]
                tau_mol += 0.5 *  tau_mol_layer[:,k,:,:]
                abackSFC[:,k,:,:] = (bsc[:,k,:,:] + backscat_mol[:,k,:,:]) * np.exp(-2*tau_aer) * np.exp(-2*tau_mol)

        return abackTOA, abackSFC
 
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

    def getPM(self,Species=None,pmsize=None,fixrh=None,aerodynamic=False):
        """
        Returns an xarray Dataset with total aerosol mass smaller than the prescribed size.

        Species:  None, str, or list. If None, all species on file,
                  otherwise subset of emissions.
    
        PMsize: float, particle diameter threshold in microns. If None, the total PM is calculated.

	Please see m2_pm25.yaml and g2g_pm25.yaml for example yaml configurations.

        """

        # All species on file or a subset
        # -------------------------------
        if Species is None:
            Species = list(self.mieTable.keys())
        if isinstance(Species,str):
            Species = [Species,]

        a = self.aer    # aerosol mixing ratio tracers

        # pre-load RH, AIRDENS, and DELP so you don't hit dask
        # repeatedly looping through  AOP calculations
        # -------------------------------------------------------
        a['DELP'].load()
        a['AIRDENS'].load()
        a['RH'].load()

        # Determine PM Threshold
        # -------------------------------
        if pmsize is None:
            rPM = None
        else: 
            rPM = float(pmsize)/2 #convert diameter to radius

        # GEOS files can be inconsistent when it comes to case
        # ----------------------------------------------------
        try:
            dp = a['DELP'] 
        except:
            dp = a['delp']

        rh = a['RH']

        # Check FIXRH option
        # --------------------------
        if fixrh is not None:
            fixrh = float(fixrh)
            if (fixrh<0.0) or (fixrh>1.0):
                raise ValueError("Your fixrh is {}, it must be between 0 - 1".format(fixrh))
            else:
                rh[:] = fixrh

        # Relevant dimensions
        # -------------------
        space = rh.shape

        pm = np.zeros(space)
        for s in Species:   # species

            if self.verbose:
                print('[] working on',s)
            
            Tracers = self.mieTable[s]['tracers']
            mie = self.mieTable[s]['mie']

            bin = 1
            for q in Tracers:

                if self.verbose:
                    print('   -',q)


                # Aerosol mass concentration in kg/m3                
                q_conc = (a['AIRDENS'] * a[q]).values
  
                # Dry aerosol density in kg m-3
                # rhod is not in all of the standard optics files, and is for now read from the yaml config 
                # rhod_ = mie.getAOP('rhod',  bin, rh, wavelength=wavelength).values
                rhod_ = self.mieTable[s]['rhod'][bin-1] 

                # Lower and upper bound of the bin's radius converted from meters to microns
                rLow_ = mie.getBinInfo('rLow', bin)*1000000 
                rUp_ = mie.getBinInfo('rUp', bin)*1000000 

                # Effective radius at the specified humidity converted from meters to microns
                rEff_ = mie.getAOP('rEff', bin, rh, wavelength=None).values*1000000 

                # Effective radius at a relative humidity of 0% converted from meters to microns
                rEff_zero = mie.getBinInfo('rEffDry', bin)*1000000 

                # If necessary, compute the aerodynamic particle radius
                # shape factor accounts for changes in the particle's dragging coefficient (see https://doi.org/10.1029/2002JD002485 for more info)
                if aerodynamic:
                    # convert rhod from kg m-3 to g cm-3
                    rLow_ = rLow_ * np.sqrt((rhod_/1000)/self.mieTable[s]['shapefactor']) 
                    rUp_ = rUp_ * np.sqrt((rhod_/1000)/self.mieTable[s]['shapefactor']) 

                # Find fraction of bin that is below the threshhold
                if rPM is None:
                    # getting total PM
                    fPM = 1.0
                else:
                    if(rUp_ <= rPM):
                        fPM = 1.0
                    else:
                        if(rLow_ < rPM):
                                # in log space get the fraction of the radius bin range covered         	
                                fPM = np.log(rPM/rLow_) / np.log(rUp_/rLow_)
                        else:
                                fPM = 0.0

                # Compute the hygroscopic growth factor based on RH 
                # this is based on a formulation from GEOS Chem 
                # (https://wiki.seas.harvard.edu/geos-chem/index.php/Particulate_matter_in_GEOS-Chem)
                # this is not the same hygroscopic growth factor that is in the GEOSmie optics files.
                rhow = 997.0  # density of water at 25 C and 1 atm in kg m-3
                growthfactor= 1 + (((np.squeeze(rEff_) / np.squeeze(rEff_zero))**3 - 1) * (rhow / rhod_))
                #Compute PM
                pm_ = q_conc * growthfactor * fPM * self.mieTable[s]['pmconversion']
                pm += pm_

                bin += 1
                

        # convert from kg m-3 to micrograms m-3
        # a more common unit for PM concentration
        # ---------------------------------------
        pm = pm*1e9


        # Attributes
        # ----------
        A = dict (PM = {'long_name':'Particulate Matter', 'units':'microgram m-3'}
                  )
        
        # Pack results into a Dataset
        # ---------------------------
        DA = dict(  PM = xr.DataArray(pm.astype('float32'),dims=rh.dims,coords=rh.coords,attrs=A['PM'])
                 )

        DA['DELP'] = dp
        DA['AIRDENS'] = a['AIRDENS']
        
        return xr.Dataset(DA)


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
    fixrh = None
    d_pm = None
    aerodynamic = False

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] aerDataset [iso_t1 iso_t2]\n"+\
                                "   aerDataset          GrADS-style ctl or a shell-style wildcard string\n"+\
                                "                       with aerosol mixing ratios, either gridded or sampled"+\
                                "   iso_t1,iso_t2       optional beginning and ending time (ISO format)",
                          version=__version__ )

    parser.add_option("-a", "--aop", dest="aop", default='ext',

              help="AOP collection, one of 'rt' or 'ext' or 'pm' (default=%s)"%aop)
    

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

    parser.add_option("--rh", dest="fixrh", default=fixrh,
              help="If specified use provided RH (0-1) in calculations")


    parser.add_option("--aerodynamic", dest="aerodynamic", default=False,
              help="If set to true, an aerodynamic diameter will be used to compute PM. This option is only valid for --aop=pm.")

    parser.add_option("-s","--size", dest="d_pm", default=None,
              help="The threshold diameter size used to compute PM in units of microns (example 2.5 for PM2.5). This option is only valid for --aop=pm.")

    parser.add_option("--noaback", dest="noaback", action="store_true",
              help="Do not calculate aerosol backscatter when --aop=ext. Note, aerosol backscatter requires a temperature vertical profile.")



    (options, args) = parser.parse_args()

    # store doaback flag
    doaback = True
    if args.noaback:
        doaback = False

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
    aer = xc.open_mfdataset(aerDataset,parallel=True,chunks='auto',engine='netcdf4')
    g = G2GAOP(aer,config=config,mieRootDir=options.rootDir,verbose=options.verbose)
    for w_ in options.wavelengths.split(','):
        w = float(w_)
        if options.aop == 'ext':
            ds = g.getAOPext(wavelength=w,fixrh=options.fixrh,doaback=doaback)
        elif options.aop == 'rt':
            ds = g.getAOPrt(wavelength=w,vector=options.vector,fixrh=options.fixrh)
        elif options.aop == 'pm':
            ds = g.getPM(pmsize=options.d_pm,fixrh=options.fixrh,aerodynamic=options.aerodynamic)
        else:
            print(options.aop)
            raise AOPError('Unknown AOP option '+options.aop)

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


    data = '/discover/nobackup/acollow/aeroeval/opticsfiles/AerosolOptics/'
    #aer_Nv = '/Users/adasilva/data/sampled/aer_Nv/CAMP2Ex-GEOS-MODISonly-aer-Nv-P3B_Model_*.nc'
    aer_Nv = '/discover/nobackup/acollow/CAMP2Ex/sampled/P3B/MODISonly/2019-09-*/CAMP2Ex-GEOS-MODISonly-aer-Nv-P3B_Model_201909*_R0.nc'

    #data = '/Users/adasilva/data/'
    #aer_Nv = '/Users/adasilva/data/sampled/aer_Nv/CAMP2Ex-GEOS-MODISonly-aer-Nv-P3B_Model_*.nc'
    #aer = xr.open_mfdataset(aer_Nv) # still having trouble with parallel


    aer = xr.open_mfdataset(aer_Nv) # still having trouble with parallel
    print('aer is loaded')
    g = G2GAOP(aer,mieRootDir=data,verbose=True)
    print('g is done')
    #rts = None # g.getAOPrt(wavelength=550,vector=False)
    #rtv = g.getAOPrt(wavelength=550,vector=True)
    #ext = None # g.getAOPext(wavelength=550)
    pm = g.getPM(pmsize=2.5)
    return (pm)

if __name__ == "__main__":

    pm = Test_g2g_aop








