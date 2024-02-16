#!/usr/bin/env python3
#

"""

   Implements API to access Version 1.0.0 of Mie LUTs produced with the GEOSmie package.
   These Mie Tables are stored in netcdf files, with the following variables:
      
   IMPORTANT: Mie Tables v0.0.0 used at GMAO up to 2024 are not dimensioined in this order.

   Coordinate Variables:
      channel          (c) channel number
      wavelengths      (c) wavelengths        [m]
      rh               (r) RH values   [fraction]
      rLow             (b) Dry upper radius [m]
      rEffDry          (b) Dry Effective radius [m]
	  rUp              (b) Dry lower radius [m]
      p                (p) Non-zero elements of phase matrix
      m                (p) Moments of phase matrix

    Data Variables:
      reff             (b,r) effective radius [m]
      bext             (b,c,r) bext values [m2 kg-1]
      bsca             (b,c,r) bsca values [m2 kg-1]
      bbck             (b,c,r) bbck values [m2 kg-1]
      g                (b,c,r) asymmetry parameter
      p11              (b,c,r) Backscatter phase function, index 0 
      p22              (b,c,r) Backscatter phase function, index 4
      pmom             (b,c,r,p,m) moments of phase function
      pback            (b,c,r,p)   moments of backscatter phase function
      gf               (b,r) hygroscopic growth factor
      rhop             (b,r) wet particle density [kg m-3]
      rhod             (b,r) dry particle density [kg m-3]
      vol              (b,r) wet particle volume [m3 kg-1]
      area             (b,r) wet particle cross section [m2 kg-1]
      refr             (b,c,r) real part of refractive index
      refi             (b,c,r) imaginary part of refractive index

      In the above the dimensions are

      c   channel
      r   relative humidity
      b   bin number
      p   number of nonzero elements in phase matrix
      m   number of moments of phase matrix

      See GEOSmie documentation for details.

      NOTE: Files record wavelengths in meters, user specifies wavelength in
            namo-meters in this package.
      

"""

import xarray as xr
import numpy  as np

__VERSION__ = '0.9.0'

supportedAOPs = ['aot',          'ssa',     'gf',      'gasym',   'g',   'growth_factor',
                 'RefIndex',     'pmom',    'area',    'volume',  'p11', 'p22', 'pback',
                 'rhod',         'rhop',    'rEff',    'bbck',    'tau',
                 'bsca',         'bext',    'refreal', 'refimag',
                 'aot_ssa_pmom',
                 'aot_ssa_gasym' ]


ALIAS = {'gf'   : 'growth_factor',
         'tau'  : 'aot',
         'gasym': 'g'}

# RH quanrization
# ---------------
_nrh = 10000   # 10,000 points
_rh_max = 0.99 # cap RH at 99%
_rh = np.linspace(0,_rh_max,_nrh) # quantized RH 
_drh = _rh_max / (_nrh - 1)
    
def _rhInterp(da,rh):
    """
    Fast RH interpolation.
    """
    if da.dims[0] != 'rh':
        raise ValueError('First dimension must be rh')
    q_da = da.interp(rh=_rh) # interpolate to high res LUT
    I = ( 0.5 + (rh.values.ravel() / _drh ) ).astype('int')
    da_ = da.values[I]

 
class MIETABLE(object):

   def __init__ (self, filename):
      """
      Loads GEOSmie created Aerosol Optical Property table for a
      single species.
      
      filename:  str, Mie Table file name

      """ 
      self.filename       = filename
      self.ds          = xr.open_dataset(filename)
      self.AOPs           = self.ds.data_vars  # intensive properties
      wavelengths         = self.ds.coords['wavelength'].values
      bin                 = self.ds.coords['bin'].values
      self.min_wavelength = min(wavelengths)
      self.max_wavelength = max(wavelengths)
      self.bins           = len(bin)

   #--
   def getDims(self):
       """
       Return dimensions of tables as a dictionary
       """
       dims = dict(self.ds.dims)
       if 'p' not in dims:
           dims['p'] = None
           dims['m'] = None
       return dims
   
   #--
   def _getAOP(self, name, bin, wavelength=None):
      """
      Return DataArray with optical property *name* for bin and, if needed,
      wavelength. No RH interpolation performed. This is internal,

      name:       str, name of the optical property
      bin:        int, 1-offset bin number
      wavelength: float, wavelength in nm
      
      """
      assert name in self.AOPs, name + ' is not found in the table ' + self.filename
      assert 1 <= bin and bin <= self.bins,  "bin " + str(bin) + " is out of range in the file " + self.filename

      if wavelength is not None:
          wavelength_ = wavelength / 1e9 # User species nm, tables use m
          
      bin_ = bin - 1
      if 'wavelength' in self.ds[name].dims:
         assert wavelength_ is not None, 'wavelength should be provided to get variable ' + name + ' in file ' + filename
         wavelength_ = min(max(wavelength, self.min_wavelength), self.max_wavelength)
         da = self.ds[name].isel({'bin':[bin_]}).interp({'wavelength': [wavelength_]})
      else:
         da = self.ds[name].isel({'bin':[bin_]})

      return da # data array with optical property for bin and 
  
#--
   def getAOP(self, name, bin, rh, q_mass=None, wavelength=None):
      """
      Returns DataArray with Aerosol Optical Property *name*.
      
      name:       str, name of the optical property. Consult module variable
                  supportedAOPs to see which propertirs are supported.
      bin:        int, bin number
      q_mass:     DataArray, aerosol column mass (Kg/m2), only needed for
                  extensive properties
      wavelength: float, wavelength in nm
      
      """
      
      assert name in supportedAOPs, "Optical Property " + name + "not supported"
      if name in ALIAS : name = ALIAS[name]
      
      if name in self.AOPs:
         aop = self._getAOP(name, bin, wavelength=wavelength)
         aop = aop.interp(rh=rh)

      elif name == 'aot' :
         assert q_mass is not None, 'aot needs q_mass as input'
         bext  = self._getAOP('bext', bin, wavelength=wavelength)
         bext_ = bext.interp(rh=rh)
         aop   = (bext_*q_mass).rename('aot')

      elif name == "ssa":
         bext = self._getAOP('bext', bin, wavelength=wavelength)           
         bsca = self._getAOP('bsca', bin, wavelength=wavelength)           
         ssa  = bsca/bext
         aop  = ssa.interp(rh=rh).rename('ssa')

      elif name == 'volume':
         rhod = self._getAOP('rhod', bin)
         gf   = self._getAOP('gf',   bin)
         vol  = gf**3/rhod
         aop  = vol.interp(rh=rh).rename('volume')

      elif name == 'area':
         rhod  = self._getAOP('rhod', bin)
         gf    = self._getAOP('gf',   bin)
         reff  = self._getAOP('rEff', bin)
         vol   = gf**3/rhod
         area  = vol/(4./3.*reff)
         aop   = area.interp(rh=rh).rename('area')

      elif name == 'RefIndex':
         refr = self._getAOP('refreal', bin, wavelength=wavelength)
         refi = self._getAOP('refimag', bin, wavelength=wavelength)
         aop  = (refr.interp(rh=rh), refi.interp(rh=rh))

      elif name == 'aot_ssa_pmom' or name == 'aot_ssa_gasym':
         assert q_mass is not None, name + 'needs q_mass as input'
         bext = self._getAOP('bext', bin, wavelength=wavelength)
         bsca = self._getAOP('bsca', bin, wavelength=wavelength)
         ssa  = (bsca/bext).interp(rh=rh).rename('ssa')
         aot  = (bext.interp(rh=rh) * q_mass).rename('aot')
         if 'pmom' in name:
            pmom = self._getAOP('pmom', bin, wavelength=wavelength)
            pmom = pmom.interp(rh=rh)
            aop  = (aot, ssa, pmom)
         elif 'gasym' in name:
            gasym = self._getAOP('g', bin, wavelength=wavelength)
            gasym = gasym.interp(rh=rh)
            aop   = (aot, ssa, gasym)

      elif name == 'p11':
         p11 = self._getAOP('pback', bin, wavelength=wavelength)
         aop = p11.isel({"p": [0]}).interp(rh=rh).rename('p11')

      elif name == 'p22':
         p22 = self._getAOP('pback', bin, wavelength=wavelength)
         aop = p22.isel({"p": [4]}).interp(rh=rh).rename('p22')
 
      return aop

   def getAOPscalar(self, bin, rh, q_mass, wavelength):
      """
      Returns tuple with 3 DataArrays with (aot,ssa,g) useful for
      scalar radiative transfer calculations.
      
      name:       str, name of the optical property
      bin:        int, bin number
      q_mass:     DataArray, aerosol column mass (Kg/m2)
      wavelength: float, wavelength in nm
      
      """

      return self.getAOP('aot_ssa_gasym', bin, rh,
                          q_mass=q_mass, wavelength=wavelength)

   def getAOPvector(self, bin, rh, q_mass, wavelength):
      """
      Returns tuple with 3 DataArrays with (aot,ssa,p) useful for
      vector radiative transfer calculations.
      
      name:       str, name of the optical property
      bin:        int, bin number
      q_mass:     DataArray, aerosol column mass (Kg/m2)
      wavelength: float, wavelength in nm
      
      """

      return self.getAOP( 'aot_ssa_pmom', bin, rh,
                          q_mass=q_mass, wavelength=wavelength)

#.....................................................................

if __name__ == "__main__":

   import numpy as np
    
   # Sample Mie Tables
   # -----------------
   #dirn   = '/discover/nobackup/projects/gmao/share/dasilva/fvInput/ExtData/chemistry/AerosolOptics/v0.0.0/x/'
   dirn   = '/Users/adasilva/data/ExtData/chemistry/AerosolOptics/v1.0.0/x/'
   Tables = [dirn + 'optics_DU.v7.nc4', dirn + 'optics_OC.v2_3.nc4']

   # Aerosol state (all species)
   # ---------------------------

   #aer_Nv = '/css/gmao/geos-it/products/Y2023/M02/D05/GEOS.it.asm.aer_inst_3hr_glo_C180x180x6_v72.GEOS5294.2023-02-05T1200.V01.nc4'
   #aer_Nv = '/Users/adasilva/data/sampled/CAMP2Ex-GEOS-MODISonly-aer-Nv-P3B_Model_20190904_R0.nc'
   aer_Nv = '/Users/adasilva/data/sampled/CAMP2Ex-GEOS-MODISonly-aer-Nv-P3B_Model_*.nc'
   aer    = xr.open_mfdataset(aer_Nv)
  
   ######
   # Dust
   #####
   # ----
   table  = Tables[0]
   mie    = MIETABLE(table)#,wavelengths)
   q_mass = aer['DU003'] * aer['delp'] / 9.81
   rh     = aer['RH']

   q_rh = np.linspace(0,0.99,1000)
   q_pmom = mie.getAOP('pmom', 3, q_rh, wavelength=550)

def xxx():
   
   varnames =['tau', 'aot', 'gasym', 'bext', 'bsca', 'ssa', 'bbck', 'rEff', 'RefIndex', 'aot_ssa_gasym']
   varnames =['pmom',]

   DU = dict()
   for var_name in varnames:
      print(var_name)
      DU[var_name]  = mie.getAOP(var_name, 3, rh, wavelength=550, q_mass=q_mass)


   #####
   # OC
   #####
   # --
   table = Tables[1]
   wavelengths = [470, 550, 670, 870]
   mie = MIETABLE(table)
   q_mass = aer['BCPHILIC'] * aer['delp'] / 9.81
   rh = aer['RH']

   OC = dict()
   for var_name in varnames:
      print(var_name)
      for wavelength in wavelengths:
         print(wavelength)
         OC[var_name]  = mie.getAOP(var_name, 1, rh, wavelength=wavelength, q_mass=q_mass)

   # Multiples
   # ---------
   aot_, ssa_, g_ = mie.getAOPscalar(1, rh, q_mass, 550)

