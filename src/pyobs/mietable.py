#!/usr/bin/env python3
#

"""

   Implements API to access Mie LUTs produced with the GEOSmie package.
   These Mie Tables are stored in netcdf files, with the following variables:
      
      wavelengths      (c) wavelengths        [m]
      rh               (r) RH values   [fraction]
      reff             (b,r) effective radius [m]
      bext             (b,r,c) bext values [m2 kg-1]
      bsca             (b,r,c) bsca values [m2 kg-1]
      bbck             (b,r,c) bbck values [m2 kg-1]
      g                (b,r,c) asymmetry parameter
      p11              (b,r,c) Backscatter phase function, index 1 
      p22              (b,r,c) Backscatter phase function, index 5
      pmom             (p,m,b,r,c) moments of phase function
      pback            (p,b,r,c)   moments of backscatter phase function
      gf               (b,r) hygroscopic growth factor
      rhop             (b,r) wet particle density [kg m-3]
      rhod             (b,r) dry particle density [kg m-3]
      vol              (b,r) wet particle volume [m3 kg-1]
      area             (b,r) wet particle cross section [m2 kg-1]
      refr             (b,r,c) real part of refractive index
      refi             (b,r,c) imaginary part of refractive index

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


class MIETABLE(object):

   def __init__ (self, filename):
      """
      Loads GEOSmie created Aerosol Optical Property table for a
      single species.
      
      filename:  str, Mie Table file name

      """ 
      self.filename    = filename
      self.mieDS       = xr.open_dataset(filename)
      self.AOPs        = self.mieDS.keys()
      wavelengths      = self.mieDS.coords['lambda'].values
      radius           = self.mieDS.coords['radius'].values
      self.min_wavelength = min(wavelengths)
      self.max_wavelength = max(wavelengths)
      self.bins           = len(radius)

   #--
   def getDims(self):
       """
       Return dimensions of tables
       """
       try:
           (p, m, b, r, c) = self.mieDS.pmom.shape
       except:
           (b, r, c) = self.mieDS.qext.shape
           p, m = None, None   # older tables do not include pmom
       return (p, m, b, r, c)
   
   #--
   def _getAOP(self, name, bin, wavelength=None):
      """
      Return DataArray with optical property *name* for bin and, if needed,
      wavelength. No RH interpolation performed. This internal,

      name:       str, name of the optical property
      bin:        int, 1-offset bin number
      wavelength: float, wavelength in nm
      
      """
      assert name in self.AOPs, name + ' is not found in the table ' + self.filename
      assert 1 <= bin and bin <= self.bins,  "bin " + str(bin) + " is out of range in the file " + self.filename

      if wavelength is not None:
          wavelength_ = wavelength / 1e9 # User species nm, tables use m
      else:
          wavelength_ = wavelength / 1e9 # User species nm, tables use m
          
      bin_ = bin - 1
      if 'lambda' in self.mieDS[name].dims:
         assert wavelength_ is not None, 'wavelength should be provided to get variable ' + name + ' in file ' + filename
         wavelength_ = min(max(wavelength, self.min_wavelength), self.max_wavelength)
         da = self.mieDS[name].isel({'radius':[bin_]}).interp({'lambda': [wavelength_]})
      else:
         da = self.mieDS[name].isel({'radius':[bin_]})

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
         print('in getAOP, aop 1', aop.shape)
         aop = aop.interp(rh=rh)
         print('in getAOP, aop 2', aop.shape)

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
         aop = p11.isel({"nPol": [0]}).interp(rh=rh).rename('p11')

      elif name == 'p22':
         p22 = self._getAOP('pback', bin, wavelength=wavelength)
         aop = p22.isel({"nPol": [4]}).interp(rh=rh).rename('p22')
 
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

   # Sample Mie Tables
   # -----------------
   dirn   = '/discover/nobackup/projects/gmao/share/dasilva/fvInput/ExtData/chemistry/AerosolOptics/v0.0.0/x/'
   Tables = [dirn + 'optics_DU.v7.nc', dirn + 'optics_OC.v2_3.nc']
   

   # Aerosol state (all species)
   # ---------------------------

   aer_Nv = '/css/gmao/geos-it/products/Y2023/M02/D05/GEOS.it.asm.aer_inst_3hr_glo_C180x180x6_v72.GEOS5294.2023-02-05T1200.V01.nc4'

   aer    = xr.open_dataset(aer_Nv).variables
  
   ######
   # Dust
   #####
   # ----
   table  = Tables[0]
   mie    = MieTABLE(table)#,wavelengths)
   q_mass = aer['DU003'] * aer['DELP'] / 9.81
   rh     = aer['RH']

   varnames =['tau', 'aot', 'gasym', 'bext', 'bsca', 'ssa', 'bbck', 'rEff', 'RefIndex', 'aot_ssa_gasym']

   for var_name in varnames:
      print(var_name)
      var  = mie.getVariable(var_name, 3, rh, wavelength=550e-9, q_mass=q_mass)

   #####
   # OC
   #####
   # --
   table = Tables[1]
   wavelengths = [470, 550, 670, 870]
   mie = MieTABLE(table)
   q_mass = aer['BCPHILIC'] * aer['DELP'] / 9.81
   rh = aer['RH']

   for var_name in varnames:
      print(var_name)
      for wavelength in wavelengths:
         print(wavelength)
         var  = mie.getVariable(var_name, 1, rh, wavelength=wavelength, q_mass=q_mass)
