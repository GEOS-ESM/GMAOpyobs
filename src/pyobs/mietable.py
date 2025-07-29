#!/usr/bin/env python3
#

"""

   Implements API to access Version 1.0.0 of Mie LUTs produced with the GEOSmie package.
   These Mie Tables are stored in netcdf files, with the following variables:
      
   IMPORTANT: Mie Tables v0.0.0 used at GMAO up to 2024 are not dimensioined in this order.

   Coordinate Variables:
      channel          (w) channel number
      wavelengths      (w) wavelengths        [m]
      rh               (r) RH values   [fraction]
      rLow             (b) Dry lower radius [m]
      rEffDry          (b) Dry Effective radius [m]
      rUp              (b) Dry upper radius [m]
      p                (p) Non-zero elements of phase matrix
      m                (p) Moments of phase matrix
      ang              (a) number of scattering angles [degrees]

    Data Variables:
      reff             (b,r) effective radius [m]
      bext             (b,w,r) bext values [m2 kg-1]
      bsca             (b,w,r) bsca values [m2 kg-1]
      bbck             (b,w,r) bbck values [m2 kg-1]
      g                (b,w,r) asymmetry parameter
      pback11          (b,w,r) Backscatter phase function, index 0 
      pback22          (b,w,r) Backscatter phase function, index 4
      pmom             (b,w,r,p,m) moments of phase function
      pback            (b,w,r,p)   moments of backscatter phase function
      gf               (b,r) hygroscopic growth factor
      rhop             (b,r) wet particle density [kg m-3]
      rhod             (b,r) dry particle density [kg m-3]
      vol              (b,r) wet particle volume [m3 kg-1]
      area             (b,r) wet particle cross section [m2 kg-1]
      refr             (b,w,r) real part of refractive index
      refi             (b,w,r) imaginary part of refractive index
      p11,12,22,33,34,44  (b,w,r,a) scattering phase matrix elements

      In the above the dimensions are

      w   wavelength
      r   relative humidity
      b   bin number
      p   number of nonzero elements in phase matrix
      m   number of moments of phase matrix
      a   number of scattering angles

      See GEOSmie documentation for details.

      NOTE: Files record wavelengths in meters, user specifies wavelength in
            nano-meters in this package.
      

"""

import xarray as xr
import numpy  as np

__VERSION__ = '0.9.0'

supportedAOPs = ['aot',          'ssa',     'gf',      'gasym',   'g',   'growth_factor',
                 'RefIndex',     'pmom',    'area',    'volume',  'pback11', 'pback22', 'pback',
                 'rhod',         'rhop',    'rEff',    'bbck',    'tau', 'sca',
                 'bsca',         'bext',    'refreal', 'refimag', 'pmatrix',
                 'p11', 'p12', 'p33', 'p34', 'p22', 'p44',
                 'aot_ssa_pmom',
                 'aot_ssa_gasym' ]


ALIAS = {'gf'   : 'growth_factor',
         'tau'  : 'aot',
         'gasym': 'g'}

# RH quantization
# ---------------
_nrh = 1000   # 1,000 points as in GOCART2G
_rh_max = 0.99 # cap RH at 99%
_rh, _drh = np.linspace(0,_rh_max,_nrh,retstep=True) # quantized RH 
    
def _iRH(rh):
    """
    Returns index of quantized RH. Values are clipped so that indices are
    in the range[0,_nRH-1].
    """
    return ( 0.5 + (rh / _drh ) ).astype('int').compute().\
             clip(min=0,max=_nrh-1) # needs .compute() for indexing

#..........................................................................

class MieTableError(Exception):
    """
    Defines general exception errors.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


#..........................................................................
class MIETABLE(object):

   def __init__ (self, filename):
      """
      Loads GEOSmie created Aerosol Optical Property table for a
      single species.
      
      filename:  str, Mie Table file name

      """ 
      self.filename       = filename
      self.ds             = xr.open_dataset(filename)
      self.AOPs           = self.ds.data_vars  # intensive properties
      
      wavelengths         = self.ds.coords['wavelength'].values
      bin                 = self.ds.coords['bin'].values
      self.min_wavelength = min(wavelengths)
      self.max_wavelength = max(wavelengths)
      self.nbins           = len(bin)


   #--
   def getBinInfo(self,name,bin):
       """
       Return one-dimension fields from the optics files that are a function of bin
       """   
       return self.ds.coords[name].values[bin-1]

   #--
   def getDims(self):
       """
       Return dimensions of tables as a dictionary
       """
       dims = dict(self.ds.sizes)
       if 'p' not in dims:
           dims['p'] = None
           dims['m'] = None
       return dims
   
   #--
   def _getAOP(self, name, bin, wavelength=None,m=None):
      """
      Return DataArray with optical property *name* for bin and, if needed,
      wavelength. No RH interpolation performed. This is internal,

      name:       str, name of the optical property
      bin:        int, 1-offset bin number
      wavelength: float, wavelength in nm
      
      """
      assert name in self.AOPs, name + ' is not found in the table ' + self.filename
      assert 1 <= bin and bin <= self.nbins,  "bin " + str(bin) \
               + " is out of range in the file " + self.filename

      if wavelength is not None:
          wavelength_ = wavelength / 1e9 # User species nm, tables use m

      bin_ = bin - 1
      if 'wavelength' in self.ds[name].dims:
         assert wavelength_ is not None, \
                'wavelength should be provided to get variable ' \
                + name + ' in file ' + filename
         wavelength_ = min(max(wavelength_, self.min_wavelength), self.max_wavelength)
         da = self.ds[name].isel({'bin':[bin_]}).interp({'wavelength': [wavelength_]})
      else:
         da = self.ds[name].isel({'bin':[bin_]})

      # limit number of moments
      if m is not None:
         if da.sizes['m'] > m:
             da = da.isel(m=slice(0,m))

      return da.squeeze()
         
#--
   def getAOP(self, name, bin, rh, q_mass=None, wavelength=None, m=None):
      """
      Returns DataArray with Aerosol Optical Property *name*.
      
      name:       str, name of the optical property. Consult module variable
                  supportedAOPs to see which propertirs are supported.
      bin:        int, bin number
      rh:         DataArray, relative humidity in [0,1]
      q_mass:     DataArray, aerosol column mass (Kg/m2), only needed for
                  extensive properties
      wavelength: float, wavelength in nm
      m:          integer, max number of pmom moments
      
      """
      
      assert name in supportedAOPs, "Optical Property " + name + " not supported"
      if name in ALIAS : name = ALIAS[name]

      rh = rh.clip(min=0,max=_rh_max)           # clip RH as in GOCART-2G

      if name in self.AOPs:
         aop = self._getAOP(name, bin, wavelength=wavelength,m=m)
         if len(aop.dims) > 1:
             aop = aop.interp(rh=_rh)[_iRH(rh)] # Faster RH interpolation
         else:
             aop = aop.interp(rh=rh)            # Regular linear interpolation

      elif name == 'aot' :
         assert q_mass is not None, 'aot needs q_mass as input'
         bext  = self._getAOP('bext', bin, wavelength=wavelength).interp(rh=rh)
         aop   = (bext*q_mass).rename('aot')

      elif name == 'sca' :
         assert q_mass is not None, 'aot needs q_mass as input'
         bsca  = self._getAOP('bsca', bin, wavelength=wavelength).interp(rh=rh)
         aop   = (bsca*q_mass).rename('sca')

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
            pmom = self.getAOP('pmom', bin, rh, wavelength=wavelength,m=m)
            aop  = (aot, ssa, pmom)
         elif 'gasym' in name:
            gasym = self.getAOP('g', bin, rh, wavelength=wavelength).rename('gasym')
            aop   = (aot, ssa, gasym)

      elif name == 'pback11':
         pback11 = self._getAOP('pback', bin, wavelength=wavelength)
         aop     = pback11.interp(rh=rh).isel({"p": 0}, drop=True).rename('pback11')

      elif name == 'pback22':
         pback22 = self._getAOP('pback', bin, wavelength=wavelength)
         aop = pback22.interp(rh=rh).isel({"p": 4}, drop=True).rename('pback22')

      elif name == 'pmatrix':
         p11 = self._getAOP('p11',bin,wavelength=wavelength)
         p12 = self._getAOP('p12',bin,wavelength=wavelength)
         p33 = self._getAOP('p33',bin,wavelength=wavelength)
         p34 = self._getAOP('p34',bin,wavelength=wavelength)
         p22 = self._getAOP('p22',bin,wavelength=wavelength)
         p44 = self._getAOP('p44',bin,wavelength=wavelength)
         aop = xr.concat((p11.interp(rh=rh),
                          p12.interp(rh=rh),
                          p33.interp(rh=rh),
                          p34.interp(rh=rh),
                          p22.interp(rh=rh),
                          p44.interp(rh=rh)),'p')
         newdim = rh.dims+('p','ang')
         aop = aop.transpose(*newdim)
         
      else:
          raise MieTableError('Unknown AOP '+name)

      return aop

#.....................................................................

if __name__ == "__main__":

   import numpy as np
    
   # Sample Mie Tables
   # -----------------
   #dirn   = '/discover/nobackup/projects/gmao/share/dasilva/fvInput/ExtData/chemistry/AerosolOptics/v0.0.0/x/'
   dirn   = '/Users/adasilva/data/ExtData/chemistry/AerosolOptics/v1.0.0/x/'
   Tables = [dirn + 'optics_DU.v15_5.nc4', dirn + 'optics_OC.v2_3.nc4']

   # Aerosol state (all species)
   # ---------------------------
   #aer_Nv = '/css/gmao/geos-it/products/Y2023/M02/D05/GEOS.it.asm.aer_inst_3hr_glo_C180x180x6_v72.GEOS5294.2023-02-05T1200.V01.nc4'
   aer_Nv = '/Users/adasilva/data/sampled/aer_Nv/CAMP2Ex-GEOS-MODISonly-aer-Nv-P3B_Model_*.nc'
   aer    = xr.open_mfdataset(aer_Nv)

   try:
       delp = aer['DELP']
   except:
       delp = aer['delp']
           
   q_mass = [aer['DU001'] * delp / 9.81, aer['OCPHILIC'] * delp / 9.81]
   rh     = aer['RH']

   # Sample variable names
   # ---------------------
   Vars = ['tau', 'aot', 'gasym', 'bext', 'bsca', 'ssa', 'bbck', 'rEff',
           'RefIndex', 'pmom', 'pback', 'pback11', 'pback22',
           'aot_ssa_gasym', 'aot_ssa_pmom']

   # Loop over Tables
   # ----------------
   for i in range(len(Tables)):
       species = Tables[i].split('_')[1][0:2]
       mie     = MIETABLE(Tables[i])
       AOP = dict()
       print('Doing',species)
       for v in Vars:
           print('-',v)
           AOP[species,v] = mie.getAOP(v, 1, rh, wavelength=550, q_mass=q_mass[i])


