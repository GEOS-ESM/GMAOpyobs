"""
Reads Level 3 grid MXD43 BRDF files.

"""

import os
import numpy       as np
import cartopy.crs as ccrs
import xarray      as xr
import pandas      as pd
    
from glob     import glob

from .mcd_sinu import MCD_SINU

# Override default variable names
Alias = dict ( LC_Type1 = 'LC_IGBP',
               LC_Type2 = 'LC_UMD',
               LC_Type3 = 'LC_LAI',
               LC_Type4 = 'LC_BGC',
               LC_Type5 = 'LC_PFT',
             )

#...........................................................................
#raise DeprecationWarning("Use MCD_SINU instead.")
class MCD12(MCD_SINU):
    """
    Implements interface to MODIS Land Cover yarly Level 3 products.
    """

#............................................................................

def _main():

      # Files on discover
      # -----------------
      fluxnet_fn = '/discover/nobackup/adasilva/fluxnet/fluxnet_stations.csv'
      mcd12q1_dn = '/css/modis/Collection6.1/L3/MCD12Q1-Landcover/2023/001'

      stations = pd.read_csv(fluxnet_fn, index_col=0)

      print(stations)

      lon = stations['lons'].values
      lat = stations['lats'].values

      # BRDF object properly initialized
      # --------------------------------

      land = MCD12(mcd12q1_dn,lon,lat,Alias=Alias)

      # Sample 1 variable at a time at obs locations
      # --------------------------------------------
      for vname in land.variables:
            print('[] Interpolating',vname)
            v = land.interp(vname)

      # Alternativaly, interpolate several variables
      # Note: omit variables to interpolate all variables
      # -------------------------------------------------
      Variables = ['LC_IGBP', 'LW', 'QC']
      #V = land.interp_many(Variables,Index=stations.index)
      V = land.interp_many(Variables)
      print('\n Dataframe:')
      print(V)

      return land

if __name__ == "__main__":

    _main()




