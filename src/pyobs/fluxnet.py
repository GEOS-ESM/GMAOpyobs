"""
Reads Level 3 grid MXD43 BRDF files.

"""
import numpy  as np
import xarray as xr
import pandas as pd

#...........................................................................

class FLUXNET(object):

    """
    Simple class for handling Joanna's Flux Net files.
    """

    def __init__(self,filen):

        self.ds = xr.open_dataset(filen)

        Name = self.ds['Name']

    def stations(self):
        """
        Return unique list of stations and coordinates as a DataFrame.
        """
        Lon  = self.ds['Longitude']
        Lat  = self.ds['Latitude']
        Name = self.ds['Name']
        
        Stations = np.unique(Name)

        lons, lats = [], [] # station coordinates
        for stn in Stations:

            I = Name==stn

            lons += [Lon[I][0],]
            lats += [Lat[I][0],]

        coords = dict( lons=np.array(lons), lats=np.array(lats) )
        return pd.DataFrame(coords,index=Stations)     
    

if __name__ == "__main__":

      fluxnet_fn = '/Users/adasilva/data/brdf/onefluxnet_daily_mcd43_c61.nc'

      stations = FLUXNET(fluxnet_fn).stations()

      print(stations)

      stations.to_csv('fluxnet_stations.csv')

      



