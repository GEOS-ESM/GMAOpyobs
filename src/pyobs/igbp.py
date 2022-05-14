"""
A Python interface to CPTEC's IGBP interface.
"""

import warnings
warnings.simplefilter('ignore', DeprecationWarning)
warnings.simplefilter('always', UserWarning)


# Vegetation types
# ----------------
TROPICAL       = 1
EXTRA_TROPICAL = 2
SAVANNA        = 3
GRASSLAND      = 4 

import sys
from numpy    import any

import IGBP_  # load f2py extension

#                      tf    xf   sv    gl
FLAMING_FRACTION = [ 0.45, 0.45, 0.75, 0.97 ]

#................................ Static Methods ...............................

class IGBP(object):
    """
    Base class to add IGBP capabilities to other obs classes.
    """

    def simpleVeg(self,Path,nonVeg=GRASSLAND):
        self.veg = getSimpleVeg(self.lon,self.lat,Path=Path,nonVeg=nonVeg)

    def detailedVeg(self,Path): 
        self.veg = getDetailedVeg(self.lon,self.lat,Path=Path)

def getSimpleVeg(lon,lat,Path,nonVeg=GRASSLAND): 
    """
        Given the pathname for the IGBP datasets, returns
    information about the *aggregated* vegetation type:  

         1  Tropical Forests        IGBP  2, 30S < lat < 30N
         2  Extra-tropical Forests  IGBP  1, 2(lat <=30S or lat >= 30N), 
                                          3, 4, 5
         3  cerrado/woody savanna   IGBP  6 thru  9
         4  Grassland/cropland      IGBP 10 thru 17

     the new attribute name is "veg". Notice that this module also
     defines the following constants:

         TROPICAL = 1
         EXTRA_TROPICAL = 2
         SAVANNA = 3
         GRASSLAND = 4

     corresponding to each aggregated vegetation type.
    """

    nobs = lon.shape[0]
    veg = IGBP_.getsimpleveg(lon,lat,Path+'/IGBP',nobs) # Fortran

    # substitute non vegetation (water, snow/ice) data with
    # another type, e.g., GRASSLAND by default
    i = (veg == 0)
    veg[i] = nonVeg # for now, waiting to hear from Saulo...

    return veg

#--
def getDetailedVeg(lon,lat,Path): 
    """
        Given the pathname for the IGBP datasets, returns
    information about the vegetation type.

    IGBP Land Cover Legend:

    Value     Description
    -----     -----------
      1       Evergreen Needleleaf Forest
      2       Evergreen Broadleaf Forest
      3       Deciduous Needleleaf Forest
      4       Deciduous Broadleaf Forest
      5       Mixed Forest
      6       Closed Shrublands
      7       Open Shrublands
      8       Woody Savannas
      9       Savannas
     10       Grasslands
     11       Permanent Wetlands
     12       Croplands
     13       Urban and Built-Up
     14       Cropland/Natural Vegetation Mosaic
     15       Snow and Ice
     16       Barren or Sparsely Vegetated
     17       Water Bodies
     99       Interrupted Areas (Goodes Homolosine Projection)
     100      Missing Data
    
    """
    nobs = lon.shape[0]
    veg = IGBP_.getdetailedveg(lon,lat,Path+'/IGBP',nobs) # Fortran
    return veg

#--
def getSimpleVegMap(Path, refine=4, res=None):
    """
        Given the pathname for the IGBP datasets, and desired 
    resolution, returns global QFED vegetation map.
    """

    # non vegetation value
    NON_VEGETATION = 0

    # output grid resolution
    refine_factor = {'a' : 1, 'b' : 2, 'c' : 4, 'd' : 8, 'e' : 16}

    if res in refine_factor:
        refine = refine_factor[res]    

    # lat&lon grid
    d_lon = 5. / refine
    d_lat = 4. / refine

    n_lon = int(360. / d_lon)
    n_lat = int(180. / d_lat + 1)

    lon = linspace(-180.0, 180.0, n_lon, endpoint=False)
    lat = linspace( -90.0,  90.0, n_lat)

    lat2d, lon2d = meshgrid(lat, lon)

    # vegetation map
    veg = zeros((n_lon, n_lat))

    for i in range(n_lon):
        veg[i, :] = getSimpleVeg(lon2d[i, :], lat2d[i, :], Path, nonVeg=NON_VEGETATION)

    # pack vegetation and geolocation data together
    map = {'data': veg, 'lat': lat, 'lon': lon}

    return map


if __name__ == "__main__":

    from optparse import OptionParser
    from numpy    import linspace, meshgrid, zeros

    from gfio     import GFIO
    from qfed.mxd14_l3 import __VERSION__ as QFED_VERSION


    # default values for the command line options
    IGBP_DIR = '/nobackup/Emissions/Vegetation/GL_IGBP_INPE'
    QFED_VEG_FILE = 'QFED.vegetation_map.x288_y181.2010.nc4'
    GRID_RESOLUTION = 'c'


    # parse command line options
    parser = OptionParser(usage="Usage: %prog [options]", version='qfed_vegtypes-1.0.b3')

    parser.add_option("-i", "--igbp", 
                      dest="igbp_dir", default=IGBP_DIR,
                      help="directory for IGBP vegetation database (default=%s)" % IGBP_DIR)
    
    parser.add_option("-o", "--output", 
                      dest="qfed_veg", default=QFED_VEG_FILE,
                      help="QFED vegetation map output NetCDF file (default=%s)" % QFED_VEG_FILE)

    parser.add_option("-r", "--resolution", 
                      dest="res", default=GRID_RESOLUTION,
                      help="horizontal resolution: 'a' for 4x5, 'b' for 2x2.5, etc. (default=%s)" % GRID_RESOLUTION)
    
    (options, args) = parser.parse_args()

    if len(args) > 0:
        parser.error("Usage: %prog [options].\nTry '%prog --help' for more information.")
    
    # generate the aggregated vegetation map
    veg_map = getSimpleVegMap(options.igbp_dir, res=options.res)

    f = GFIO()
       
    n_ymd = 20100609
    n_hms = 120000

    v_name  = ['vegetation', ]
    v_title = ['QFED v%3.1f vegetation map' % QFED_VERSION, ]
    v_units = ['', ]

    title = 'QFED v%3.1f Vegetation Map' % QFED_VERSION
    source = 'NASA/GSFC/GMAO GEOS-5 Aerosol Group'
    contact = ('%s; %s') % ('arlindo.dasilva@nasa.gov', 'anton.darmenov@nasa.gov')

    # create the output file
    try:
        f.create(options.qfed_veg,
                 v_name,
                 n_ymd,
                 n_hms,
                 lon=veg_map['lon'],
                 lat=veg_map['lat'],
                 vtitle=v_title,
                 vunits=v_units,
                 title=title,
                 source=source,
                 contact=contact)
    except:
        print("ERROR: Could not create QFED vegetation file.")
        sys.exit(1)

    # write the data
    try:
        f.write(v_name[0], n_ymd, n_hms, veg_map['data'])
    except:
        print("ERROR: Could not save the QFED vegetation file.")

        try:
            f.close()
        except:
            pass

        sys.exit(2)

    # close the file
    try: 
        f.close()
    except: 
        pass

    print("INFO: The QFED v%3.1f vegetation map was saved to '%s'\n" % (QFED_VERSION, options.qfed_veg))

