from netCDF4 import Dataset
import numpy as np
import math as m

def convert_to_cart(lon,lat):
    cos_lat=np.cos(lat)
        return(
           np.cos(lon)*cos_lat,
           np.sin(lon)*cos_lat,
           np.sin(lat))

class Geolocation(object):

   def __init__(self,gridfile):
       ncFid = Dataset(gridfile, mode='r', format='NETCDF4')
       if ("nf" in ncFid.dimensions):
          self.lon_corners=ncFid.variables['corner_lons'][:]
          self.lat_corners=ncFid.variables['corner_lats'][:]
          self.lon=ncFid.variables['lons'][:]
          self.lat=ncFid.variables['lats'][:]
          self.npts=self.lon_corners.shape[1]
       elif ("grid_size" in ncFid.dimensions):
          self.npts = len(ncFid.dimensions['grid_size'])
          self.npts=self.npts//6
          self.npts=int(np.sqrt(self.npts))+1
          self.lon_corners=np.empty([6,self.npts,self.npts])
          self.lat_corners=np.empty([6,self.npts,self.npts])
          self.lon=np.empty([6,self.npts-1,self.npts-1])
          self.lat=np.empty([6,self.npts-1,self.npts-1])

          fsz = (self.npts-1)*(self.npts-1)
          noff = 0

          tmplon=ncFid.variables['grid_corner_lon'][:]
          tmplat=ncFid.variables['grid_corner_lat'][:]
          noff=0
          for n in range(6):
              for j in range(self.npts-1):
                  for i in range(self.npts-1):
                       self.lon_corners[n,j,i]=tmplon[noff,0]
                       self.lon_corners[n,j,i+1]=tmplon[noff,1]
                       self.lon_corners[n,j+1,i+1]=tmplon[noff,2]
                       self.lon_corners[n,j+1,i]=tmplon[noff,3]

                       self.lat_corners[n,j,i]=tmplat[noff,0]
                       self.lat_corners[n,j,i+1]=tmplat[noff,1]
                       self.lat_corners[n,j+1,i+1]=tmplat[noff,2]
                       self.lat_corners[n,j+1,i]=tmplat[noff,3]
                       noff=noff+1

          tmplon=ncFid.variables['grid_center_lon'][:]
          tmplat=ncFid.variables['grid_center_lat'][:]
          noff=0
          for n in range(6):
              for j in range(self.npts-1):
                  for i in range(self.npts-1):
                       self.lon[n,j,i]=tmplon[noff]
                       self.lat[n,j,i]=tmplat[noff]
                       noff=noff+1

       else:
          raise Exception("Unsupported grid file type")

       self.lon_corners=self.lon_corners*m.pi/180.0
       self.lat_corners=self.lat_corners*m.pi/180.0
       self.lon=self.lon*m.pi/180.0
       self.lat=self.lat*m.pi/180.0
       # if we want to store as xyz we may have to transpose first that convert_to_cart
       # puts the xyz coords at the front
       #self.xyzc=np.array(convert_to_cart(self.lon_corners,self.lat_corners))
       #self.xyz=np.array(convert_to_cart(self.lon,self.lat))

   def getIndices(self,lon,lat):
       npoints = lon.size
       indices=[]

 
       # fortran goes here


       return (indices)
