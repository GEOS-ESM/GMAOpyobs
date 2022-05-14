#!/bin/env python
"""
   Implements Python interface to the OMI SO2 Level 2 product (OMSO2).   
"""

import h5py
from   numpy    import ones
from   datetime import date, datetime, timedelta

class OMSO2_L2(object):

      """Class for OMI SO2 objects."""


      def __init__ (self,file_name,keep=24):
        """
        Creates a OMI SO2 object defining the following attributes:

             lon            ---  longitudes in degrees        (nobs)
             lat            ---  latitudes in degrees         (nobs)

             solar_zenith   ---  solar zenith angle           (nobs)
             relat_azymuth  ---  relative azymuth angle       (nobs)
             sensor_zenith  ---  sensor zenith angle          (nobs)

       Ravi: please document opther attributes...

       The optional parameter *keep* is used to specify the number of scan
       lines (from the lefy of the swath) to keep. This is needed for
       coping with the row anomaly problem.

        """
        

#       SDS to Keep
#       -----------

        SDS = {'HDFEOS/SWATHS/OMI Total Column Amount SO2/Geolocation Fields':
                  ('Latitude','Longitude','RelativeAzimuthAngle','SolarAzimuthAngle',
                   'SolarZenithAngle','TerrainHeight','Time','Date'),
               'HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields':
                  ('RadiativeCloudFraction','CloudPressure','ColumnAmountO3',
                   'ColumnAmountSO2_PBL','ColumnAmountSO2_TRL','TerrainPressure',
                   'UVAerosolIndex','QualityFlags_PBL','QualityFlags_TRL')
               }

        # Short names
        # -----------
        Short_Name = dict(
                     CloudPressure = 'ocp',
                    ColumnAmountO3 = 'o3',
               ColumnAmountSO2_PBL = 'so2_pbl',
               ColumnAmountSO2_TRL = 'so2_trl',
                          Latitude = 'lat' ,
                         Longitude = 'lon' ,
                  QualityFlags_PBL = 'qa_pbl',
                  QualityFlags_TRL = 'qa_trl',
            RadiativeCloudFraction = 'rcf',
              RelativeAzimuthAngle = 'relat_azymuth',
                 SolarAzimuthAngle = 'solar_zenith',
                  SolarZenithAngle = 'solar_zenith',
                     TerrainHeight = 'zs',
                   TerrainPressure = 'ps',
                    UVAerosolIndex = 'ai',
                     )

        # Extract Strating Julian Time for this orbit
        # -------------------------------------------
        start = 2448988.5
        DATE_START = datetime(1993,1,1,0,0,0)

        # Open the OMI SO2 file and loop over the datasets and
        # extract GEOLOCATION and Data fields
        f = h5py.File(file_name)

        for sds in list(SDS.keys()):
          g = f.get(sds)

          for v in SDS[sds]:

            if v == 'Time': 

              time = g.get(v)
              im = time.shape[0]
              nymd  = ones(im)
              nhms  = ones(im)

              for i in range(im):

                dt = time[i]
                ns = timedelta(seconds=dt)
                omi_datetime = DATE_START + ns

                year    = omi_datetime.year
                month   = omi_datetime.month
                day     = omi_datetime.day
                hour    = omi_datetime.hour
                minutes = omi_datetime.minute
                seconds = omi_datetime.second

                nymd[i] = int('%4d%02d%02d' % (year,month,day))
                nhms[i] = int('%02d%02d%02d' % (hour,minutes,seconds))

              self.time = time[:]
              self.nymd = nymd
              self.nhms = nhms

            elif v != 'Date':             
                data = g.get(v)
                try:
                      name = Short_Name[v]
                except:
                      name = v
                self.__dict__[name] = data[:,0:keep]

            self.nlines, self.nscans = self.lat.shape

#---

      def writeODS(self,filename=None,dir='.',expid='omso2',channels=None,
                 Verb=1):
        """
        Writes the un-gridded OMI SO2 object to an ODS file. 
        """

        if filename is None:
           filename = '%s/%s.obs.%d_%02dz.ods'%(dir,expid,self.nymd,self.nhms/10000)

        if channels is None:
            channels = list(self.channels)

        # Create and populated ODS object
        # -------------------------------
        ns = self.nobs
        nobs = len(channels) * ns
        ods = ODS(nobs=nobs, kx=self.kx, kt=ktAOD)
        i = 0
        for ch in channels:
            I = list(range(i,i+ns))
            j = channels.index(ch)
            ods.ks[I]  = i+1
            ods.lat[I] = self.lat[:]
            ods.lon[I] = self.lon[:]
            ods.time[I] = self.time[:]
            ods.lev[I] = ch
            ods.qch[I] = self.qa_flag[:].astype('int')
            ods.obs[I] = self.aod[:,j]
            i += ns

        # Exclusion flag
        # --------------
        iGood = (ods.qch>0) & (ods.obs<10.)
        ods.qcx[:] = 1     # All bad...
        ods.qcx[iGood] = 0 # ... unless good

        ods.write(filename,self.nymd,self.nhms)
        
        if Verb >=1:
            print("[w] Wrote file "+filename)

#---
      def writeg(self,filename=None,dir='.',expid='omso2',refine=8,res=None,
               channels=None,Verb=1):
       """
        Writes gridded MODIS measurements to file.

         refine  -- refinement level for a base 4x5 GEOS-5 grid
                       refine=1  produces a   4  x  5    grid
                       refine=2  produces a   2  x2.50   grid
                       refine=4  produces a   1  x1,25   grid
                       refine=8  produces a  0.50x0.625  grid
                       refine=16 produces a  0.25x0.3125 grid
        Alternatively, one can specify the grid resolution with a
        single letter:

         res     -- single letter denoting GEOS-5 resolution,
                       res='a'  produces a   4  x  5    grid
                       res='b'  produces a   2  x2.50   grid
                       res='c'  produces a   1  x1,25   grid
                       res='d'  produces a  0.50x0.625  grid
                       res='e'  produces a  0.25x0.3125 grid

                   NOTE: *res*, if specified, supersedes *refine*.

         Verb -- Verbose level:
                 0 - really quiet (default)
                 1 - Warns if invalid file is found
                 2 - Prints out non-zero number of fires in each file.


       """

#      Output grid resolution
#      ----------------------
       if res is not None:
           if res=='a': refine = 1 
           if res=='b': refine = 2
           if res=='c': refine = 4
           if res=='d': refine = 8
           if res=='e': refine = 16

#      Lat lon grid
#      ------------
       dx = 5. / refine
       dy = 4. / refine
       im = int(360. / dx)
       jm = int(180. / dy + 1)

       glon = linspace(-180.,180.,im,endpoint=False)
       glat = linspace(-90.,90.,jm)

       if channels is None:
           channels = self.channels
       levs = array(channels)

       nch = len(channels)
       nymd = self.nymd
       nhms = self.nhms

       vtitle = [ 'Aerosol Optical Depth',
                  'Aerosol Optical Depth (Revised)',
                  'Aerosol Optical Depth (Fine Mode)',
                  'Cloud Fraction' ]

       vname  = ['tau', 'tau_', 'tau_fine', 'cloud' ]
       vunits = [ '1',    '1',     '1',       '1',  ]
       kmvar  = [ nch,    nch,     nch,        0    ]

       title = 'Gridded MODIS Measurements'
       source = 'NASA/GSFC/GMAO GEOS-5 Aerosol Group'
       contact = 'arlindo.dasilva@nasa.gov'

       if filename is None:
           filename = '%s/%s.obs.%d_%02dz.nc4'%(dir,expid,self.nymd,self.nhms/10000)

       # Create the file
       # ---------------
       f = GFIO()
       f.create(filename, vname, nymd, nhms,
                lon=glon, lat=glat, levs=levs, levunits='nm',
                vtitle=vtitle, vunits=vunits,kmvar=kmvar,amiss=MISSING,
                title=title, source=source, contact=contact)

       # Subset AOD at specified channels
       # --------------------------------
       I = []
       for ch in channels:
           i = list(self.channels).index(ch)
           I = I + [i,]
       aod = self.aod[:,I]
       aod_fine = self.aod_fine[:,I]
       
       # The Revised AOD may not exist
       # -------------------------------
       try:
           aod_ = self.aod_[:,I]
       except:
           aod_ = MISSING * ones(aod.shape) # will compress like a charm

       # Grid variable and write to file
       # -------------------------------
       f.write('tau', nymd, nhms, 
               binobs3d(self.lon,self.lat,aod,im,jm,MISSING) )
       f.write('tau_', nymd, nhms, 
               binobs3d(self.lon,self.lat,aod_,im,jm,MISSING) )
       f.write('tau_fine', nymd, nhms, 
               binobs3d(self.lon,self.lat,aod_fine,im,jm,MISSING) )
       f.write('cloud', nymd, nhms, 
               binobs2d(self.lon,self.lat,self.cloud,im,jm,MISSING) )
           
#       try:
#           f.close()
#       except:
#           pass

       if Verb >=1:
           print("[w] Wrote file "+filename)

#....................................................................

if __name__ == "__main__":

    q = OMSO2('OMI-Aura_L2-OMSO2_2011m0101t0111-o34379_v003-2011m0101t074655.he5')
    

