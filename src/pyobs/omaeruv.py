"""
   Implements Python interface to the OMI aerosol product OMAERUV.
"""

import os
from types    import *

import h5py
from   numpy    import ones, concatenate, linspace, reshape
from   datetime import date, datetime, timedelta

from .aura import AURA_L2, orbits, MISSING

MISSING = -1.267651e+30

CHANNELS = [354., 388., 471.]

kxOMAERUV = 314
ktAOD  = 45
ktAAOD = 49 

SDS = {'HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields':
                    ('Latitude','Longitude','RelativeAzimuthAngle','ViewingZenithAngle',
                     'SolarZenithAngle','TerrainPressure','Time','GroundPixelQualityFlags'),
       'HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields':
                     ('FinalAerosolLayerHeight','FinalAerosolOpticalDepth',
                      'FinalAerosolAbsOpticalDepth',
                      'ImaRefractiveIndex','FinalAlgorithmFlags','UVAerosolIndex',
                      'SurfaceAlbedo','Reflectivity',
                      'NormRadiance','AerosolType')
       }


class OMAERUV_L2(AURA_L2):

    """
    Implements interface to the OMI UV aerosol product (OMAERUV).
    """
      
    def __init__ (self,Path,SDS=SDS,keep=None,Verbose=0,only_good=True):
        """
        Creates an AURA object defining the attributes corresponding
        to the SDS's on input.

        The optional parameter *keep* is used to specify the number of scan
        lines (from the left of the swath) to keep. This is needed for
        coping with the row anomaly problem.

        """
        AURA_L2.__init__ (self,Path,SDS,keep,Verbose,only_good)
        
        self.channels = CHANNELS
        self.nch = len(self.channels)
        self.kx = kxOMAERUV


#---
    def writeODS(self,filename=None,dir='.',expid=None,channels=None):
        """
        Writes the un-gridded AURA object to an ODS file. 
        """
        from pyods    import ODS # must be inported before pyhdf

        if expid is None:
            expid = self.expid

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
    def ibits(self, arg, pos, len):
        return (arg >> pos) & ((1 << len) - 1)

#---
    def writeg(self,syn_time,nsyn=8,filename=None,dir='.',expid='omaeruv',refine=4,res=None,Verb=1):
       """
        Writes gridded OMI measurements to file.

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

       from gfio     import GFIO
       from binObs_  import binobs2d, binobs3d

       # Determine synoptic time range
       # -----------------------------
       dt = timedelta(seconds = 12. * 60. * 60. / nsyn)
       t1, t2 = (syn_time-dt,syn_time+dt)

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

       nch = self.nch
       nymd = 10000 * syn_time.year + 100 * syn_time.month  + syn_time.day
       nhms = 10000 * syn_time.hour + 100 * syn_time.minute + syn_time.second

       vtitle = [ 'Albedo',
                  'UV Aerosol Index',
                  'Radiance', 
                  'Total Aerosol Optial Depth',
                  'Absorption Aerosol Optial Depth' ]

       vname  = ['albedo','ai', 'rad', 'tauext', 'tauabs' ]
       vunits = [ '%',    '1',  '1',   '1',      '1' ]
       kmvar  = [nch,      0,    nch,    nch,   nch ]

       title = 'Gridded OMI Measurements and VLIDORT Simulations from GEOS-5'
       source = 'NASA/GSFC/GMAO GEOS-5 Aerosol Group'
       contact = 'arlindo.dasilva@nasa.gov'

       if filename is None:
           filename = '%s/%s.obs_l3a.%d_%02dz.nc4'%(dir,expid,nymd,nhms/10000)

#      Create the file
#      ---------------
       f = GFIO()
       f.create(filename, vname, nymd, nhms,
                lon=glon, lat=glat, levs=self.channels, levunits='nm',
                vtitle=vtitle, vunits=vunits,kmvar=kmvar,amiss=MISSING,
                title=title, source=source, contact=contact)

       # QA filtering
       # ------------
       I_bad = (self.FinalAlgorithmFlags!=0) # bad data
       
       # Time filter of data
       # -------------------
       lon = self.lon.ravel()
       lat = self.lat.ravel()
       albedo = _timefilter(self.time,t1,t2,self.albedo,I_bad).reshape((-1,nch))
       ai = _timefilter(self.time,t1,t2,self.ai,I_bad).ravel()
       radiance = _timefilter(self.time,t1,t2,self.radiance,I_bad).reshape((-1,nch))
       aod = _timefilter(self.time,t1,t2,self.aod,I_bad).reshape((-1,nch))
       aaod = _timefilter(self.time,t1,t2,self.aaod,I_bad).reshape((-1,nch))

       
#      Grid variable and write to file
#      -------------------------------
       f.write('albedo', nymd, nhms, binobs3d(lon,lat,albedo,  im,jm,MISSING) )
       f.write('ai',     nymd, nhms, binobs2d(lon,lat,ai,      im,jm,MISSING) )
       f.write('rad',    nymd, nhms, binobs3d(lon,lat,radiance,im,jm,MISSING) )
       f.write('tauext', nymd, nhms, binobs3d(lon,lat, aod,    im,jm,MISSING) ) 
       f.write('tauabs', nymd, nhms, binobs3d(lon,lat,aaod,    im,jm,MISSING) ) 

       if Verb >=1:
           print("[w] Wrote file "+filename)
           
#....................................................................

def _timefilter ( t, t1, t2, a, I_bad ):
    filler = MISSING * ones(a.shape[1:])
    b = a.copy()
    for i in range(len(t)):
        if (t[i]<t1) or (t[i]>=t2):
            b[i] = filler
    if len(b.shape) == 3:
        b[I_bad,:] = MISSING
    elif len(b.shape) == 2:
        b[I_bad] = MISSING
    else:
        raise IndexError("Invalid rank=%d for time filtering"%len(b.shape))
    return b

#..........................................................................

if __name__ == "__main__":

    omi_fn = '/nobackup/OMI/Level2/OMAERUV/2008/06/30/OMI-Aura_L2-OMAERUV_2008m0630t0233-o21055_v003-2009m0305t115716.he5'
    g5_fn = '/nobackup/ARCTAS/opendap/arctas.ddf'
    npz_fn = 'aer_Nv.npz'

    onlyVars = ['ps', 'delp', 'LWI', 'AIRDENS', 'RH',
                'du001', 'du002', 'du003', 'du004', 'du005',
                'ss001', 'ss002', 'ss003', 'ss004', 'ss005',
                'BCphobic', 'BCphilic', 'OCphobic', 'OCphilic']
    
    omi = OMAERUV_L2(omi_fn,Verbose=1)
    omi.sampleFile(g5_fn,onlyVars=onlyVars,npzFile=npz_fn,Verbose=1)

def hold():
    syn_time = datetime(2008,0o6,30,12,0,0)
    files = orbits("/nobackup/OMI/Level2","OMAERUV",syn_time,Verbose=0)
    q = OMAERUV_L2(files,Verbose=1)

