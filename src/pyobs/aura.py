#!/bin/env python
"""
   Implements Python interface to several AURA datasets.
"""

import os
from types    import *
from glob     import glob

import h5py
from   numpy    import ones, concatenate, savez, load, array, tile
from   datetime import date, datetime, timedelta

MISSING = -1.267651e+30
kxOMAERUV = 314
kxOMSO2   = 323 

ktAOD  = 45
ktAAOD = 49 

SDS = dict (

      OMSO2 = {'HDFEOS/SWATHS/OMI Total Column Amount SO2/Geolocation Fields':
                  ('Latitude','Longitude','RelativeAzimuthAngle','SolarAzimuthAngle',
                   'SolarZenithAngle','TerrainHeight','Time'),
                'HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields':
                  ('RadiativeCloudFraction','CloudPressure','ColumnAmountO3',
                   'ColumnAmountSO2_PBL','ColumnAmountSO2_TRL','TerrainPressure',
                   'UVAerosolIndex','QualityFlags_PBL','QualityFlags_TRL')
               },

      OMAERUV = {'HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields':
                    ('Latitude','Longitude','RelativeAzimuthAngle','ViewingZenithAngle',
                     'SolarZenithAngle','TerrainPressure','Time','GroundPixelQualityFlags'),
                  'HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields':
                     ('FinalAerosolLayerHeight','FinalAerosolOpticalDepth',
                      'FinalAerosolAbsOpticalDepth',
                      'ImaRefractiveIndex','FinalAlgorithmFlags','UVAerosolIndex',
                      'SurfaceAlbedo','Reflectivity',
                      'NormRadiance','AerosolType')
                },
             )
ALIAS = dict (
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
                 SolarAzimuthAngle = 'solar_azimuth',
                  SolarZenithAngle = 'solar_zenith',
                     TerrainHeight = 'zs',
                   TerrainPressure = 'ps',
                    UVAerosolIndex = 'ai',
                ViewingZenithAngle = 'sensor_zenith',
           FinalAerosolLayerHeight = 'zs',
          FinalAerosolOpticalDepth = 'aod',
       FinalAerosolAbsOpticalDepth = 'aaod',
               FinalAlgorithmFlags = 'qa_flag',
                     SurfaceAlbedo = 'albedo',
                      NormRadiance = 'radiance',
           GroundPixelQualityFlags = 'qa_grd',
              )

#........................................................................

class AURA_L2(object):

    """
    Base class for generic AURA objects.
    """
      
    def __init__ (self,Path,SDS,keep=None,Verbose=0,only_good=True):
        """
        Creates an AURA object defining the attributes corresponding
        to the SDS's on input.

        The optional parameter *keep* is used to specify the number of scan
        lines (from the left of the swath) to keep. This is needed for
        coping with the row anomaly problem.

        """
        
        # Initially are lists of numpy arrays for each granule
        # ----------------------------------------------------
        self.verb = Verbose
        self.keep = keep
        self.SDS = SDS

        # Variable names
        # --------------
        self.Names = []
        for group in list(self.SDS.keys()):
            for name in self.SDS[group]:
                self.Names.append(name)
        self.Names += ['nymd','nhms']

        # Create empty lists for SDS to be read from orbit file;
        #  each element of the list contains data for one orbit
        # ------------------------------------------------------
        for name in self.Names:
            self.__dict__[name] = []
        self.time = [] # to hold datetime objects
        
        # Read each orbit, appending them to the list
        # -------------------------------------------
        if type(Path) is ListType:
            if len(Path) == 0:
                self.nobs = 0
                print("WARNING: Empty AURA object created")
                return
        else:
            Path = [Path, ]
           
        self._readList(Path)

        # Make each attribute a single numpy array
        # ----------------------------------------
        for name in self.Names:
            try:
                self.__dict__[name] = concatenate(self.__dict__[name])
            except:
                print("Failed concatenating "+name)

        # Determine index of "good" observations
        # --------------------------------------
        pass # to do

        # Keep only "good" observations
        # -----------------------------
        if only_good:
            pass

        # Make aliases for compatibility with older code 
        # ----------------------------------------------
        Alias = list(ALIAS.keys())
        for name in self.Names:
            if name in Alias:
                self.__dict__[ALIAS[name]] = self.__dict__[name] 

        # ODS friendly attributes
        # -----------------------
        pass # todo

        self.time = array(self.time)
        self.tyme = self.time # for consistency with other packages

#---
    def _readList(self,List):
        """
        Recursively, look for files in list; list items can
        be files or directories.
        """
        for item in List:
            if os.path.isdir(item):      self._readDir(item)
            elif os.path.isfile(item):   self._readOrbit(item)
            else:
                print("%s is not a valid file or directory, ignoring it"%item)

#---
    def _readDir(self,dir):
        """Recursively, look for files in directory."""
        for item in os.listdir(dir):
            path = dir + os.sep + item
            if os.path.isdir(path):      self._readDir(path)
            elif os.path.isfile(path):   self._readOrbit(path)
            else:
                print("%s is not a valid file or directory, ignoring it"%item)

#---
    def _readOrbit(self,filename):
        """Reads one OMRUV/OMSO2 granule with Level 2 aerosol data."""

        # Reference time
        # --------------
        REF_DATE = datetime(1993,1,1,0,0,0)

        # Open the AURA file and loop over the datasets,
        # extracting GEOLOCATION and Data fields
        # ----------------------------------------------
        if self.verb:
            print("[] working on <%s>"%filename)
        f = h5py.File(filename,mode='r')
            
        for group in list(self.SDS.keys()):

          g = f.get(group)

          for v in self.SDS[group]:

            if v == 'Time': 

              Time = g.get(v)
              nobs = Time.shape[0]
              nymd  = ones(nobs).astype('int')
              nhms  = ones(nobs).astype('int')

              for i in range(nobs):

                t_secs = Time[i]
                n_secs = timedelta(seconds=t_secs)
                t = REF_DATE + n_secs

                yy, mm, dd = (t.year,t.month,t.day)
                h, m, s = (t.hour,t.minute,t.second)

                nymd[i] = 10000 * yy + 100 * mm + dd
                nhms[i] = 10000 * h  + 100 * m  + s
                self.time.append(t)
                
              self.Time.append(Time[:]) # time as on file
              self.nymd.append(nymd)
              self.nhms.append(nhms)

            else:

                data = g.get(v)

                if self.keep != None:
                    self.__dict__[v].append(data[:,0:self.keep])
                else:
                    self.__dict__[v].append(data[:,:])

    def sampleFile(self, inFile, npzFile=None, onlyVars=None, Verbose=True):
        """
        Interpolates all variables of inFile and optionally
        save them to file *npzFile*. On input, I is an optional filtering
        index.
        """
        from gfio import GFIOctl, GFIOHandle

        # Instiate grads and open file
        # ----------------------------
        f = GFIOctl(inFile)

        # Check if all variables in file are upper case
        # ---------------------------------------------
        lc = False
        for v in f.vname:
            if v != v.upper(): lc = True

        self.sample = GFIOHandle(inFile)
        if onlyVars is None:
            onlyVars = f.vname

        # Flatten coordinates
        # -------------------
        nt, nr = self.lon.shape
        tyme = tile(self.time,(nr,1)).T
        lons = self.lon.ravel()
        lats = self.lat.ravel()
        tymes = tyme.ravel()
             
        # Loop over variables on file
        # ---------------------------
        for v in onlyVars:
            if Verbose:
                print("<> Sampling ", v)
            if not lc:
                v_ = v.upper()
            else:
                v_ = v
            var = f.sample(v_,lons,lats,tymes,Verbose=Verbose)
            if len(var.shape) == 1:
                self.sample.__dict__[v] = var.reshape((nt,nr))
            elif len(var.shape) == 2:
                var = var.T # shape should be (nobs,nz)
                self.sample.__dict__[v] = var.reshape((nt,nr,-1))
            else:
                raise IndexError('variable <%s> has rannk = %d'%len(var.shape))

        if npzFile is not None:
            savez(npzFile,**self.sample.__dict__)            

    def sampleLoadz(self,npzFile):
        """
        Loads sample from npz file.
        """
        from grads.gahandle import GaHandle
        self.sample = GaHandle(npzFile)
        npz = load(npzFile)
        for v in list(npz.keys()):
            self.sample.__dict__[v] = npz[v]

    def sampleReduce(self,I,npzFile=None):
        """
        Reduce sampled arrays based on the logical indices I;
        notice that shape(I) = (nt,nr).
        """
        sample = self.sample
        self.iGood = I # save for later relating sample to swath coordinates
        for v in sample.__dict__:
            if v != 'name':
               rank = len(sample.__dict__[v].shape)
               if rank > 1:
                  sample.__dict__[v] = sample.__dict__[v][I]
        if npzFile is not None:
            savez(npzFile,iGood=self.iGood,**self.sample.__dict__)  

    def __copy__(self,Indices=None):
        """Implements copy operation; if Indices are specified, a subset
        based on these indices are returned."""


        nobs = self.lat.shape[0]
        Alias = list(ALIAS.keys())

#     Copy data attributes
#     --------------------
        if nobs > 0:

           if Indices is None:
              I = list(range(nobs))
           else:
              I = Indices
         
           for name in self.Names:
               rank = len(self.__dict__[name].shape)
               if rank > 1:
                  self.__dict__[name] = self.__dict__[name][I]

               if name in Alias:   # alias too
                  self.__dict__[ALIAS[name]] = self.__dict__[name]
#---

def orbits ( path, prod, syn_time, nsyn=8, Verbose=0 ):
    """
    Returns a list of OMI orbits for a given product at given synoptic time.
    On input,

    path      ---  mounting point for the MxD04 Level 2 files
    prod      ---  product: OMAERUV, OMSO2
    syn_time  ---  synoptic time (timedate format)

    nsyn      ---  number of synoptic times per day (optional)

    """

    # Determine synoptic time range
    # -----------------------------
    dt = timedelta(seconds = 12. * 60. * 60. / nsyn)
    t1, t2 = (syn_time-dt,syn_time+dt)

    # print "[*] ", t1,"|", t2

    today = syn_time
    yesterday = today - timedelta(hours=24)

    Files = []
    for t in (yesterday,today):
        
        yy, mm, dd = (t.year,t.month,t.day)
        doy = t.timetuple()[7]

        if prod == 'OMAERUV':
            dirn = "%s/%s/%4d/%02d/%02d"%(path,prod,yy,mm,dd)
        elif prod == 'OMSO2':
            dirn = "%s/%s/%4d/%03d"%(path,prod,yy,doy)
        else:
            raise ValueError("Invalid product <%s>"%prod)

        Files += glob("%s/OMI-Aura_L2-%s*.he5"%(dirn,prod))

    Orbits = []
    for f in Files:
        dirn, filen = os.path.split(f)
        tokens = filen.split('-')
        beg_yy = int(tokens[2].split('m')[0].split('_')[1])
        beg_mm = int(tokens[2].split('m')[1][0:2])
        beg_dd = int(tokens[2].split('m')[1][2:4])
        beg_h = int(tokens[2].split('t')[1][0:2])
        beg_m = int(tokens[2].split('t')[1][2:4])
        end_h = int(tokens[4].split('t')[1].split('.')[0][0:2])
        end_m = int(tokens[4].split('t')[1].split('.')[0][2:4])
        end_yy = int(tokens[4].split('m')[0])
        end_mm = int(tokens[4].split('m')[1][0:2])
        end_dd = int(tokens[4].split('m')[1][2:4])
        orb = tokens[3].split('_')[0][1:]
        t_beg = datetime(beg_yy,beg_mm,beg_dd,beg_h,beg_m,0)
        t_end = t_beg + timedelta(minutes=90)
        # t_end = datetime(end_yy,end_mm,end_dd,end_h,end_m,0)
        if (t_beg>=t1 and t_beg<t2) or (t_end>=t1 and t_end<t2):
            #print "[x] ", t_beg, '|', t_end, tokens
            Orbits += [f,]
            if Verbose:
                print("[] ", f)

    return Orbits

#............................................................................

class OMSO2_L2(AURA_L2):

    """
    Implements interface to the OMI SO2 product (OMSO2).
    """
      
    def __init__ (self,Path,SDS=SDS['OMSO2'],keep=24,Verbose=0,only_good=True):
        """
        Creates an AURA object defining the attributes corresponding
        to the SDS's on input.

        The optional parameter *keep* is used to specify the number of scan
        lines (from the left of the swath) to keep. This is needed for
        coping with the row anomaly problem.

        """
        AURA_L2.__init__ (self,Path,SDS,keep,Verbose,only_good)

#........................................................................

def _timeConv(str):
    year, tmp = str.split('m')
    month = tmp[0:2]
    day = tmp[2:4]
    hour = tmp[5:7]
    min = tmp[7:10]
    t = datetime(year=int(year),month=int(month),day=int(day),
                 hour=int(hour), minute=int(min))

    return t
    
def parseFilename ( path ):
    """
    Parses AURA file name.
    """
    
    tokens = os.path.basename(path).split('-')
    instr = tokens[0]
    sat, lev = tokens[1].split('_')
    prod, vtime = tokens[2].split('_')
    orb, ver = tokens[3].split('_')
    ptime = tokens[4]
    t =_timeConv(vtime)

    return instr, sat, lev, prod, vtime, orb,ver, ptime, t

#........................................................................

if __name__ == "__main__":

#    q = OMSO2_L2('OMI-Aura_L2-OMSO2_2010m0705t1640-o31767_v003-2010m0706t012444.he5')
###  q = OMSO2_L2('/nobackup/OMI/Level2/SO2/2010/186',Verbose=1)
#    q = OMAERUV_L2('OMI-Aura_L2-OMSO2_2010m0705t1640-o31767_v003-2010m0706t012444.he5')
#    q = OMAERUV_L2('/nobackup/OMI/Level2/2010/07/05',Verbose=1)

     syn_time = datetime(2010,0o6,30)
     for syn_hour in (0,3,6,9,12,15,18,21):
         syn_time = datetime(2010,0o6,30,syn_hour,0,0)
         print("Synoptic time: ", syn_time) 

