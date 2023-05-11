#!/bin/env python
"""
   Implements Python interface to the HSRL L2 data.
"""

import h5py
from   numpy    import ones, zeros, interp, NaN, isnan, array
from   datetime import datetime, timedelta

from .config import strTemplate

SDS_DIAL = {'Nav_Data':
                  ('gps_alt','gps_date','gps_lat','gps_lon','gps_time'),
            'Data_Products':
                  ('1064_bsc_cloud_screened',
                   '1064_ext',
                   '1064_total_attn_bsc_cloud_screened',
                   '1064_Single_Scattering_Albedo',
                   '532_bsc_cloud_screened',
                   '532_ext',
                   '532_total_attn_bsc_cloud_screened',
                   '532_Single_Scattering_Albedo',
                   'Altitude',
                   'Dust_Mixing_Ratio',
                   'ground_altitude',
                   'SeaSalt_Mixing_Ratio',
                   'SO4_Mixing_Ratio',
                   'OC_Mixing_Ratio',  
                   'BC_Mixing_Ratio'  ),
              'State_Parameters':
                  ( 'Pressure',
                    'Temperature',          
                    'Relative_Humidity',
                  )
               }

SDS_HSRL_daq = {'Engineering/header': ('date',), 
             'ApplanixIMU':
                        ('gps_alt','gps_lat','gps_lon','gps_time'),
             'DataProducts' : (
                         'Altitude',
                         '1064_bsc_cloud_screened',
                         '1064_ext',
                         '1064_aer_dep',
                         '532_bsc_cloud_screened',
                         '532_ext',
                         '532_aer_dep',
                         '355_bsc_cloud_screened',
                         '355_ext',
                         '355_aer_dep',
                         ),
              'State':
                         ( 'Pressure',
                           'Temperature',          
                           'Relative_Humidity',
                         )
               }

SDS_HSRL2 = {'header': ('date',),
             'ApplanixIMU':
#             'ER2_IMU':
                        ('gps_alt','gps_lat','gps_lon','gps_time'),
             'DataProducts':
                         (
                         'Altitude',
                         '1064_bsc_cloud_screened',
                         '1064_ext',
                         '1064_aer_dep',
                         '532_bsc_cloud_screened',
                         '532_ext',
                         '532_aer_dep',
                         '355_bsc_cloud_screened',
                         '355_ext',
                         '355_aer_dep',
                         ),
              'State':
                         ( 'Pressure',
                           'Temperature',          
                           'Relative_Humidity',
                         )
               }

NAV = ( 'Altitude','date', 'gps_date','gps_lat','gps_lon','gps_time')

# Short names
# -----------
Short_Name = dict(
                      gps_alt = 'lev',
                     gps_date = 'datex',
                      gps_lat = 'lat',
                      gps_lon = 'lon',
                     gps_time = 'time',
                     Altitude = 'z',
              ground_altitude = 'zs',
                     )

Short_Name['1064_bsc_cloud_screened'] = 'bsc_1064'
Short_Name['1064_ext'] = 'ext_1064'
Short_Name['1064_aer_dep'] = 'dep_1064'
Short_Name['1064_total_attn_bsc_cloud_screened'] = 'absc_1064'
Short_Name['1064_Single_Scattering_Albedo'] = 'ssa_1064'

Short_Name['532_bsc_cloud_screened'] = 'bsc_532'
Short_Name['532_ext'] = 'ext_532'
Short_Name['532_aer_dep'] = 'dep_532'
Short_Name['532_total_attn_bsc_cloud_screened'] = 'absc_532'
Short_Name['532_Single_Scattering_Albedo'] = 'ssa_532'

Short_Name['355_bsc_cloud_screened'] = 'bsc_355'
Short_Name['355_ext'] = 'ext_355'
Short_Name['355_aer_dep'] = 'dep_355'


class HSRL(object):

      """Implements the HSRL class."""


      def __init__ (self,hsrl_filename,Nav_only=False,verbose=True,freq=3.0,
                    zmax=10,SDS=SDS_HSRL2):
        """
        Creates an HSRL object defining the following attributes:

             lon            ---  longitudes in degrees        (nobs)
             lat            ---  latitudes in degrees         (nobs)

        """

        self.verb = verbose

        
        # Open the HSRL file and loop over the datasets,
        # extracting navigation and data
        # ----------------------------------------------
        f = h5py.File(hsrl_filename,mode='r+')
        if self.verb:
           print("[] Opening HSRL file <%s>"%hsrl_filename)

        # Read the variables proper
        # -------------------------
        for sds in list(SDS.keys()):
          g = f.get(sds)
          for v in SDS[sds]:
            if Nav_only and v not in NAV:
              continue
            if self.verb:
                  print("   + Reading <%s>"%v)
            data = g.get(v)
            try:
              name = Short_Name[v]
            except:
              name = v
            self.__dict__[name] = data

        # Read coords and make these numpy arrays
        # ---------------------------------------
        self.lon  = self.lon[:].flatten()
        self.lat  = self.lat[:].flatten()
        self.z = self.z[:].flatten()
        self.time = self.time[:].flatten()

        self.nt = self.lat.shape[0]
        self.nz = self.z.shape[0]

        # Reduce vertical levels
        # ----------------------
        self.K = self.z/1000 <= zmax
        self.z = self.z[self.K]
        
        # Handle incosistency of date across HSRL datasets
        # -----------------------------------------------
        if self.nt != self.date.shape[0]:
          date_ = self.date[0,0]
          yy = int(date_)/10000
          mm = (int(date_) - yy * 10000)/100
          dd = int(date_) - yy*10000 - mm*100
          dates = '%02d/%02d/%4d'%(mm,dd,yy)
          self.date = array([dates for i in range(self.nt)])
          
        # Create datetime
        # ---------------
        self.Time = []
        for i in range(self.nt):
           dt = timedelta(seconds = int(self.time[i]* 60. * 60.+0.5))
           mm, dd, yy = self.date[i].split('/')
           self.Time += [datetime(int(yy), int(mm), int(dd)) + dt,]
        self.Time = array(self.Time)
        self.tyme = self.Time.reshape((self.nt,1))
        
        # Find bracketing synoptic times
        # ------------------------------
        dt = timedelta(seconds=int(freq * 60. * 60.)) # 3 hour
        tmin, tmax = (self.Time[0], self.Time[-1]) 
        t_ = datetime(tmin.year,tmin.month,tmin.day)  # beg of day 
        s0 = int((tmin-t_).seconds / (freq*60*60) )   
        T  = t_ + s0 * dt                             # lower synoptic hour
        self.syn = [T,]
        while ( T<tmax ):
          T += dt
          self.syn += [ T, ]
        self.dt_syn = dt

        # Find time indices within synoptic brackets
        # ------------------------------------------
        self.IA = [] # index and weight for time interpolation
        for t in self.syn[:-1]:
          a = ones(self.nt)
          I = (a==1.)
          for n in range(self.nt):
             I[n] = (self.Time[n]>=t)&(self.Time[n]<t+self.dt_syn)
             a[n] = (self.Time[n] - t).seconds
          a = a / float(self.dt_syn.seconds) 
          self.IA += [(I,a),]
        
        # Save file handle for later
        # --------------------------
        self.h5 = f

#--
      def simul(self, asm_filename,
                      aer_filename,
                      ext532_filename,
                      ext1064_filename):

         """Simulate HSRL observables from GEOS-5"""

         # Read in met fields
         # ------------------
         Vars = ( 'T', 'H', 'RH', 'PHIS', 'PL' )
         self._attachVar(asm_filename,Vars)
         self.zs[:] = self.PHIS[:]/9.8     # model surface height

         # Vertical Interpolation: met fields
         # ----------------------------------
         self.Pressure[:,:] = self.zInterp(self.PL)
         self.Temperature[:,:] = self.zInterp(self.T)
         self.Relative_Humidity[:,:] = self.zInterp(self.RH)

         # Read in aerosol mixing ratio variables
         # --------------------------------------
         Vars = ( 'DU001', 'DU002', 'DU003', 'DU004', 'DU005', 
                  'SS001', 'SS002', 'SS003', 'SS004', 'SS005',
                  'BCPHOBIC', 'BCPHILIC', 'OCPHOBIC', 'OCPHILIC',
                  'SO4' )
         self._attachVar(aer_filename,Vars)

         # Combine bins
         # ------------
         DU = self.DU001 + self.DU002 + self.DU003 + self.DU004 + self.DU005 
         SS = self.SS001 + self.SS002 + self.SS003 + self.SS004 + self.SS005 
         BC = self.BCPHILIC + self.BCPHOBIC
         OC = self.OCPHILIC + self.OCPHOBIC

         # Vertical interpolation: aerosol mixing ratio
         # --------------------------------------------
         self.Dust_Mixing_Ratio[:,:] = self.zInterp(DU)
         self.SeaSalt_Mixing_Ratio[:,:] = self.zInterp(SS)
         self.BC_Mixing_Ratio[:,:] = self.zInterp(BC)
         self.OC_Mixing_Ratio[:,:] = self.zInterp(OC)
         self.SO4_Mixing_Ratio[:,:] = self.zInterp(self.SO4)

         # Extinction coefficients: 532 nm
         # -------------------------------
         Vars = ( 'extinction', 'ssa', 'backscat', 'aback_toa' )
         self._attachVar(ext532_filename,Vars)
         self.ext_532[:,:] = self.zInterp(self.extinction)
         self.ssa_532[:,:] = self.zInterp(self.ssa)
         self.bsc_532[:,:] = self.zInterp(self.backscat)
         self.absc_532[:,:] = self.zInterp(self.aback_toa)

         # Extinction coefficients: 1064 nm
         # --------------------------------
         Vars = ( 'extinction', 'ssa', 'backscat', 'aback_toa' )
         self._attachVar(ext1064_filename,Vars)
         self.ext_1064[:,:] = self.zInterp(self.extinction)
         self.ssa_1064[:,:] = self.zInterp(self.ssa)
         self.bsc_1064[:,:] = self.zInterp(self.backscat)
         self.absc_1064[:,:] = self.zInterp(self.aback_toa)

#--
      def zInterp(self,v5):
            """Vertically interpolates the GEOS-5 variable v5
            to the HSRL heights."""
            
            if len(v5.shape) != 2:
                  raise ValueError('variable to be interpolated must rave rank 2')

            nt, nz = v5.shape
             
            if self.nt != nt:
                  raise ValueError('inconsistent time dimension')

            if self.H.shape[1] != nz:
                  raise ValueError('inconsistent GEOS-5 vertical dimension')

            v = ones((self.nt,self.nz)) # same size as HSRL arrays
            for t in range(self.nt):
                  v[t,:] = interp(self.z,self.H[t,:],v5[t,:],left=NaN)

            return v

#--
      def save(self):
            """Save the modified file to disk, overwriting it."""
            self.h5.flush()
            self.h5.close()

#--
      def curtain(self,v, tit=None, vmin=None, vmax=None,figFile=None,
                  mask=None, mask_as=None, Log=False,xlab='Time (Hours UTC)',

                  ticks=None, format=None, aspect=0.205):
            """Plots a 2D curtain plot."""

            from matplotlib.pyplot import imshow, xlabel, ylabel, title, gca, savefig
            from matplotlib.colors import LogNorm

            z = self.z / 1000.
            t = self.time
            ext = ( t[0], t[-1], z[0], z[-1] )

            nt, nz = v.shape
            if nz==len(self.K):
                  v_ = v[:,self.K]
            else:
                  v_ = v[:,:]

            if mask_as != None:
                  mask = isnan(mask_as[:,:])
            if mask != None:
                  v_[mask] = NaN

            if vmin is not None:
                  v_[v_<vmin] = NaN
            if vmax is not None:
                  v_[v_>vmax] = NaN
                  
            gca().set_axis_bgcolor('black')

            if Log:
                  imshow(v_.T,origin='lower',extent=ext,aspect=aspect,
                         norm=LogNorm(vmin=vmin,vmax=vmax))
            else:
                  imshow(v_.T,origin='lower',extent=ext,aspect=aspect,vmin=vmin,vmax=vmax)

            # Labels
            # ------
            xlabel(xlab)
            ylabel('Height (km)')
            if tit != None: 
                  title(tit)

            _colorbar()
            # colorbar(shrink=0.78,ticks=ticks,format=format)

            if figFile is not None:
                  savefig(figFile,dpi=180)
            
#--
      def coords(self,text_filename):
          """Prints lat/lon/time to a text file."""
          f = open(text_filename,'w')
          for i in range(self.lon.shape[0]):
                f.write("%f,%f,%s\n"%(self.lon[i],self.lat[i],self.Time[i].isoformat()))
          f.close()

#--
      def _attachVar (self,g5_filename,Vars):
         """
         Reads a GEOS-5 variable and interpolate to HSRL track;
         If g5_filename is a GrADS template then a spatial & time
         interpolation is performed. Otherwise, only a spatial
         interpolation is performed.
         """
         if g5_filename.count('%') > 0:
               self._attachVarXYT(g5_filename,Vars)
         else:
               self._attachVarXY(g5_filename,Vars)

#--
      def _attachVarXY (self,g5_filename,Vars):
         """
         Reads a GEOS-5 variable and interpolate to HSRL track;'
         spatial interpolation only.
         """
         from gfio import GFIO
         g = GFIO(g5_filename)
         if self.verb:
            print("[] Opening GEOS-5 file <%s>"%g5_filename)
         for var in Vars:
            if self.verb:
               print("   - Reading <%s>"%var)
            v = g.interp(var,self.lon,self.lat)
            if len(v.shape) == 2:
               v = v[-1::-1,:].T # flip vertical and transpose
            self.__dict__[var] = v

      def _attachVarXYT (self,g5_filename,Vars):
         """
         Reads a GEOS-5 variable and interpolate to HSRL track, space
         and time interpolation.
         """
         from gfio import GFIO
 
         # Compute index and open file for each syn time
         # ---------------------------------------------
         G = []
         for t in self.syn:
           filename = strTemplate(g5_filename,dtime=t) # expand 
           if self.verb:
             print("[] Opening GEOS-5 file <%s>"%filename)
           g = GFIO(filename)
           G += [(g,t),]
           km = g.km # most likely 72

         # Loop over variables
         # -------------------
         for var in Vars:
         
           V = []
           for g, t in G:

             if self.verb:
               print("   - Reading <%s> at "%var, t)

             v = g.interp(var,self.lon,self.lat) 
             if len(v.shape) == 2:
                v = v[-1::-1,:].T # flip vertical and transpose

             V += [v,]

           # Interpolate in each synoptic interval
           # -------------------------------------
           s = 0
           v = NaN * zeros(V[0].shape)          
           for I, a in self.IA:
             a_ = a[I]
             if len(v.shape) == 1:
               v[I] = (1-a_) * V[s][I] + a_ * V[s+1][I]
             else:
               for k in range(km):
                  v[I,k] = (1-a_) * V[s][I,k] + a_ * V[s+1][I,k]
             s += 1

           self.__dict__[var] = v
               
#---------------------------------------------------------------------------------
def _colorbar():                                                                
    """ Draw a colorbar """                                                     
                                                                                
    from matplotlib.pyplot import gca, colorbar, axes
    import unicodedata
    
    ax = gca()                                                                  
    bbox = ax.get_position()                                                    
    l,b,w,h = bbox.bounds # mpl >= 0.98                                         
    cax = axes([l+w, b, 0.04, h]) # setup colorbar axes.
    cbar = colorbar(cax=cax)
    labels = [unicodedata.normalize('NFKD', item.get_text()).encode('ascii','ignore') for item in cbar.ax.get_yticklabels()]
    #labels = ['%3.3f' % float(item) for item in labels]
    labels = ['%g' % float(item) for item in labels]
    
    cbar.ax.set_yticklabels(labels)        
    cbar.ax.set_axis_bgcolor('white') 
    axes(ax)  # make the original axes current again

    #....................................................................

if __name__ == "__main__":

    from matplotlib.pyplot import clf, subplot, figure, show

    h = HSRL('/Users/adasilva/iesa/oracles/hsrl/HSRL2_ER2_20160916_R1.h5',
             zmax=8,
             Nav_only=False)

def Hold():

    figure(figsize=(11,5.5))
    subplot(211)
    h.curtain(h.ext_532,Log=False,vmin=0.05,tit='Aerosol Extinction at 532nm',
              figFile='test.png',xlab=' ' )
    subplot(212)
    h.curtain(h.ext_532,Log=False,vmin=0.05,tit='Aerosol Extinction at 532nm',
              figFile='test.png' )
    show()
    
def Hold2():

#    print "Loading HSRL:"
#    h = HSRL('discoveraq-HSRL_UC12_20110705_RA_L2_sub.h5')
    print("Loading GEOS-5:")
    g = HSRL('discoveraq-geos5-HSRL_MODEL_20110705_RA_L2_sub.h5')

    g.simul('e572_fp.inst3_3d_asm_Nv.%y4%m2%d2_%h200z.nc4',
            'e572_fp.inst3_3d_aer_Nv.%y4%m2%d2_%h200z.nc4',
            'e572_fp.inst3_3d_ext532_Nv.%y4%m2%d2_%h200z.nc4',
            'e572_fp.inst3_3d_ext1064_Nv.%y4%m2%d2_%h200z.nc4' )

    g.simul('e572_fp.inst3_3d_asm_Nv.20110705_1800z.nc4',
            'e572_fp.inst3_3d_aer_Nv.20110705_1800z.nc4',
            'e572_fp.inst3_3d_ext532_Nv.20110705_1800z.nc4',
            'e572_fp.inst3_3d_ext1064_Nv.20110705_1800z.nc4' )

    g._attachVarT('e572_fp.inst3_3d_asm_Nv.%y4%m2%d2_%h200z.nc4',
                  ('H', 'RH', 'T',) )

    g.save()
    
