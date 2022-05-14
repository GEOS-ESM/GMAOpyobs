"""
   Class for reading MINX ASCII files with MISR plume height information.
"""

import os

from types import *
from numpy import loadtxt, ones, median, array, load, savez
from datetime import datetime
from glob  import glob

from numpy import percentile as prctile

from .config  import Config
from .kde     import calc_kde1d

Months = ('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')

def date2nymd(s):
    S = s.split('/')
    mm = int(S[0])
    dd = int(S[1])
    yy = 2000+int(S[2]) 
    return (yy,mm,dd)

def time2nhms(s):
    S = s.split(':')
    h = int(S[0])
    m = int(S[1])
    return (h,m)

def gatime(date,time):
    yy,mm,dd = date2nymd(date)
    return "%sZ%d%s%d"%(time,dd,Months[mm-1],yy)

# Columns in the MINX files
# -------------------------
Vars1 = ('Location','Longitude','Latitude','Block','Sample','Line','Distance','Orientation',\
        'Elevation','PlumeHeight_NoWind','PlumeHeight', 'Uacross', 'Ualong', \
        'alb_b','alb_g','alb_r','alb_n', 'alb_toa', \
        'tau_b','tau_g','tau_r','tau_n', \
        'ssa_b','ssa_g','ssa_r','ssa_n', \
        'f_small','f_medium','f_large', 'f_sphere', 'angstrom',
        'pow_MW', 'refl', 'B21', 'B31', 'E21', 'E31', )

Vars2 = ('Location','Longitude','Latitude','Block','Sample','Line',
         'Distance','Orientation','Elevation',
         'PlumeHeight_NoWind','PlumeHeight', 'filtered', 
         'Uacross', 'Ualong', 'Utotal',
         'tau_b','tau_g','tau_r','tau_n', 
         'ssa_b','ssa_g','ssa_r','ssa_n', 
         'f_small','f_medium','f_large', 'f_sphere', 'angstrom',
         'pow_MW',)


Fmts1 = (  'i4',       'f4',       'f4',     'i4',    'i4',   'i4', )   \
     + 31 * ('f4',)

Fmts2 = (  'i4',       'f4',       'f4',     'i4',    'i4',   'i4', )   \
     + 23 * ('f4',)

MISSING = -99.99

# New, short names
# ----------------
Alias = {}
for name in Vars1+Vars2:
    Alias[name] = name # by default same name

Alias['Longitude'] = 'lon'
Alias['Latitude'] = 'lat'
Alias['Elevation'] = 'zs'
Alias['PlumeHeight_NoWind'] = 'z0'
Alias['PlumeHeight'] = 'z'

class MINX(object):
    """Reads CSV files into Numpy arrays"""

    def __init__ (self,fname):
        """
        Loads a CSV file, given the name and format of each column.
        """
        
#       Load useful attributes
#       ----------------------
        cf = Config(fname)
        RC = dict()
        for rc in list(cf.keys()): RC[rc.upper()] = rc # handle change in case
        self.region = cf(RC['REGION NAME'])
        self.date = cf(RC['DATE ACQUIRED'])
        self.time = cf(RC['UTC TIME'])
        self.area = float(cf(RC['AREA (SQ KM)']))
        self.perimeter = float(cf(RC['PERIMETER LENGTH (KM)']))
        try:
            self.frp = float(cf(RC['POWER OF FIRE IN MW']))
        except:
            self.frp = MISSING
        self.zm = float(cf(RC['BEST MEDIAN HT (M ASL)']))
        self.zt = float(cf(RC['BEST TOP HT (M ASL)']))
        # self.qa = cf(RC['REGION DATA QUALITY'])

#       Python time
#       -----------
        tyme_ = [int(x) for x in self.date.split('-') + self.time.split(':')]
        self.tyme = datetime(*tyme_)
        
#       Determine number of rows to skip
#       --------------------------------
        n_skip = 0
        for i in range(len(cf.Lines)):
            line = cf.Lines[i]
            if 'Pt#' in line:
                n_skip = i+2


        if ( len(cf.Lines[n_skip+2].split()) == len(Vars2) ):
            Vars = Vars2
            Fmts = Fmts2
        else:
            Vars = Vars1
            Fmts = Fmts1

#       Read the data
#       -------------
        data = loadtxt(fname, 
                       dtype={'names':Vars,'formats':Fmts},
                       skiprows=n_skip )
        N = len(data)

#       Save data columns as attributes, with nicer names
#       -------------------------------------------------
        for i in range(len(Vars)):
            if Fmts[i]=='f4' or Fmts[i]=='i4':
                v = ones(N)
                for j in range(N):
                    v[j] = data[j][i]
            else:
                v = []
                for j in range(N):
                    v.append(data[j][i])

            self.__dict__[Vars[i]] = v

#       Keep only valid plume heights, create aliases
#       ---------------------------------------------
        self.iGood = (self.PlumeHeight-self.Elevation>0) & (self.PlumeHeight_NoWind-self.Elevation>0)
        for i in range(len(Vars)):
            v = self.__dict__[Vars[i]] # [iGood]
            self.__dict__[Vars[i]] = v
            if Alias[Vars[i]] != Vars[i]:
                self.__dict__[Alias[Vars[i]]] = v # alias

        self.N = self.z.size 

#---
    def add3Vars(self,ga,fname,vars=('u10m','v10m','slp')):
        """
        Given a grads object having the correct MERRA file as default, writes out 
        a CSV file with the 3 variables.
        """
        f = open(fname,"w")
        print("Location,Date,DOY,Time,Longitude,Latitude,%s,%s,%s"%vars, file=f)
        print("Location,Date,DOY,Time,Longitude,Latitude,%s,%s,%s"%vars)

        N = self.N
        U = ones(N)
        V = ones(N)
        P = ones(N)

        for i in range(N):

            d = self.data[i]
            x = self.lon[i]
            y = self.lat[i]
            t = self.time[i]

            ga('set lon %f %f'%(x-1.,x+1.),Quiet=True)
            ga('set lat %f %f'%(y-1.,y+1.),Quiet=True)
            ga('set time %s'%t,Quiet=True)

            u, levs = ga.interp(vars[0], lons=(x,),lats=(y,))
            v, levs = ga.interp(vars[1], lons=(x,),lats=(y,))
            p, levs = ga.interp(vars[2], lons=(x,),lats=(y,))

            r = (d[0],d[1],d[2],d[3],d[4],d[5],u.data,v.data,p.data)
            print("%s,%s,%d,%s,%6.3f,%6.3f,%f,%f,%f"%r, file=f)
            print("%s,%s,%d,%s,%6.3f,%6.3f,%f,%f,%f"%r)

            U[i] = u.data
            V[i] = v.data
            P[i] = p.data

        self.__dict__[vars[0]] = U
        self.__dict__[vars[1]] = V
        self.__dict__[vars[2]] = P

        f.close()

#....................................................................

class MINXs(object):
    """
    Reads 1 or more MINX files and create an object with the average
    properties of each plume.
    """

    def __init__(self,minxFiles,N_min=5):
        """
        Given a list of MINX plume files, retain the main average
        properties for each plume.
        """

        if type(minxFiles) == StringType:
            #minxFiles = sorted(glob(minxFiles))
            minxFiles = glob(minxFiles) # should be sorted, must be fixed after AGU (ams)

        Atts = ['file','nobs','qa','lon', 'lat', 'lon_f', 'lat_f','zs','z', 'zm', 'zt','tyme', 'frp']

        for a in Atts:
            self.__dict__[a] = []

        for minxFile in minxFiles:

            m = MINX(minxFile)

            I = m.iGood
            N = len(m.lon[I])
            if N<N_min: continue         # must have enough valid points in plume

            k = m.pow_MW.argmax()        # Coordinates of stronger fire
            if m.pow_MW[k]<0: 
                print('-- No fires for plume <%s>'%os.path.basename(minxFile))
                continue

            self.file.append(minxFile)
            self.tyme.append(m.tyme)
            self.frp.append(m.frp)
            self.nobs.append(m.N)
            #self.qa.append(m.qa='GOOD')

            self.zs.append(m.zs[k])       # Terrain height near stronger fire
                                          # Rationale: we will model PR at fire location

            #h = m.z[I]-m.zs[I]           # IMPORTANT: above surface!!!!
            h = m.z[I]-m.zs[k]            # IMPORTANT: above surface at stronger fire
            j = h.argsort()[N/2]
            self.lon.append(m.lon[j])     # coordinates of median fire     
            self.lat.append(m.lat[j])
            self.zm.append(h[j]) 
            
            prc = prctile(h,p=(0,5.,50,95.,100.))
            self.zt.append(prc[3])

            bins, P = calc_kde1d(h) 
            self.z.append(bins[P.argmax()])

            self.lon_f.append(m.lon[k]) # coordinates nearest stronger fire
            self.lat_f.append(m.lat[k])
            
        for a in Atts:
            self.__dict__[a] = array(self.__dict__[a]) 

        self.N = self.lon.size
        self.sample = None

#---
    def getMERRA(self, npzFile=None):
        """
        Get Met fields from MERRA.
        """
        from grads import GrADS

        ga = GrADS(Window=False,Echo=False)
    
        ga.open('http://goldsmr2.sci.gsfc.nasa.gov:80/dods/MAT1NXFLX') # 2D fluxes
        self.sampleFile(ga,onlyVars=('ustar', 'bstar', 'pblh', 'hflux', 'eflux', 'rhoa', 'tsh'))
        ga('close 1')

        ga.open('http://goldsmr1.sci.gsfc.nasa.gov:80/dods/MAT3FVCHM') # 3D fields
        ga('set z 1 72')
        self.sampleFile(ga,onlyVars=('u','v','t','qv','delp'), npzFile=npzFile)
        
#---
    def getMERRAero(self, npzFile=None, dir='/home/adasilva/iesa/MERRAero'):
        """
        Get Met fields from MERRA.
        """
        from grads import GrADS

        ga = GrADS(Window=False,Echo=False)
    
        ga.open(dir+'/geosgcm_surf.ddf')
        self.sampleFile(ga,onlyVars=('ustar', 'bstar', 'pblh', 'shfx', 'lhfx', 'rhos', 'ts'))
        ga('close 1')

        ga.open(dir+'/inst3d_prog_v.ddf')
        ga('set z 1 72')
        self.sampleFile(ga,onlyVars=('u','v','t','qv'))
        ga('close 1')

        ga.open(dir+'/inst3d_aer_v.ddf')
        ga('set z 1 72')
        self.sampleFile(ga,onlyVars=('delp',), npzFile=npzFile)

#---
    def getFP(self, npzFile=None, dir='/home/adasilva/iesa/aerosol/experiments/seac4rs_01'):
        """
        Get Met fields from GEOS-5 Forward Processing.
        """
        from grads import GrADS

        ga = GrADS(Window=False,Echo=False)
    
        ga.open(dir+'/tavg1_2d_flx_Nx.ddf')
        self.sampleFile(ga,onlyVars=('ustar', 'bstar', 'pblh', 'hflux', 'eflux', 'rhoa', 'tsh'))
        ga('close 1')

        ga.open(dir+'/inst3_3d_asm_Nv.ddf')
        ga('set z 1 72')
        self.sampleFile(ga,onlyVars=('u','v','t','qv','o3','delp'),
                        npzFile=npzFile)
        ga('close 1')

#---
    def sampleFile(self, ga, npzFile=None, onlyVars=None, Verbose=True):
        """
        Interpolates all variables of defualt file in ga and optionally
        save them to file *npzFile*
        """
        
        if self.sample == None:
            self.sample = Sample() # first time
            
        if onlyVars is None:
            fh = ga.query('file',Quiet=True)
            onlyVars = fh.vars

        tymes = self.tyme
        lons = self.lon
        lats = self.lat

        # Loop over variables on file
        # ---------------------------
        for v in onlyVars:
            if Verbose:
                print("<> Sampling ", v)
            var = ga.sampleXYT(v,lons,lats,tymes,Verbose=Verbose)
            self.sample.__dict__[v] = var.data

        if npzFile is not None:
            savez(npzFile,**self.sample.__dict__)            
       
#---
    def sampleLoadz(self,npzFile):
        """
        Loads sample from npz file.
        """
        self.sample = Sample()
        npz = load(npzFile)
        for v in list(npz.keys()):
            self.sample.__dict__[v] = npz[v]

#---
    def drawMap(self,bbox=None,I=None):
        """
        Draw plume locations on map.
        """
        from mpl_toolkits.basemap import Basemap

        if bbox is None:
            bbox = self.lon.min()-5, self.lat.min()-5, self.lon.max()+5, self.lat.max()+5
        if I is None:
            I = list(range(self.N))

        # Basemap
        # -------
        m = Basemap(projection='cyl',
            llcrnrlon=bbox[0],urcrnrlon=bbox[2],
            llcrnrlat=bbox[1],urcrnrlat=bbox[3],
            rsphere=6371200.,resolution='l',area_thresh=10000)

        m.bluemarble()
        m.drawcoastlines()
        m.drawstates()
        m.drawcountries()

        # plot fires
        # ----------
        m.plot(self.lon[I],self.lat[I],'ro',label='Below PBL')
    
        return m

#....................................................................
class Sample(object):
    def __init__(self):
        return

#....................................................................

if __name__ == "__main__":

#    m = MINXs('/home/adasilva/workspace/misrPlumes/western_fires_2013/*.txt')
    m = MINXs('/Users/adasilva/workspace/misrPlumes/western_fires_2013/*.txt')

    #m.getFP(npzFile='seac4rs_01.npz')

def hold():


    #m = MINXs('/Users/adasilva/workspace.local/misrPlumes/canada2008/Plumes_O450*.txt')
    #m.getMERRA(npzFile='merra_O450.npz')

    m = MINXs('/Users/adasilva/workspace.local/misrPlumes/canada2008/Plumes_*.txt')

#    m.getMERRA(npzFile='merra.npz')
    
    aqua = CSV('AquaWindLocation.csv')
    terra = CSV('TerraWindLocation.csv')



