"""
   Classes for reading CSV files the Maritime Aerosol Network (MAN) data.
"""

MISSING = -999.
from datetime import datetime
from numpy import loadtxt, ones, savez, pi, log

Months = ('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')

def date2nymd(s):
    if ('/' in s):
        S = s.split('/')
        yy = 2000+int(S[2])
        mm = int(S[0])
        dd = int(S[1])
    elif (':' in s):
        S = s.split(':')
        yy = int(S[2])
        mm = int(S[1])
        dd = int(S[0])

    return (yy,mm,dd)

def time2nhms(s):
    S = s.split(':')
    h = int(S[0])
    m = int(S[1])
    return (h,m)

def gatime(date,time):
    yy,mm,dd = date2nymd(date)
    h, m = time2nhms(time)
    return "%d:%dZ%d%s%d"%(h,m,dd,Months[mm-1],yy)

def _getCols(fname):
    """
    From CSV Spreadsheet.
    """
    f = open(fname,'r')
    cols = f.readline()
    f.close()
    #cols = 'Date,Time,Air Mass,Latitude,Longitude,tau_340,tau_380,tau_440,tau_500,tau_675,tau_870,tau_1020,tau_1640,Water Vapor(cm),440-870_Angstrom_Exponent,STD_Latitude,STD_Longitude,STD_340,STD_380,STD_440,STD_500,STD_675,STD_870,STD_1020,STD_1640,STD_Water_Vapor(cm),STD_440-870_Angstrom_Exponent,Number_of_Observations,Last_Processing_Date'
    return cols

def _getNames(cols):
    
    Names = cols.split(',')
    NewName = {}
    for name in Names:
        NewName[name] = name # by default same name as Rob's

    NewName['Longitude'] = 'lon'
    NewName['Latitude'] = 'lat'

    NewName['AOD_340nm'] = 'aTau340'
    NewName['AOD_380nm'] = 'aTau380'
    NewName['AOD_440nm'] = 'aTau440'
    NewName['AOD_500nm'] = 'aTau500'
    NewName['AOD_675nm'] = 'aTau675'
    NewName['AOD_870nm'] = 'aTau870'
    NewName['AOD_1020nm'] = 'aTau1020'
    NewName['AOD_1640nm'] = 'aTau1640'

    NewName['Date(dd:mm:yyyy)'] = 'Date'
    NewName['Time(hh:mm:ss)']   = 'Time'

    return (Names, NewName)

#----
Vars = ( 'Longitude', 'Latitude', 'Date', 'Time',
         'AOD_340nm', 'AOD_380nm', 'AOD_440nm', 'AOD_500nm',
         'AOD_675nm', 'AOD_870nm', 'AOD_1020nm', 'AOD_1640nm' )

class MAN(object):
    """Base class for Martime Aerosol Network (MAN)."""

    def __init__ (self,
                  fname='All_MAN_series_Level2_June21_2010.csv',
                  Vars=Vars):

        cols = _getCols(fname)
        Names, NewName = _getNames(cols)

        self.filename = fname
            
        iVars = ()
        formats = ()
        converters = {}
        for name in Vars:
            try:
                i = Names.index(name)
            except:
                raise ValueError("cannot find <%s> in file"%name)

            iVars = iVars + (i,)
            if name=='Date':
                formats = formats + ('S8',)
            elif name=='Date(dd:mm:yyyy)':
                formats = formats + ('S10',)
            elif 'Time' in name:
                formats = formats + ('S8',)
            else:
                converters[i] = lambda s: float(s or MISSING)
                formats = formats + ('f4',)

#       Read the data
#       -------------
        data = loadtxt(fname, delimiter=',',
                       dtype={'names':Vars,'formats':formats},
                       converters = converters,
                       skiprows=1, usecols=iVars)
        N = len(data)

        self.N = N

#       Save data columns as attributes, with nicer names
#       -------------------------------------------------
        for i in range(len(Vars)):
            if formats[i]=='f4':
                v = ones(N)
                for j in range(N):
                    v[j] = data[j][i]
            else:
                v = []
                for j in range(N):
                    v.append(data[j][i])

            self.__dict__[Vars[i]] = v

            if NewName[Vars[i]] != Vars[i]:
                self.__dict__[NewName[Vars[i]]] = v # alias

#       Interpolate AOD to 553
#       ----------------------
        alpha = (log(553.) - log(500.) ) / ( log(675.) - log(500.) )
        i = (self.aTau500>0) & (self.aTau675>0)
        self.aTau550 =  MISSING * ones(N)
        self.aTau550[i] =  (1.-alpha) * self.aTau500[i] + alpha * self.aTau675[i]

#       Create grads time
#       -----------------
        self.time = []
        for i in range(N):
               self.time.append(gatime(self.Date[i],self.Time[i]))

#       Create python datetime
#       ----------------------
        self.tyme = []
        for i in range(N):
            yy,mm,dd = date2nymd(self.Date[i])
            h, m = time2nhms(self.Time[i])
            self.tyme.append(datetime(yy,mm,dd,h,m))
#---
    def addVar(self,ga,outfile=None,expr='tau',vname=None, clmYear=None):
        """
        Given a grads object having the correct file as default,
        writes out a CSV file with the specified expression/variable. When *clmYear* is
        specified, the actual year in the time attribute is replaced
        with he climatological year *clmYear*. 
        """

        N = self.N
        U = ones(N)
        U[:] = MISSING

        if vname == None:
            vname = expr

        for i in range(N):

            x = self.lon[i]
            y = self.lat[i]

            if clmYear == None:
                t = self.time[i]
            else:
                t = self.time[i][:-4] + str(clmYear) # replace year

            ga('set lon '+str(x))
            ga('set lat '+str(y))
            ga('set time %s'%t,Quiet=True)

            U[i] = ga.expr(expr).data # nearest neighbor

            if i%100 == 0:
                print("%9s %9s %8.3f %8.3f %10.3f  ...%8.3f%%"\
                  %(self.Date[i],self.Time[i],x,y,U[i],100.*i/float(N)))

        self.__dict__[vname] = U

        if outfile is not None:
            version = 1
            meta = [ version, self.filename, vname, expr ]
            savez(outfile,meta=meta,lon=self.lon,lat=self.lat,
                  date=self.Date,time=self.Time,var=U)

#....................................................................

def _estimate():

    from grads import GrADS
    
    man = MAN()
    
    ga = GrADS(Echo=0,Window=False)

#    fh = ga.open('nnr_001.modo.ddf')
#    man.addVar(ga,'nnr_001.modo.tau.npz',expr='tau')
#    man.addVar(ga,'nnr_001.modo.tau_.npz',expr='tau_')

    ga('reinit')
    fh = ga.open('nnr_001.mydo.ddf')
    man.addVar(ga,'nnr_001.mydo.tau.npz',expr='tau')
    man.addVar(ga,'nnr_001.mydo.tau_.npz',expr='tau_')

if __name__ == "__main__":

    from pylab import plot, xlabel, ylabel, title, savefig, legend
    from numpy import linspace, load
    from scipy.stats import linregress

    ident = 'mydo'

    if ident == 'modo':
        prod = 'Terra'
    else:
        prod = 'Aqua'

    man = MAN()
    mxdo  = load('nnr_001.'+ident+'.tau.npz')['var']
    mxdo_  = load('nnr_001.'+ident+'.tau_.npz')['var']
    anet = man.aTau550

    i = (anet>0) & (mxdo>0) & (mxdo_>0)

    lmxdo = log(mxdo[i]+0.01)
    lmxdo_ = log(mxdo_[i]+0.01) 
    lanet = log(anet[i]+0.01)

    slope, intercept, r_value, p_value, std_err = linregress(lanet,lmxdo)
    slope_, intercept_, r_value_, p_value_, std_err_ = linregress(lanet,lmxdo_)

    x = linspace(-5.,1.,100)
    y = slope * x + intercept
    y_ = slope_ * x + intercept_

    plot(lanet,lmxdo_,'ro')
    plot(lanet,lmxdo, 'bo')
    plot(x,y,'b',linewidth=3,label=prod)
    plot(x,y_,'r',linewidth=3,label='Neural Net')
    plot(x,x,'k',linewidth=3,label='1:1')
    xlabel('Maritime Aerosol Network')
    ylabel('Retrievals')
    title(r'Log($\tau_{550}$+0.01)')
    legend(loc='upper left')
    savefig('man.'+ident+'.scat.png')
    savefig('man.'+ident+'.scat.eps')
