"""

Very minimalistic ODS reader.

"""

VARS = ['lon', 'lat', 'lev', 'time', 'obs', 'qcexcl','kt']
EXTRA = ['kt', 'kx', 'ks', 'qchist', 'xm', 'xvec', 'omf', 'oma']

from netCDF4 import Dataset

from numpy import array, exp
from datetime import timedelta
from dateutil.parser import parse as isoparser

class ODSreader(object):

    def __init__ (self, filename, Extra=None):
        """
        Read select attributes from an ODS file.
        """ 

        nc = Dataset(filename)

        nc.set_auto_mask(False)
        nc.set_auto_scale(False)
        
        Vars = VARS
        if Extra is not None:
             Vars += Extra

        self._getVars(nc,Vars)

        self.nobs = len(self.lon)

        self.kx_names = []
        for x in nc.variables['kx_names'][:]:
            try:
                self.kx_names += [bytearray(x).decode().rstrip(),]
            except:
                break
        self.kt_names = []
        for x in nc.variables['kt_names'][:]:
            try:
                self.kt_names += [bytearray(x).decode().rstrip(),]
            except:
                break
                
        # Create python time
        # ------------------
        tunits = nc.variables['time'].units
        onemin = timedelta(minutes=1)
        toff  = isoparser(tunits.split(' since ')[1].replace(' ','T'))
        t = array([toff + self.time[i] * onemin for i in range(self.nobs)])

        self.time = t # array of python datetime objects

    def _getVar(self,nc,name):
        """
        Retrieve and scale one variable.
        """
        
        v = nc.variables[name]
        q = v[:,:]

        try:
            q = q * v.scale_factor
        except:
            pass
        try:
            q = q + v.add_offset
        except:
            pass

        return array(q[:])
        
    def _getVars(self,nc,Vars):
        """
        Read and flatten variables on file. 
        Assuems 1 synoptic time per file.
        """

        # Use lon for mask
        # ----------------
        lon = self._getVar(nc,'lon')
        I = (lon>=-180)&(lon<=180)
        self.lon = lon[I]

        Vars.remove('lon')
        for name in Vars:
            q = self._getVar(nc,name)
            self.__dict__[name] = q[I]
            
#-------------

if __name__ == "__main__":

    m = ODSreader('d5124_m2_jan10.aod.obs.20161231_2100z.ods',Extra=('kx','kt','omf','oma'))

    # Convert from log(aod+0.1) to aod
    # --------------------------------
    aod_o = exp(m.obs[:]) - 0.01              # observed AOD
    aod_f = exp(m.obs[:]-m.omf[:]) - 0.01     # forecast AOD: F = O - (O-F)
    aod_a = exp(m.obs[:]-m.oma[:]) - 0.01     # analysis AOD: A = O - (O-A)
    
    # Print first 100 values
    # Note: kx is data source index, kt is data type index
    # ----------------------------------------------------
    for i in range(100):
        print(("%20s %20s %6.2f %6.2f %s %7.4f %7.4f %7.4f"%\
        (m.kx_names[m.kx[i]-1],m.kt_names[m.kt[i]-1],m.lon[i],m.lat[i],m.time[i].isoformat(),aod_o[i],aod_f[i],aod_a[i])))
