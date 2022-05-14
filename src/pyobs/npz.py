"""
Simple class for reading NPZ files.
"""

from __future__ import print_function

from numpy import load, savez, concatenate, shape, sort
from glob  import glob

MISSING = -999.

class NPZ(object):

    def __init__ (self,npzFiles,Verbose=False):
        """
        Loads npz files, returning variables as attributes.
        On input npzFiles can contain a single file, a list of file,
        or a string with wild characters like '*.npz'.
        """

        if isinstance(npzFiles,str):
            npzFiles = sorted(glob(npzFiles))

        if Verbose:
            print('[] Loading ', npzFiles[0])

        # Single file
        # -----------
        if len(npzFiles) == 1:
            f = load(npzFiles[0])
            for v in f.keys():
                self.__dict__[v] = f[v]

        # Multiple files
        # --------------
        else:
            
            # For many files, process first file in list
            # ------------------------------------------
            f = load(npzFiles[0])
            V = dict()
            for v in f.keys():
                if len(shape(f[v])) == 0: 
                    V[v] = [[f[v],],]
                else:
                    V[v] = [f[v],]

            # Append the other files
            # ----------------------
            for npzFile in npzFiles[1:]:
                if Verbose:
                    print( '[] Loading ', npzFile)
                f = load(npzFile)
                for v in V:
                    if len(shape(f[v])) == 0: 
                        V[v].append([f[v],])
                    else:
                        V[v].append(f[v])
                f.close()

            # Concatenate
            # -----------
            for v in V:
                self.__dict__[v] = concatenate(V[v])

    def __getitem__(self,name):
        return self.__dict__[name]
        
    def sampleVar(self,ctlfile,vname,aname=None,I=None,**kwds):
        """
        Interpolates variable *vname* from a gridded GFIO collection to obs location/time.
        On input, *aname* is the attribute name to be defined; if not specified it defaults
        to *vname*.
        """
        from gfio import GFIOctl
        f = GFIOctl(ctlfile)
        if aname == None:
            aname = vname
        if I is None:
            self.__dict__[aname] = f.sample(vname,self.lon,self.lat,self.time,**kwds)
        else:
            self.__dict__[aname][I] = MISSING * ones(len(self.lon))
            self.__dict__[aname][I] = f.sample(vname,self.lon[I],self.lat[I],self.time[I],**kwds)
            
    def savez(self,outFile):
        """
        Save all atributes to npz file.
        """
        savez(outFile,**self.__dict__)
        
