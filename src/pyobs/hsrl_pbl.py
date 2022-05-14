#!/usr/bin/env python

"""
Simple class to read HSRL PBL ASCII files.

"""

from glob  import glob
from numpy import array, sort
from datetime import datetime, timedelta

class HSRL_PBL(object):

    def __init__ (self, filename,ARCTAS=False):
        """
        Given a filename (or a file name patter such as *.txt),
        read a text file where first row has name of columns,
        first colume has time in hours withinna day, and return
        an object with each column as an attribute in the form
        of a numpy array.
        """

        Files = sort(glob(filename))

        # Empty list for each column
        # --------------------------
        Header = open(Files[0]).readline().split()
        for name in Header:
            self.__dict__[name] = []

        # Read and parse files
        # --------------------
        for f in Files:
            cdate = f.split('_')
            if len(cdate[2])<len(cdate[1]):
                cdate = cdate[1]
            else:
                cdate = cdate[2]
            day = datetime(int(cdate[0:4]),int(cdate[4:6]),int(cdate[6:]))
            print("[ ] working on ", day)
            Lines = open(f).readlines()
            for line in Lines[1:]:
                Line = line.replace('\r\n','').split()
                time = day + timedelta(seconds=int(float(Line[0])*3600))
                self.__dict__['Time'].append(time)
                i = 1
                for name in Header[1:]:
                    self.__dict__[name].append(float(Line[i]))
                    i += 1
                         
        # Create arrays
        # -------------
        for name in Header:
            try:
                self.__dict__[name] = array(self.__dict__[name])
            except:
                print("<> Failed creating array for ", name)

    def addVar(self,ga,expr='pblh',vname=None,Verbose=False):

        """
        Given a grads object having the correct file as default,
        creates an attribute with *expr* interpolated to
        (time,lat,lon).

        """

        if vname == None:
            vname = expr

        self.__dict__[vname] = ga.sampleXYT(expr,
                                            self.Longitude,
                                            self.Latitude,
                                            self.Time,
                                            Verbose=Verbose)
