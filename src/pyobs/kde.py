"""
     Simple functions for plotting KDE.
"""

import os

from types             import *

from numpy             import sum, zeros, ones, sqrt, std, mean, unique, \
                              concatenate, where, linspace, meshgrid, exp, savez, \
                              percentile
from matplotlib        import cm, ticker

from matplotlib.pyplot import plot, title, xlabel, ylabel, figure, imshow, \
                              grid, colorbar, draw, axes, legend

from scipy             import stats, mgrid, c_, reshape, random, rot90

MISSING = 1.e15

#..............................................................
class aodFormat(ticker.Formatter):
    def __call__(self,x,pos=None):
        y = exp(x)-0.01
        return '%4.2f'%y

#..............................................................
def _cat (X, Y):
    """
    Given 2 arrays of same shape, returns array of shape (2,N),
    where N = X.size = Y.size
    """
    xy = concatenate((X.ravel(),Y.ravel())) # shape is (N+N)
    return reshape(xy,(2,X.size))         # shape is (2,N)

class aodFormat(ticker.Formatter):
    def __call__(self,x,pos=None):
        y = exp(x)-0.01
        return '%4.2f'%y

#---           
def calc_kde1d(Xin,N=256,range=None,Verbose=False,name=None):
    """
    Calculates 1D KDE. On input,

       Xin   ---  input data array
       N     ---  number of bins to evaluate KDE
       range ---  range of Xin values to work with

    Example:

       bins, P = calc_kde1d(obs,range=(-2,2))
       
    """

    try:
        X = Xin.data[Xin.mask==False].ravel()    
    except:
        X = Xin
    if range==None:
        q = (0.0, 25.0, 50.0, 75.0, 100.0)
        prc = prctile(X,q)
        range = [prc[0],prc[4]]
    bins = linspace(range[0],range[1],N)
    if Verbose:
        if name != None:
            print(name)
            print('Evaluating 1D kernel with %d observations'%len(X))
    kernel = stats.kde.gaussian_kde(X)
    if Verbose:
        print('Evaluating 1D KDE with %d bins'%N)
    P = kernel(bins)
    return (bins,P)

#---
def plot_kde1d(bins,P,Title=None,Xlabel=None,Ylabel=None):

    plot(bins,P)
    if Title != None:
        title(Title)
    if Xlabel != None:
        xlabel(Xlabel)
    if Ylabel != None:
        xlabel(Ylabel)

#---           
def calc_kde2d(x_values,y_values,x_range=None,y_range=None,
               Nx=256, Ny=256,
               npz=None, Verbose=True, name=None):

    if Verbose:
        if name != None:
            print("[] ", name)
        print("Starting the 2D kernel density estimation with %d data points..."\
              %len(x_values))

    kernel = stats.kde.gaussian_kde(_cat(x_values,y_values))

    if x_range==None:
        q = (0.0, 25.0, 50.0, 75.0, 100.0)
        prc = prctile(x_values,q)
        x_range = [prc[0],prc[4]]

    if y_range==None:
        y_range = x_range
        
    x_bins = linspace(x_range[0],x_range[1],Nx)
    y_bins = linspace(y_range[0],y_range[1],Ny)

    if Verbose:
        print("Evaluating 2D kernel on grid with (Nx,Ny)=(%d,%d) ..."%(Nx,Ny))

    X, Y = meshgrid(x_bins,y_bins) # each has shape (Ny,Nx)
    P = kernel(_cat(X,Y))           # shape is (Ny*Nx)
    P = reshape(P,X.shape)
        
    if Verbose:
        print("X, Y, P shapes: ", X.shape, Y.shape, P.shape)

    #   Save to file
    #   ------------
    if npz != None:
        print("Saving 2D KDE to file <"+npz+"> ...")
        savez(npz,pdf=P,x_bins=x_bins,y_bins=y_bins)
    
    return (x_bins,y_bins,P)

#---
def plot_kde2d( x_bins, y_bins, P, centroid=False, formatter=None,dpi=None,
                regression=None,Title=None,xLabel=None,yLabel=None):


    #   Plot results with 2 different colormaps
    #   ---------------------------------------
    fig = figure(dpi=dpi)
    ax = fig.add_axes([0.1,0.1,0.75,0.75])
    if formatter != None:
        ax.xaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)
    imshow(P, cmap=cm.gist_earth_r, origin='lower', 
           extent=(x_bins[0],x_bins[-1],y_bins[0],y_bins[-1]) )
    grid()
    plot([x_bins[0],x_bins[-1]],[y_bins[0],y_bins[-1]],'k')

    if Title != None:
        title(Title)
    if xLabel != None:
        xlabel(xLabel)
    if yLabel != None:
        ylabel(yLabel)

    # Regression
    # ----------
    if regression!=None:
        plot_linregress(regression[0],regression[1],x_bins)

    # Centroid
    # --------
    if centroid:
        X, Y = meshgrid(x_bins,y_bins) # each has shape (Ny,Nx)
        y_c = sum(P*Y,axis=0)/sum(P,axis=0)
        plot(x_bins,y_c,'r',label='Centroid')

    # Tighter colorbar
    # ----------------
    draw()
    bbox = ax.get_position()
    l,b,w,h = bbox.bounds # mpl >= 0.98
    cax = axes([l+w+0.02, b, 0.04, h]) # setup colorbar axes.
    colorbar(cax=cax) # draw colorbar

#---
def plot_linregress(x_values,y_values,x_bins):
    slope, intercept, r, prob, see = stats.linregress(x_values,y_values)
    yy = slope * x_bins + intercept
    plot(x_bins,yy,'k',label='Regression')
    legend(loc='upper left')

#---           
def print_stats(name,x=None):
    "Prints simple stats"
    if type(name) is not StringType:
        x = name
        name = 'mean,stdv,rms,min,25%,median,75%,max: '
    if name == '__header__':
        print('')
        n = (80 - len(x))/2
        print(n * ' ' + x)
        print(n * ' ' + len(x) * '-')
        print('')
        print('   Name       mean      stdv      rms      min     25%    median     75%      max')
        print(' ---------  -------  -------  -------  -------  -------  -------  -------  -------')
    elif name == '__sep__':
        print(' ---------  -------  -------  -------  -------  -------  -------  -------  -------')
    elif name == '__footer__':
        print(' ---------  -------  -------  -------  -------  -------  -------  -------  -------')
        print('')
    else:
        ave = x.mean()
        std = x.std()
        rms = sqrt(ave*ave+std*std)
        q = (0.0, 25.0, 50.0, 75.0, 100.0)
        prc = prctile(x,q)
        print('%10s  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  '%\
            (name,ave,std,rms,prc[0],prc[1],prc[2],prc[3],prc[4]))
#---           

def writeg_del(obs, ana, bkg, lon, lat, vname, units,
               filename, title, time,
               refine=8,res=None,Verbose=True):
       """
         Grid & Write Obs, Ana and Bkg to file.

         vname   -- variable name for long name, e.g., 'AOD'
         units   -- variable unit

         filename - output file name
         title   -- title for file
         time    -- time stamp for file
         
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
       from binObs_  import binobs2d, binrms2d


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

       nymd = 10000 * time.year + 100 * time.month  + time.day
       nhms = 10000 * time.hour + 100 * time.minute + time.hour

       Vtitle = [ vname + ' Observations (mean)',
                  vname + ' Analysis (mean)',
                  vname + ' Background (mean)',
                  vname + ' Observations (RMS)',
                  vname + ' Analysis (RMS)',
                  vname + ' Background (RMS)',
                ]

       Vname  = ['obs','ana', 'bkg',
                 'obs_rms','ana_rms', 'bkg_rms']
       
       Vunits = 6 * [ units, ]
       Kmvar  = 6 * [ 0, ]
       
       source = 'NASA/GSFC/GMAO GEOS-5 Aerosol Group'
       contact = 'arlindo.dasilva@nasa.gov'

#      Create the file
#      ---------------
       f = GFIO()
       f.create(filename, Vname, nymd, nhms,
                lon=glon, lat=glat, levs=[1,], levunits='1',
                vtitle=Vtitle, vunits=Vunits,kmvar=Kmvar,amiss=MISSING,
                title=title, source=source, contact=contact)

#      Grid variable and write to file
#      -------------------------------
       f.write('obs',     nymd, nhms, binobs2d(lon,lat,obs,im,jm,MISSING) )
       f.write('ana',     nymd, nhms, binobs2d(lon,lat,ana,im,jm,MISSING) )
       f.write('bkg',     nymd, nhms, binobs2d(lon,lat,bkg,im,jm,MISSING) )
       f.write('obs_rms', nymd, nhms, binrms2d(lon,lat,obs,im,jm,MISSING) )
       f.write('ana_rms', nymd, nhms, binrms2d(lon,lat,ana,im,jm,MISSING) )
       f.write('bkg_rms', nymd, nhms, binrms2d(lon,lat,bkg,im,jm,MISSING) )

       if Verbose:
           print("[w] Wrote file "+filename)
           

