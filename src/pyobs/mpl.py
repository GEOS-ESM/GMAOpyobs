#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script reads MPLnet data (Level 1.5) into memory and displays a curtain of bext. 
See example of use at the bottom:

A NetCDFData class to store the data and attributes

A read_netcdf function that:
    Takes a filename (required) and
    (optional) list of variables to read (default variables if none specified)
    Reads data & attributes for each variable
    Stores everything in the class structure
    Verb=0(default), no print,=1 variables names and sizes, =2 var names and attributes

Main Subroutines:
    read_file : actual read a single file
    get_path  : set path where files are stored (needed if working in different clusters)
    f_mpl_plot2d : plot data
----------------------------------------
Versioning:
Mar/24/2025 - Delivered to Github
Mar/10/2025 - based on v2, this code implements the actual AOCH calculation because right now it only computes AOD
Mar/7/2025 - Based pn     cbar = plt.colorbar(extend='max') , it now selects subsets by time range set by the user and then takes averages and computes AOCH
Feb/24/2025 - Same as v2 but the plotting routine is now inside a module. it Works!
Feb/19/2025 - Ver2: same as V1 with the addition of simple plot creation

Created on Fri Jan 17 10:51:10 2025

@author: sgassoumd + chat.gsfc
"""
from netCDF4 import Dataset
import numpy as np
import platform
import sys
from datetime import datetime, timedelta # need to get date from filename

# import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors

#from matplotlib.colors import ListedColormap
from scipy import ndimage


class MPL_L15:
    def __init__(self):
        """Initialize the class with empty attributes"""
        self.time = None
        self.time_attributes = {}

        self.altitude = None
        self.altitude_attributes = {}

        self.extinction = None
        self.extinction_attributes = {}

        self.extinction_err = None
        self.extinction_err_attributes = {}

        self.qa_extinction = None
        self.qa_extinction_attributes = {}

    @classmethod
    def read_file(cls, filename, variables=None, verb=0):
        """
        Read data from NetCDF file

        Parameters:
        filename (str)  : Path to the NetCDF file
        variables (list): List of variables to read. If None, reads default variables
        verb (int): Verbosity level
                   0: No printing
                   1: Print variable name, type, and dimensions
                   2: Print level 1 info plus attributes

        Returns:
        NetCDFData: Class instance containing the requested data and attributes
        """
        # Default variables to read if none specified
        default_variables = ['time', 'altitude', 'extinction', 'extinction_err']

        # Use default variables if none provided
        variables_to_read = variables if variables is not None else default_variables

        # Initialize the class
        instance = cls()

        try:
            # Open the NetCDF file
            with Dataset(filename, 'r') as nc:
                # Read each requested variable
                for var_name in variables_to_read:
                    if var_name in nc.variables:
                        # Get the variable
                        var = nc.variables[var_name]

                        # Read the data
                        var_data = var[:]

                        # Read all attributes
                        var_attributes = {attr: var.getncattr(attr)
                                       for attr in var.ncattrs()}

                        # Print information based on verbosity level
                        if verb >= 1:
                            print("\nVariable Information:")
                            print(f"Name: {var_name}")
                            print(f"Type: {var.dtype}")
                            print(f"Dimensions: {var.dimensions}")
                            print(f"Shape: {var.shape}")

                            if verb >= 2:
                                print("\nAttributes:")
                                for attr, value in var_attributes.items():
                                    print(f"  {attr}: {value}")
                            print("=" * 50)

                        # Store data and attributes in the class
                        setattr(instance, var_name, var_data)
                        setattr(instance, f"{var_name}_attributes", var_attributes)
                    else:
                        print(f"Warning: Variable {var_name} not found in file")

        except Exception as e:
            print(f"Error reading NetCDF file: {str(e)}")
            return None

        # convert time to datetime
        dt_Offset = 2400000.500 # Julian date of 1858-11-17
        instance.tyme = np.array([datetime(1858, 11, 17) + timedelta(julian_date-dt_Offset) for julian_date in instance.time])

        return instance

    def get_variable_names(self):
        """Return list of available variables"""
        return [attr for attr in self.__dict__.keys() if not attr.endswith('_attributes')]

    def get_variable_info(self, variable_name):
        """Print information about a specific variable"""
        if hasattr(self, variable_name):
            data = getattr(self, variable_name)
            attrs = getattr(self, f"{variable_name}_attributes", {})
            print(f"\nVariable: {variable_name}")
            print(f"Shape: {data.shape}")
            print(f"Type: {data.dtype}")
            print("Attributes:", attrs)
        else:
            print(f"Variable {variable_name} not found")

def get_path(now_os, now_computer):
    # Your existing get_path function remains the same
    if now_os=='win32': # PC
        base_path    = 'C:/Users/sgasso/Downloads/'
        pth_fig_out='C:/Users/sgasso/OneDrive - NASA/ToShare/2025/GEOS/Pyfigures/'
    elif now_os == 'darwin': # My Laptop
#        base_path='/Users/sgasso/Downloads/'
        base_path='/Volumes/ExtData1/Surface/MPL/V3_L15/'
        pth_fig_out='/Users/sgasso/Library/CloudStorage/OneDrive-NASA/ToShare/2025/AOCH/PyFigures/'
    elif now_os == 'linux' and "calculon" in now_computer:
        base_path = '/nobackup/CALIPSO/Level1.5/Santiago/'
    elif now_os == 'linux' and "discover" in now_computer:
        base_path = '/nobackup/CALIPSO/Level1.5/Santiago/'
    else:
        print('Current operating system no recognized.')
        print('Cannot set path to MPL  files. Terminate Run')
        sys.exit()

    return base_path,pth_fig_out

#### Read list of MPL sites
def read_mpl_site_list(filename):
#Expected format: comma separated
#line 1 :Name_MPL_sites , name,lat, lon, altiude(km)
#line 2 : GSFC           ,38.9930,-76.8400,0.050
#.............
#
    sites = []
    with open(filename, 'r') as file:
        header = file.readline()
        for line in file:
            data = [item.strip() for item in line.split(',')]
            site = (
                data[0],
                float(data[1]),
                float(data[2]),
                float(data[3])
            )
            sites.append(site)
    return sites #this is a tuple


def f_mpl_plot2d(param, xtime, z, title_strings, colorange=(1e-6, 1e-3), scale='lin', window_ave=0, esY=10):
    """
    Create 2D plot of MPL data

    Parameters:
    param (array): Data to plot (extinction)
    xtime (array): Time array
    z (array): Height array
    title_strings (list): [title_str, colbar_str]
    colorange (tuple): Min and max values for colorbar
    scale (str): Scale type ('lin' or 'log')
    window_ave (int): Window size for averaging
    esY (int): Spacing for time axis
    """
    # import matplotlib.pyplot as plt
    # import matplotlib.colors as colors
    # import numpy as np
    # from scipy import ndimage

    plot_title = title_strings[0]
    colorbar_title = title_strings[1]
    esX = 1  # Fixed X-axis spacing for height dimension

    # Data preparation
    lastZ = len(param[0,:])
    var0 = param[::esX,::esY]
    plot_time = xtime[::esY]
    plot_z = z[::-esX]

    # Apply filtering
    if window_ave == 0:
        var0[var0 < 1e-12] = np.nan
        var = var0
    else:
        var = ndimage.uniform_filter(var0, size=window_ave, mode='constant')
        var[var <= 0] = np.nan

    # Prepare plot variables
    plot_var = np.log10(var)
    vmin = np.log10(colorange[0])
    vmax = np.log10(colorange[1])
    space_col = 30

    # Color setup
    base_cmap = plt.get_cmap('jet')
    levels = np.linspace(vmin, vmax, space_col + 1)
    calipso_colors = [base_cmap(i) for i in np.linspace(0, 1, space_col-10)]
    gray_colors = [colors.to_rgba(str(min(i/10, 1))) for i in range(1, 12)]
    calipso_colors.extend(gray_colors)
    cmap = colors.ListedColormap(calipso_colors)
    norm = colors.BoundaryNorm(levels, len(calipso_colors))

    # Create plot
    fig, ax = plt.subplots(figsize=(16, 6))
    plt.pcolor(plot_var, norm=norm, cmap=cmap)
    cbar = plt.colorbar(extend='max')


    # Colorbar setup for log scale
    if scale == 'log':
        bar_ticks = np.array([1e-6, 2e-6, 3e-6, 5e-6, 7e-6,
                             1e-5, 2e-5, 3e-5, 5e-5, 7e-5,
                             1e-4, 2e-4, 3e-4, 5e-4, 7e-4,
                             1e-3, 2e-3, 3e-3, 5e-3, 7e-3, 1e-2])
        log_bar = np.log10(bar_ticks)
        log_bar_with_ends = np.concatenate(([log_bar[0]], log_bar, [log_bar[-1]]))
        cbar.ax.minorticks_off()
        cbar.set_ticks(log_bar_with_ends)
        tick_labels = [f"{10**label:.1e}" for label in log_bar_with_ends]
        cbar.set_ticklabels(tick_labels)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label(colorbar_title, rotation=270, labelpad=15)

    # Axis setup

    yt = np.arange(len(plot_z))
    if yt[0] != 0:
        yt = np.insert(yt, 0, 0)
    if yt[-1] != len(plot_z) - 1:
        yt = np.append(yt, len(plot_z) - 1)
    nYTicks = 20

    iyt = np.round(np.linspace(0, len(yt) - 1, nYTicks)).astype(int)

    # X axis labels
    # x_h_start=plot_time[0].hour
    # x_h_end  =plot_time[-1].hour + 1

#    xt = np.arange(x_h_start, x_h_end) * 60
    xt = np.arange(0,len(xtime),10) # every X min
    xvals = []
    xstrs = []
    for i in xt:
        xvals.append(i)
        # xstrs.append('%s' % int(i/60))
        xstrs.append(xtime[i].strftime("%H:%M"))

    # Set axis labels and ticks
    plt.xticks(xvals, xstrs, fontsize=8)
    plt.yticks(iyt[::-1], [f"{c:.1f}" for c in plot_z[iyt]], fontsize=8)

    # Final labels and title
    plt.ylabel('Height (km)')
    plt.xlabel('Time (UTC)')
    if window_ave == 0:
        plt.title(plot_title)
    else:
        plt.title(f"{plot_title} , Moving average = {window_ave}")

    return fig, ax

### Method 1: Average over time, then integrate over altitude
def method1(varin, varin_err, z_in):
    Nlayers, Ntimes = varin.shape

    # Average over time (second dimension)
    valid_data = (varin > 0) & (varin < 1e+20)
    avg_ext_in_time = np.zeros(Nlayers)
    avg_ext_in_time_err = np.zeros(Nlayers)
    valid_counts = np.sum(valid_data, axis=1)
    # Loop over each layer and compute time average ext(i_layer)
    for i in range(Nlayers):
        valid_indices = valid_data[i, :]
        if np.any(valid_indices):
            avg_ext_in_time[i] = np.mean(varin[i, valid_indices])
            avg_ext_in_time_err[i] = np.sqrt(np.sum(varin_err[i, valid_indices]**2)) / valid_counts[i]

    # Integrate over altitude (first dimension)
    valid_altitude = (avg_ext_in_time > 0)
    if np.sum(valid_altitude) > 1:
        AOCH = np.sum(avg_ext_in_time[valid_altitude]*z_in[valid_altitude])/np.sum(avg_ext_in_time[valid_altitude])
        ##### Error in AOCH ,
        ##### propagating error for ratio sum(ext_ave_layer*zi)/sum(ext_ave_layer)
        a= np.sum((avg_ext_in_time_err[valid_altitude]*z_in[valid_altitude])**2)
        b= np.sum(avg_ext_in_time_err[valid_altitude]**2)
        c= np.sum(avg_ext_in_time[valid_altitude])
        d= AOCH * c
        AOCH_err = AOCH * np.sqrt((a/(d**2))+(b/(c**2)))
    else:
        AOCH = np.nan
        AOCH_err = np.nan

    total_obs = Ntimes
    valid_obs = np.sum(valid_data, axis=0)
    obs_used = np.sum(valid_altitude)

    return float(AOCH), float(AOCH_err), total_obs, valid_obs, obs_used

# Method 2: Integrate over altitude, then average over time
def method2(varin, varin_err, z_in):
    Nlayers, Ntimes = varin.shape

    integrated_alt = np.zeros(Ntimes)
    integrated_alt_err = np.zeros(Ntimes)
    obs_used_in_integral = 0

    ### Loop over time(second dimension) and compute AOCH in each time step
    for i in range(Ntimes):
        valid_data = (varin[:, i] > 0) & (varin[:, i] < 1e+20)
        if np.sum(valid_data) > 1:
            integrated_alt[i] = np.sum(varin[valid_data, i]* z_in[valid_data])/np.sum(varin[valid_data, i])
            ##### Error in AOCH ,
            ##### propagating error for ratio sum(ext_ave_layer*zi)/sum(ext_ave_layer)
            a  = np.sum((varin_err[valid_data, i]* z_in[valid_data])**2)
            b  = np.sum(varin_err[valid_data, i]**2)
            c  = np.sum(varin[valid_data, i])
            d  = integrated_alt[i]  * c
            integrated_alt_err[i] = integrated_alt[i]  * np.sqrt((a/(d**2))+(b/(c**2)))
            #integrated_alt_err[i] = np.sqrt(np.sum((varin_err[valid_data, i]* z_in[valid_data])**2))/np.sum(z_in[valid_data])
            obs_used_in_integral += 1
        else:
            integrated_alt[i] = np.nan
            integrated_alt_err[i] = np.nan

    ##### Now average over time all AOCH computed at each time
    valid_integrated = ~np.isnan(integrated_alt)
    if np.any(valid_integrated):
        avg_integrated    = np.mean(integrated_alt[valid_integrated])
        avg_integrated_err = np.sqrt(np.sum(integrated_alt_err[valid_integrated]**2)) / np.sqrt(np.sum(valid_integrated))
    else:
        avg_integrated = np.nan
        avg_integrated_err = np.nan

    total_obs = Ntimes
    valid_obs = np.sum((varin > 0) & (varin < 1e+20), axis=0)

    return float(avg_integrated), float(avg_integrated_err), total_obs, valid_obs, obs_used_in_integral


######----------------------------------------------


if __name__ == "__main__":
    # Option 1: Direct call to the class method
    current_os    = platform.system()
    computer_name = platform.node()
    current_pth,pth_fig = get_path(current_os.lower(),computer_name)
    filename = "MPLNET_V3_L15_AER_20230820_MPL44258_GSFC.nc4"
    full_pathname = current_pth + filename

    # Define variables to read
    variables_to_read = ['time', 'altitude', 'extinction','extinction_err']

    ##### Set time frame of satellite overpass
    # User-provided start and end times
    t_start0 = 14 # start time in UTC
    t_end0   = 15 # End time in UTC

    # Read the data
    data = MPL_L15.read_file(full_pathname, variables=variables_to_read, verb=0)

    # ### now print some misc data directly
    # print(data.extinction.shape)
    # # Access attributes
    # print(data.time_attributes)
    # print(data.extinction_attributes.get('units'))

    #------------------------------------------------
    ### now do a quick plot along the track
    ###
    ### Get date from file
    # time is in the Julian Date format of number of days since January 1st, 4713 BC , noon UTC
    # See http://quasar.as.utexas.edu/BillInfo/JulianDatesG.html
    # and online calculator https://ssd.jpl.nasa.gov/tools/jdc/#/cd
    # where integer par is Ndays since this date and decimal time is hh:mm:ss since noon UTC
    # so for example : 2440588.000000 is A.D. 1970 January 1	12:00:00.0
    # this number can be substracted from data.time and use it as new reference point
    # altenrnativel the python module astropy can deal with it. NOTE the python module
    # datetime does not deal with this reference point .
    # I used a different approach based on the date in the filename.

    x     =data.time
    yyyy,mm,dd=[filename[18:22],filename[22:24],filename[24:26]]
    base_date = datetime(int(yyyy),int(mm),int(dd),0,0,0) ## Your base date
    xtime = []
    for jd in x:
        # Get just the fractional part since we know the date
        day_fraction = jd % 1
        if day_fraction >= 0.5:
            day_fraction = day_fraction - 0.5
        else:
            day_fraction = day_fraction+0.5
        # Convert fractional day to minutes
        #seconds = int(day_fraction * 86400)
        minutes = int(day_fraction * 1440)
        # Create datetime by adding the time to base date
        dt = base_date + timedelta(minutes=minutes)
        xtime.append(dt)
    ### Because MPL data is reported every 60 min,
    ### one could create a time array of size 60*24 = 1440 and it should be the same.
    # Convert the list to a numpy array
    xtime = np.array(xtime)


    ### Assign inputs to module
    site_name = filename[36:-4]
    z_in    =data.altitude[0,:]*1000 # convert to meters

    ### Now select variable of interest and transpose so it is compatible
    ### with plotting routine
    varin       =1e-3*data.extinction[0,:].T # 1440 x 400 to 400 x 1440
    varin_err   =1e-3*data.extinction_err[0,:].T # 1440 x 400 to 400 x 1440
    strings_plot= ['MPL Extinction 532nm, Site ' + site_name,'1/m']


    t_start = datetime(int(yyyy), int(mm), int(dd), t_start0, 0)
    t_end   = datetime(int(yyyy), int(mm), int(dd),   t_end0, 0)

    # Find the indexes where the values span the time frame
    indexes = np.where((xtime >= t_start) & (xtime <= t_end))[0]

    # Extract the respective values using the indexes
    var_sub = varin[:,indexes]
    var_sub_err=varin_err[:,indexes]

    ### now proceed to compute AOCH
    # print("Indexes:", indexes)


    # Method 1 : average over time each layer and the integrate over column
    # Method 2 : Integrate each column and then take average AOCHs
    # Apply both methods
    # result1, error1 = method1(varin, varin_err, z_in)
    # result2, error2 = method2(varin, varin_err, z_in)

    # Apply both methods
    result1, error1, total_obs1, valid_obs1, used_in_integral1 = method1(var_sub, var_sub_err, z_in)
    result2, error2, total_obs2, valid_obs2, used_in_integral2 = method2(var_sub, var_sub_err, z_in)

    # Print results for Method 1
    print("Method 1 : Average Bext in each layer and then compute AOCH")
    if np.isnan(result1):
        print("  AOCH (km): No valid data for integration")
    else:
        print(f"  AOCH (km): {result1/1000:.2f} ± {error1/1000:.2f}")
    print(f"  Total Profiles: {total_obs1}")
    print(f"  Stats of N layers with EXT obs: min={valid_obs1.min()}, max={valid_obs1.max()}, mean={valid_obs1.mean():.2f}")
    print(f"  N layers(out of 400) used in final mean: {used_in_integral1}")

    print("\nMethod 2: Compute AOCH in each profiles and then take mean:")
    if np.isnan(result2):
        print("  AOCH (km): No valid data for integration")
    else:
        print(f"  AOCH (km): {result2/1000:.2f} ± {error2/1000:.2f}")
    print(f"  Total Profile: {total_obs2}")
    print(f"  Stats of N layers w/EXT obs per profile: min={valid_obs2.min()}, max={valid_obs2.max()}, mean={valid_obs2.mean():.2f}")
    print(f"  N Profiles used in final value : {used_in_integral2}")


    ### Plot
    # f_mpl_plot2d(varin[0:100,indexes], xtime[indexes], z_in[0:100]/1000, strings_plot,
    #           colorange=(1e-6, 1e-2), scale='log', window_ave=0, esY=1)
