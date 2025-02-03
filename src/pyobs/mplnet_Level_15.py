#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script reads MPLnet data into memory:

A NetCDFData class to store the data and attributes

A read_netcdf function that:
    Takes a filename (required) and
    (optional) list of variables to read (default variables if none specified)
    Reads data & attributes for each variable
    Stores everything in the class structure
    Verb=0(default), no print,=1 variables names and sizes, =2 var names and attributes

Main Subroutines:
    read_file : actual read single file
    get_path  : set path where files are stored (needed if working in different clusters)


Created on Fri Jan 17 10:51:10 2025

@author: sgasso + chat.gsfc
"""
from netCDF4 import Dataset
# import numpy as np
import platform
import sys

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
    if now_os=='win32':
        base_path    = 'C:/Users/sgasso/Downloads/'
        pth_fig_out='C:/Users/sgasso/OneDrive - NASA/ToShare/2025/GEOS/Pyfigures/'
    elif now_os == 'darwin':
        base_path='/Users/sgasso/Downloads/'
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

######----------------------------------------------


if __name__ == "__main__":
    # Option 1: Direct call to the class method
    current_os    = platform.system()
    computer_name = platform.node()
    current_pth,pth_fig = get_path(current_os.lower(),computer_name)
    filename = "MPLNET_V3_L15_AER_20190603_MPL44258_GSFC.nc4"
    full_pathname = current_pth + filename

    # Define variables to read
    variables_to_read = ['time', 'altitude', 'extinction']

    # Read the data
    data = MPL_L15.read_file(full_pathname, variables=variables_to_read, verb=2)

    ### now access data directly
    print(data.time.shape)
    print(data.altitude.shape)
    print(data.extinction.shape)

    # Access attributes
    print(data.time_attributes)
    print(data.extinction_attributes.get('units'))