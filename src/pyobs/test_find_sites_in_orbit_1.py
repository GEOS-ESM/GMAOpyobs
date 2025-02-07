#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Jan/28/2025 - This code takes the output of test_read_level2_geo_2.py thon
(the structure called structured_data) , then it loops in every orbit and find

@author: sgasso
"""
import os,sys,platform
import glob
import h5py as h5
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta

from tropomaer import TROPOAER_L2

def find_hdf5_files(folder_path):
    """Find all HDF5 files in the specified folder"""
    hdf5_files = []
    for file in Path(folder_path).glob('*.hdf5'):  # Also add '*.h5' if needed
        hdf5_files.append(str(file))
    return hdf5_files

def process_hdf5_file(filename, use_indices=False, i_range=None, j_range=None, use_spacing=False, spacing=(1,1)):
    """Process individual HDF5 file with index-based subsetting"""
    try:
        with h5.File(filename, "r") as idf:
            print(filename)
            ID_geo = idf['GEODATA']
            # FROM HDF5 file: Get Reference time of the measurements.
            # The reference time is set to yyyy-mm-ddT00:00:00 UTC,
            # where yyyy-mm-dd is the day on which the measurements of a
            # particular data granule start.
            time=ID_geo['time'][()]
            # Now get first and last element in Time array, milliseconds since start of day
            dt  =np.array([ID_geo['delta_time'][0],ID_geo['delta_time'][-1]])
            if use_indices:
                i_start, i_end = i_range
                j_start, j_end = j_range

                if use_spacing:
                    # Apply both index ranges and different spacings for each dimension
                    x_spacing, y_spacing = spacing
                    lat = ID_geo['latitude'][i_start:i_end:x_spacing, j_start:j_end:y_spacing]
                    lon = ID_geo['longitude'][i_start:i_end:x_spacing, j_start:j_end:y_spacing]
                else:
                    # Apply only index ranges
                    lat = ID_geo['latitude'][i_start:i_end, j_start:j_end]
                    lon = ID_geo['longitude'][i_start:i_end, j_start:j_end]
            else:
                if use_spacing:
                    # Apply only spacing to full array with different spacings
                    x_spacing, y_spacing = spacing
                    lat = ID_geo['latitude'][::x_spacing, ::y_spacing]
                    lon = ID_geo['longitude'][::x_spacing, ::y_spacing]
                else:
                    # Read full arrays
                    lat = ID_geo['latitude'][()]
                    lon = ID_geo['longitude'][()]

            return lat, lon,time,dt
    except Exception as e:
        print(f"Error processing file {filename}: {str(e)}")
        return None, None,None,None


def get_array_dimensions(filename):
    ####Get dimensions from HDF5 file"""
    try:
        with h5.File(filename, "r") as idf:
            n_scanlines = len(idf['scanline'])
            n_pixels    = len(idf['ground_pixel'])
            #Note that we use .decode() because the attribute is stored as a
            # bytes object.This converts it to a regular Python string
            time_start_str=idf.attrs['time_coverage_start'].decode()
            time_end_str  =idf.attrs['time_coverage_end'].decode()
            # Convert to datetime object
            time_start_dt = datetime.strptime(time_start_str, "%Y-%m-%dT%H:%M:%SZ")
            time_end_dt   = datetime.strptime(time_end_str, "%Y-%m-%dT%H:%M:%SZ")
            return n_scanlines, n_pixels,time_start_dt,time_end_dt
    except Exception as e:
        print(f"Error reading dimensions from file {filename}: {str(e)}")
        return None, None,None,None



def get_path(now_os, now_computer):
    # Your existing get_path function remains the same
    if now_os=='win32':
        base_path    = 'C:/Users/sgasso/Downloads/'
        pth_fig_out='C:/Users/sgasso/OneDrive - NASA/ToShare/2025/GEOS/Pyfigures/'
    elif now_os == 'darwin':
        base_path= '/Volumes/ExtData1/SatData/Tropomi/Level2/2023/359/'
        pth_fig_out='/Users/sgasso/Library/CloudStorage/OneDrive-NASA/ToShare/2025/AOCH/PyFigures/'
    elif now_os == 'linux' and "calculon" in now_computer:
        base_path = '/nobackup/TROPOMAER/'
        pth_fig_out = ''
    elif now_os == 'linux' and "discover" in now_computer:
        base_path = '/nobackup/CALIPSO/Level1.5/Santiago/'
    else:
        print('Current operating system no recognized.')
        print('Cannot set path to MPL  files. Terminate Run')
        sys.exit()
    return base_path,pth_fig_out
    
if __name__ == "__main__":
    current_working_directory = os.getcwd()
    print('\nThis code is running from directory :', current_working_directory)
    # Get user inputs
    current_os    = platform.system()
    computer_name = platform.node()
    print(f'This code is executed in computer {computer_name} \n')
    current_pth,pth_fig = get_path(current_os.lower(),computer_name)
    # filename='TROPOMI-Sentinel-5P_L2-TROPOMAER_2023m1224t234903-o32115_v01-2023m1231t055100.nc'
    #filename='TROPOMI-Sentinel-5P_L2-TROPOMAER_2023m1225t013033-o32116_v01-2023m1231t050855.nc'
    #folder_path = '/Volumes/ExtData1/SatData/Tropomi/Level2/2023/359/'
    
    #### This list was created by running test_read_level2_geo_2.py
    #### Using a local dictionary now but eventually this will be the output of test_read_level2_geo_2.py
    matching_orbits_list={\
    'TROPOMI-Sentinel-5P_L2-TROPOMAER_2023m1225t132102-o32123_v01-2023m1231t105113.nc':\
     [['Santa_Cruz_Tenerife',28.472,-16.247,0.052,datetime(2023, 12, 25, 13, 15),datetime(2023, 12, 25, 15,15)]],\
    'TROPOMI-Sentinel-5P_L2-TROPOMAER_2023m1225t164402-o32125_v01-2023m1230t212055.nc': \
     [['GSFC',38.993,-76.84,0.05,datetime(2023, 12, 25, 17, 5),datetime(2023, 12, 25, 19, 5)]],\
    'TROPOMI-Sentinel-5P_L2-TROPOMAER_2023m1225t182532-o32126_v01-2023m1230t205012.nc': \
    [['GSFC',38.993,-76.84,0.05,datetime(2023, 12, 25, 17, 5),datetime(2023, 12, 25, 19, 5)]]
    }

    # Process each orbit and its associated sites
    for orbit_file, sites_list in matching_orbits_list.items():
        print(f"With orbit: {orbit_file}")
        ### 
        Files = current_pth+'2023/359/'+orbit_file
        # Extract only GEO data
        # TROP = TROPOAER_L2(Files,GEO=True,Verbose=1)
        TROP = TROPOAER_L2(Files,GEO=True,Verbose=1)
        
        sys.exit()
        # Process each site within this orbit
        for site_data in sites_list:
            site_name = site_data[0]    # Site name
            lat = site_data[1]          # Latitude
            lon = site_data[2]          # Distance
            dist = site_data[3]         # Distance
            start_time = site_data[4]   # Start datetime
            end_time = site_data[5]     # End datetime
            
            print(f"  Processing site: {site_name}")
            #OK 
            #1) Load lat/long
            #2) compute distance or use where (check time difference)
            
            
            # Your processing operations here
            # Example:
            # process_site_data(orbit_file, site_name, lat, lon, start_time, end_time, ...)
    
