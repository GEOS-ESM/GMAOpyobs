#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Jan/28/2025 - it works. it outputs a structure with following format
[tropo_orbit_filename, MPL_site_name, site_lat,site_lon,site_alt,start_time,end_time]
that can be used late for finding the matches.  

This code test reading lat,lon and time tags from a L2 Tropo file and compares
with time of tytpical from the site of interest.


Created on Fri Jan 24 11:56:58 2025

@author: sgasso
"""
import os,sys,platform
import glob
import h5py as h5
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta

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

#### Read list of MPL sites
def read_mpl_site_list(filename,day0):
#Expected format: comma separated
#line 1 :Name_MPL_sites , name,lat, lon, altiude(km),time_start, time_end
#line 2 : GSFC           ,38.9930,-76.8400,0.050,    hh:mm , hh:mm
#.............
#
    sites = []
    with open(filename, 'r') as file:
        header = file.readline()
        for line in file:
            data = [item.strip() for item in line.split(',')]
           # Convert time strings to time objects
            time1 = datetime.strptime(data[4], "%H:%M").time()
            time2 = datetime.strptime(data[5], "%H:%M").time()
           # Combine date with times
            datetime1 = datetime.combine(day0, time1)
            datetime2 = datetime.combine(day0, time2)
            site = [
                data[0],
                float(data[1]),
                float(data[2]),
                float(data[3]),
                datetime1,
                datetime2
            ]
            sites.append(site)
    return sites


def find_matching_orbits(orbits_list, sites_list):
    """
    Find orbits that have matching sites based on time overlaps.
    
    Parameters:
    orbits_list: list of tuples (filename, start_time, end_time)
    sites_list: list of tuples (name, lat, lon, altitude, start_time, end_time)
    
    Returns:
    Dictionary with orbit filenames as keys and lists of matching site names as values
    """
    orbit_matches = {}

    # For each orbit
    for orbit in orbits_list:
        filename, orbit_start, orbit_end = orbit
        matching_sites = []

        # Check each site for this orbit
        for site in sites_list:
            site_name, _, _, _, site_start, site_end = site
            
            # Check if site observation time overlaps with orbit time
            if (orbit_start <= site_end and orbit_end >= site_start):
                matching_sites.append(site_name)
        
        # Only include orbits that have matching sites
        if matching_sites:
            orbit_matches[filename] = matching_sites

    return orbit_matches

######----------------------------------------------


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

    file_list = glob.glob(current_pth+'/2023/359/*.nc')

    orbits_list=[]
    for full_pathname in file_list:
       
      ### Get some global info from file
      Nlines, Ncols,t_start,t_end = get_array_dimensions(full_pathname)
      # print(t_start.strftime("%Y-%m-%d %H:%M:%S"))
      # print(t_end.strftime("%Y-%m-%d %H:%M:%S"))
      # Store in array if needed
      orbit_time = [t_start, t_end]
#      # print('Orbit: ',full_pathname[-80:])
#      # print('     Start,End time',t_start.strftime("%H:%M:%S"),t_end.strftime("%H:%M:%S"))
      orbits_list.append((full_pathname[-80:],orbit_time[0],orbit_time[1]))

    ### Save date for converting time in MPL sites
    sel_date=orbits_list[1][1].date() # it returns only the day

    ### now load list of MPL sites
    pth_site_list=current_working_directory + '/' + 'list_mpl_sites.txt'
    mpl_sites_list = read_mpl_site_list(pth_site_list,sel_date)
    # print("\nPython List of Dictionaries:")
    #for site in mpl_sites_list:
    #    print(site)
    
    #### Print Sites header
    print("\n------------------------------")
    print(f"{'Location':<20} {'Latitude':>8} {'Longitude':>9} {'Altitude':>7} {'Start Time':>15}   {'End Time':>15}")
    # print("-" * 85)  # Separator line
    for entry in mpl_sites_list:
       print(f"{entry[0]:<20} {entry[1]:>8.3f} {entry[2]:>9.3f} {entry[3]:>7.3f} {entry[4].strftime('%Y-%m-%d %H:%M'):>20} {entry[5].strftime('%Y-%m-%d %H:%M'):>20}")
       
    ### Compare the orbit time frame with site time frame
    matching_orbits = find_matching_orbits(orbits_list, mpl_sites_list)

    # Print results
    print("------------------------------\n")
    print("\nMatching sites for each orbit:")
    for orbit, sites in matching_orbits.items():
        print(f"Orbit: {orbit}")
        print("  Matching sites:")
        for site in sites:
            print(f"  - {site}")
        print("")
        
    # # Create structured list
    # structured_data = []
    # for orbit, sites in matching_orbits.items():
        # for site in sites:
            # structured_data.append([orbit, site])
    # # Print the structured data
    # print("\nStructured Data:")
    # print("-" * 85)
    # print(f"{'Orbit':<70} {'Site':<20}")
    # print("-" * 85)
    # for entry in structured_data:
        # print(f"{entry[0]:<70} {entry[1]:<20}")
    
    # #### Create structured list
    # structured_data = []
    # for orbit, site_names in matching_orbits.items():
        # for site_name in site_names:
            # # Find the complete site information from mpl_sites_list
            # site_info = next(site for site in mpl_sites_list if site[0] == site_name)
            # structured_data.append([
                # orbit,             # orbit filename
                # site_info[0],      # site name
                # site_info[1],      # latitude
                # site_info[2],      # longitude
                # site_info[3],      # altitude
                # site_info[4],      # start time
                # site_info[5]       # end time
            # ])
    # ##### Print the structured data
    # print("\nStructured Data:")
    # print("-" * 140)
    # print(f"{'Orbit':<70} {'Site':<20} {'Lat':>8} {'Lon':>9} {'Alt':>7} {'Start Time':>15} {'End Time':>15}")
    # print("-" * 140)
    # for entry in structured_data:
        # print(f"{entry[0]:<70} {entry[1]:<20} {entry[2]:>8.3f} {entry[3]:>9.3f} {entry[4]:>7.3f} {entry[5].strftime('%H:%M'):>15} {entry[6].strftime('%H:%M'):>15}")
        
    # Create dictionary of orbits and their matching sites
    orbit_sites = []
    for orbit, site_names in matching_orbits.items():
        orbit_num=int(orbit[51:56])
        for site_name in site_names:
            # Find the complete site information from mpl_sites_list
            site_info = next(site for site in mpl_sites_list if site[0] == site_name)
            location_info=[orbit_num,orbit]+site_info
            orbit_sites.append(location_info)

            # If orbit already exists as key, append the site info
            # If not, create new list with site info
#             if orbit_num in orbit_sites:
#                 orbit_sites[orbit_num].append(location_info)
#             else:
#                 orbit_sites[orbit_num] = [location_info]

#             if orbit_num in orbit_sites:
#            orbit_sites.append(location_info)
#             else:
#                 orbit_sites = [location_info]
#     sys.exit()
#     # Print the dictionary contents
#     print("\nOrbit-Sites Dictionary:")
#     print("-" * 80)
#     for orbit, sites in orbit_sites.items():
#         print(f"\nOrbit: {orbit}")
#         print(f"{'Site':<20} {'Lat':>8} {'Lon':>9} {'Alt':>7} {'Start Time':>15} {'End Time':>15}")
#         print("-" * 80)
#         for site in sites:
#             print(f"{site[0]:<20} {site[1]:>8.3f} {site[2]:>9.3f} {site[3]:>7.3f} {site[4].strftime('%H:%M'):>15} {site[5].strftime('%H:%M'):>15}")

    print("\nFormatted Orbit-Site Information:")
    print("-" * 100)  # Print a separator line
    print(f"{'Orbit':^8} {'Site':^20} {'Latitude':^10} {'Longitude':^10} {'Altitude':^10} {'Start Time':^20} {'End Time':^20}")
    print("-" * 100)  # Print a separator line
    
    for entry in orbit_sites:
        print(f"{entry[0]:^8} {entry[2]:^20} {entry[3]:^10.3f} {entry[4]:^10.3f} {entry[5]:^10.3f} {entry[6].strftime('%Y-%m-%d %H:%M'):^20} {entry[7].strftime('%Y-%m-%d %H:%M'):^20}")
    
    print("-" * 100)  # Print a separator line