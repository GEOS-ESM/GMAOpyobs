#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Apr/8/2025 Added an output path for output files

Mar/3/2025 . Set of subroutines called by collocate_tropo_mpl_2.py . Code requires the module pyobs.tropomi_l2_reader

@author: sgasso
"""
import os,sys,platform
import glob
import h5py as h5
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta

from src.pyobs.tropoaer import TROPOAER_L2 # need to read TROPOMI data
#from tropoaer import TROPOAER_L2 # need to read TROPOMI data


#------------ Functions related to reach the files 

def find_hdf5_files(folder_path):
    ##### Find all HDF5 files in the specified folder
    hdf5_files = []
    for file in Path(folder_path).glob('*.hdf5'):  # Also add '*.h5' if needed
        hdf5_files.append(str(file))
    return hdf5_files

def get_path(now_os, now_computer):
    #### Module to setup path to where the datafiles are. Needed if working in multiple computers
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
def read_mpl_site_list(filename):
#Expected format: comma separated
#line 1 :Name_MPL_sites , name,lat, lon, altiude(km),time_start, time_end
#line 2 : GSFC           ,38.9930,-76.8400,0.050,    hh:mm , hh:mm
#.............
# Last two columns are ignored 
    sites = []
    with open(filename, 'r') as file:
        header = file.readline()
        for line in file:
            data = [item.strip() for item in line.split(',')]
            site = [
                data[0],
                float(data[1]),
                float(data[2]),
                float(data[3]),
            ]
            sites.append(site)
    return sites


# def get_orbit_site_matches(list_mpl_sites,yyyy, mm=None, dd=None, julian=None, verb=0):
# #--------------------------------------------------------------------------------------------
# #     Inputs:
# #     list_mpl_sites : 
# #     yyyy : int , Year
# #     mm : int, optional ,  Month (1-12)
# #     dd : int, optional ,  Day (1-31)
# #     julian : int, optional , Julian day (alternative to mm/dd)
# #     verb : int, optional ,  Verbose output flag (default: 0,0: minimal print ,1: some print 2: print all)
# #         
# #     Returns: List containing files with files in selected directory and the path to the folder
# #--------------------------------------------------------------------------------------------
    
    # if not yyyy:
        # raise ValueError('Input year must be provided')
    # if not mm and not julian:
        # raise ValueError('Either Month and day or Julian day must be provided')

    # # Set correct paths according to the current computer and OS
    # current_working_directory = os.getcwd()
    # current_os           = platform.system()
    # computer_name        = platform.node()
    # current_pth, pth_fig = get_path(current_os.lower(), computer_name)
    
    # # Convert mm/dd to julian if needed
    # if mm:
        # date_obj = datetime(yyyy, mm, dd)
        # julian = (date_obj - datetime(yyyy, 1, 1)).days + 1
    # julian_str = f"{julian:03d}"

    # if verb>-1:print(f'   Year is {yyyy}, julian {julian_str}, {date_obj:%Y-%m-%d}')

    # # Get list of files for selected date
    # path_2_folder_with_orbits = current_pth+str(yyyy)+'/'+julian_str +'/'
    # file_list = glob.glob(path_2_folder_with_orbits+'*.nc')
    
    # # Process orbit information
    # orbits_list = []
    # for full_pathname in file_list:
        # orbits_list.append(full_pathname[-80:])

    # return orbits_list,path_2_folder_with_orbits

def get_orbit_site_matches(list_mpl_sites, yyyy, mm=None, dd=None, julian=None, verb=0):
    #--------------------------------------------------------------------------------------------
    #     Inputs:
    #     list_mpl_sites : 
    #     yyyy   : int , Year
    #     mm     : int, optional ,  Month (1-12)
    #     dd     : int, optional ,  Day (1-31)
    #     julian : int, optional , Julian day (alternative to mm/dd)
    #     verb   : int, optional ,  Verbose output flag (default: 0,0: minimal print ,1: some print 2: print all)
    #         
    #     Returns: List containing files with files in selected directory and the path to the folder
    #--------------------------------------------------------------------------------------------
    
    if not yyyy:
        raise ValueError('Input year must be provided')
    if not mm and not julian:
        raise ValueError('Either Month and day or Julian day must be provided')

    # Set correct paths according to the current computer and OS
    current_working_directory = os.getcwd()
    current_os           = platform.system()
    computer_name        = platform.node()
    current_pth, pth_fig = get_path(current_os.lower(), computer_name)
    
    # Convert mm/dd to julian if needed
    if mm:
        date_obj = datetime(yyyy, mm, dd)
        julian = (date_obj - datetime(yyyy, 1, 1)).days + 1
    julian_str = f"{julian:03d}"

    if verb>-1:print(f'   Year is {yyyy}, julian {julian_str}, {date_obj:%Y-%m-%d}')

    # Get list of files for selected date
    path_2_folder_with_orbits = current_pth+str(yyyy)+'/'+julian_str +'/'
    
    # Check if folder exists
    if not os.path.isdir(path_2_folder_with_orbits):
        if verb>0:
            print(f"Folder for julian day {julian_str} does not exist. Returning empty outputs.")
        return [], path_2_folder_with_orbits  # Return empty list and the path (even though it doesn't exist)
    
    # Get list of files for selected date
    file_list = glob.glob(path_2_folder_with_orbits+'*.nc')
    
    # Filter files by size (non-zero size)
    valid_files = []
    for file_path in file_list:
        try:
            file_size = os.path.getsize(file_path)
            if file_size > 0:  # Keep only non-zero size files
                valid_files.append(file_path)
            elif verb>1:
                print(f"Skipping zero-size file: {os.path.basename(file_path)}")
        except OSError:
            if verb>1:
                print(f"Error accessing file: {os.path.basename(file_path)}")
    
    # Update file_list with only valid files
    file_list = valid_files
    
    # Process orbit information
    orbits_list = []
    for full_pathname in file_list:
        orbits_list.append(full_pathname[-80:])

    return orbits_list, path_2_folder_with_orbits
    
#------------------------- 
#---------- functions related to finding sites in each orbit
#-------------------------

def find_hdf5_files(folder_path):
   #### Find all HDF5 files in the specified folder
    hdf5_files = []
    for file in Path(folder_path).glob('*.hdf5'):  # Also add '*.h5' if needed
        hdf5_files.append(str(file))
    return hdf5_files

def find_nearest_point_haversine(lat_array, lon_array, lat0, lon0):
    # Earth's radius in kilometers
    R = 6371.0
    
    # First find approximate nearest point using simpler calculation
    lat_diff = np.abs(lat_array - lat0)
    lon_diff = np.abs(lon_array - lon0)
    diff = lat_diff + lon_diff
    
    # Get index of minimum value
    idx = np.argmin(diff)
    row, col = np.unravel_index(idx, lat_array.shape)
    
    # Now compute accurate Haversine distance only for this point
    lat1 = np.radians(lat0)
    lon1 = np.radians(lon0)
    lat2 = np.radians(lat_array[row, col])
    lon2 = np.radians(lon_array[row, col])
    
    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    dist0 = R * c  # Distance in kilometers
    return row, col,dist0
    
def create_output_files(output_path,output_array, year, userstring):
    header = 'yyyy-mm-dd hh:mm:ss.d,orbit_number,source_filename     ,line,column'
    # get the key with the site name
    for key in output_array.keys():
        # Create filename using the specified format
        filename = f"{key}_{year}_{userstring}.txt"
        print('...... Saving to file ', filename)
        # Open file for writing
        with open(output_path+filename, 'w') as f:
          # Write header first
            f.write(header + '\n')
            # Loop through each tuple in the array
            for row in output_array[key]:
                # Convert all elements to strings and join with commas
                line = ','.join(str(item) for item in row)
                # Write line to file with newline
                f.write(line + '\n')
                
def process_matching_orbits(orbits_list, \
                           current_pth, SEL_DATE, mpl_sites_list, verb=True):

    REF_DATE=datetime(2010, 1, 1) # Tropomi reference date

    # Initialize output dictionary
    output_array = {}
    for site in mpl_sites_list:
        output_array[site[0]] = []
        
    #### now loop over each element in orbits_list and check if the orbit needs to be 
    #### loaded or is already in memory.
    current_orbit = None
    current_file = None
    start_time2 = datetime.now() # start timing the loop
    
    for entry in orbits_list:

        filename_tro = entry
        orbit_num    = filename_tro[51:56]

        # Check if we need to load a new file or use existing orbit
        if orbit_num != current_orbit:
            if verb>1:print(f"   Loading new file for orbit {orbit_num}: {filename_tro}")
            file_pth = current_pth+filename_tro
            # Extract only GEO data
            TROP = TROPOAER_L2(file_pth,GEO=True,Verbose=0)
            # saving for next loop
            current_orbit = orbit_num 
        else:
            if verb>1:print(f"   Using existing file for orbit {orbit_num}")
        # For current orbit, loop over every site and find distance
        for site in mpl_sites_list:
            # Now get site data
            site_name = site[0]
            lat0      = site[1]
            lon0      = site[2]
            if verb: print(f"   Processing site: {site_name}")
            line, col,distance = find_nearest_point_haversine(TROP.lat, TROP.lon, lat0, lon0)
            # check distance and determine if site is there. 
            if distance > 10 :
                if verb>1:print(f"          Closest TROPO pixel is {distance:.2f} km away from site. \n     <------ Collocation skipped ---->")
            else:
                found_lat = TROP.lat[line, col]
                found_lon = TROP.lon[line, col]
                refday_sec= int(TROP.time[0])
                found_time= float(TROP.dtime[line])
                date1     = (REF_DATE + timedelta(seconds=int(TROP.time[0]))) # reference time from global time.
                date1_str = date1.strftime('%Y-%m-%d %H:%M:%S.f')[:-4]
                date2     = date1 + timedelta(milliseconds=int(TROP.dtime[0])) # time since begining of orbit
                date2_str = date2.strftime('%Y-%m-%d %H:%M:%S')
                if verb>0:
                    ### Print 
                    print(f"    TROPOMI line,column {line},{col} , error in distance collocation {distance:.2f}  ")
                    print(f"         Time scanline at site: {date2}")
                    print(f"         For   lat,lon  {found_lat:.4f},{found_lon:.4f} ")
                    print(f"    MPL site located at {lat0:.4f},{lon0:.4f}")
                
                if not output_array[site_name]:  # checks if the list is empty
                    output_array[site_name] = [(date2, orbit_num, filename_tro, line, col)]
                else:
                    output_array[site_name].append((date2, orbit_num, filename_tro, line, col))
                    
    end_time = datetime.now()
    time_diff = (end_time - start_time2).total_seconds()
    if verb>0: print(f"    Time this orbit: {time_diff:.2f} seconds")
    
    return output_array