#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Feb/12/2025 . Set of subroutines called by test_collocate_tropo_mpl_1.py


@author: sgasso
"""
import os,sys,platform
import glob
import h5py as h5
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta

from tropomi_l2_reader import TROPOAER_L2


def get_path(now_os, now_computer):
    #####  Your existing get_path function remains the same
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

def find_hdf5_files(folder_path):
   #### Find all HDF5 files in the specified folder
    hdf5_files = []
    for file in Path(folder_path).glob('*.hdf5'):  # Also add '*.h5' if needed
        hdf5_files.append(str(file))
    return hdf5_files

def process_hdf5_file(filename, use_indices=False, i_range=None, j_range=None, use_spacing=False, spacing=(1,1)):
    #### Process individual HDF5 file with index-based subsetting"""
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
    
def create_output_files(output_array, year, userstring):
    header = 'yyyy-mm-dd hh:mm:ss.d,orbit_number,source_filename     ,line,column'
    for key in output_array.keys():
        # Create filename using the specified format
        filename = f"{key}_{year}_{userstring}.txt"
        print('...... Saving to file ', filename)
        # Open file for writing
        with open(filename, 'w') as f:
          # Write header first
            f.write(header + '\n')
            # Loop through each tuple in the array
            for row in output_array[key]:
                # Convert all elements to strings and join with commas
                line = ','.join(str(item) for item in row)
                # Write line to file with newline
                f.write(line + '\n')
def process_matching_orbits(matching_orbits_list, \
                         current_pth, SEL_DATE, mpl_sites_list, verb=True):
    REF_DATE=datetime(2010, 1, 1) # Tropomi reference date
    # Initialize output dictionary
    output_array = {}
    for site in mpl_sites_list:
        output_array[site[0]] = []
        
    #### now loop over each element in matching_orbits_list and check if the orbit needs to be 
    #### loaded or is already in memory.
    current_orbit = None
    current_file = None
    start_time = datetime.now() # start timing the loop
    
    for entry in matching_orbits_list:
        orbit_num    = entry[0]
        filename_tro = entry[1]
        # Check if we need to load a new file or use existing orbit
        if orbit_num != current_orbit:
            if verb:print(f"\nLoading new file for orbit {orbit_num}: {filename_tro}")
            file_pth = current_pth+filename_tro
            # Extract only GEO data
            TROP = TROPOAER_L2(file_pth,GEO=True,Verbose=0)
            # saving for next loop
            current_orbit = orbit_num 
        else:
            if verb:print(f"\nUsing existing file for orbit {orbit_num}")

        # Now get site data
        site_name = entry[2]
        lat0      = entry[3]
        lon0      = entry[4]
        obs_time1  = entry[6]
        obs_time2  = entry[7]
        if verb: 
            print(f"   Processing site: {site_name}")
            print(f"   Operating on data for {site_name} at ({lat0:.4f}, {lon0:.4f}), OverpassTimeFrame {obs_time1.strftime('%H:%M')}  ,  {obs_time2.strftime('%H:%M')}")
        line, col,distance = find_nearest_point_haversine(TROP.lat, TROP.lon, lat0, lon0)
        
        if distance > 10 :
            if verb:print(f"          Closest TROPO pixel is {distance:.2f} km away from site. Collocation skipped")
        else:
            found_lat = TROP.lat[line, col]
            found_lon = TROP.lon[line, col]
            refday_sec= int(TROP.time[0])
            found_time= float(TROP.dtime[line])
            date1    = (REF_DATE + timedelta(seconds=int(TROP.time[0])))
            date1_str= date1.strftime('%Y-%m-%d %H:%M:%S.f')[:-4]
            date2     = date1 + timedelta(milliseconds=int(TROP.dtime[0]))
            print(date2)
            date2_str = date2.strftime('%Y-%m-%d %H:%M:%S')
            if verb:
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
    time_diff = (end_time - start_time).total_seconds()
    if verb: print(f"\nExecution time: {time_diff:.2f} seconds")
    
    return output_array