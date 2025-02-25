#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Feb12 - This is working

Feb/6/2025 - based on version _3, this code incorporates some improvements towards final code
such as scanning for whether the site contains target site, operationalize more respect to previous version.

This code works
Feb/1/2025 _ same as _v2 but now it upload the actual TROPOMI file and carries out
the distance calculation
Jan/28/2025 - This code takes the output of test_read_level2_geo_2.py thon
(the structure called structured_data) , then it loops in every orbit and find
line and row indexes in the orbit and saves them

@author: sgasso
"""
import os,sys,platform
import glob
import h5py as h5
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta

from tropomi_l2_reader import TROPOAER_L2

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



def get_path(now_os, now_computer):
    #####  Your existing get_path function remains the same
    if now_os=='win32':
        base_path    = 'C:/Users/sgasso/Downloads/'
        pth_fig_out='C:/Users/sgasso/OneDrive - NASA/ToShare/2025/GEOS/Pyfigures/'
    elif now_os == 'darwin':
        base_path= '/Volumes/ExtData1/SatData/Tropomi/Level2/2023/359/'
        pth_fig_out='/Users/sgasso/Library/CloudStorage/OneDrive-NASA/ToShare/2025/AOCH/PyFigures/'
    elif now_os == 'linux' and "calculon" in now_computer:
        base_path = '/nobackup/TROPOMAER/2023/359/'
        pth_fig_out = ''
    elif now_os == 'linux' and "discover" in now_computer:
        base_path = '/nobackup/CALIPSO/Level1.5/Santiago/'
    else:
        print('Current operating system no recognized.')
        print('Cannot set path to MPL  files. Terminate Run')
        sys.exit()
    return base_path,pth_fig_out
    
# def find_nearest_point(lat_array, lon_array, lat0, lon0):
    # # Compute absolute differences
    # lat_diff = np.abs(lat_array - lat0)
    # lon_diff = np.abs(lon_array - lon0)
    
    # # Combined difference
    # diff = lat_diff + lon_diff
    
    # # Get index of minimum value
    # idx = np.argmin(diff)
    
    # # Convert to 2D indices
    # row, col = np.unravel_index(idx, lat_array.shape)
    
    # return row, col

#### Read list of MPL sites
def read_mpl_site_list(filename):
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
#             datetime1 = datetime.combine(day0, time1)
#             datetime2 = datetime.combine(day0, time2)
            site = [
                data[0],
                float(data[1]),
                float(data[2]),
                float(data[3]),
                time1,
                time2
            ]
            sites.append(site)
    return sites


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
    # print(dist0)
    return row, col,dist0
    
def create_output_files(output_array, year, userstring):
    # Define header
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
def process_matching_orbits(matching_orbits_list, current_working_directory, current_pth, SEL_DATE, REF_DATE,mpl_sites_list):
    ### now load list of MPL sites
#     pth_site_list = current_working_directory + '/' + 'list_mpl_sites.txt'
#     mpl_sites_list = read_mpl_site_list(pth_site_list, SEL_DATE)
    
    # Initialize output dictionary
    output_array = {}
    for site in mpl_sites_list:
        output_array[site[0]] = []
        
    #### now loop over each element in matching_orbits_list and check if the orbit needs to be 
    #### loaded or is already in memory.
    current_orbit = None
    current_file = None
    start_time = datetime.now()
    
    for entry in matching_orbits_list:
        orbit_num    = entry[0]
        filename_tro = entry[1]
        # Check if we need to load a new file or use existing orbit
        if orbit_num != current_orbit:
            print(f"\nLoading new file for orbit {orbit_num}: {filename_tro}")
            file_pth = current_pth+filename_tro
            # Extract only GEO data
            TROP = TROPOAER_L2(file_pth,GEO=True,Verbose=0)
            # saving for next loop
            current_orbit = orbit_num 
        else:
            print(f"\nUsing existing file for orbit {orbit_num}")

        # Now get site data
        site_name = entry[2]
        lat0      = entry[3]
        lon0      = entry[4]
        obs_time1  = entry[6]
        obs_time2  = entry[7]

        print(f"   Processing site: {site_name}")
        print(f"   Operating on data for {site_name} at ({lat0:.4f}, {lon0:.4f}), OverpassTimeFrame {obs_time1.strftime('%H:%M')}  ,  {obs_time2.strftime('%H:%M')}")
        line, col,distance = find_nearest_point_haversine(TROP.lat, TROP.lon, lat0, lon0)
        
        if distance > 10 :
            print(f"          Closest TROPO pixel is {distance:.2f} km away from site. Collocation skipped")
        else:
            found_lat = TROP.lat[line, col]
            found_lon = TROP.lon[line, col]
            refday_sec= int(TROP.time[0])
            found_time= float(TROP.dtime[line])
            date1    = (REF_DATE + timedelta(seconds=int(TROP.time[0])))
            date1_str= date1.strftime('%Y-%m-%d %H:%M:%S.f')[:-4]
            date2     = date1 + timedelta(milliseconds=int(TROP.dtime[0]))
            date2_str = date2.strftime('%Y-%m-%d %H:%M:%S.f')[:-4]

            ### Print 
            print(f"    TROPOMI line,column {line},{col} , error in distance collocation {distance:.2f}  ")
            print(f"         Time scanline at site: {date2_str}")
            print(f"         For   lat,lon  {found_lat:.4f},{found_lon:.4f} ")
            print(f"    MPL site located at {lat0:.4f},{lon0:.4f}")
            
            if not output_array[site_name]:  # checks if the list is empty
                output_array[site_name] = [(date2_str, orbit_num, filename_tro, line, col)]
            else:
                output_array[site_name].append((date2_str, orbit_num, filename_tro, line, col))
                
    end_time = datetime.now()
    time_diff = (end_time - start_time).total_seconds()
    print(f"\nExecution time: {time_diff:.2f} seconds")
    
    return output_array


if __name__ == "__main__":
    REF_DATE=datetime(2010, 1, 1) # Reference date for Time Array in Tropomi SDS variable
    current_working_directory = os.getcwd()
    print('\nThis code is running from directory :', current_working_directory)
    # Get user inputs
    current_os    = platform.system()
    computer_name = platform.node()
    print(f'This code is executed in computer {computer_name} \n')
    current_pth,pth_fig = get_path(current_os.lower(),computer_name)
    
    #### Get user defined text file with MPL sites to process
    pth_site_list  = current_working_directory + '/' + 'list_mpl_sites.txt'
    mpl_sites_list = read_mpl_site_list(pth_site_list)
    #### The next list was created by running test_read_level2_geo_2.py
    #### Using a local dictionary now but eventually this will be the output of test_read_level2_geo_2.py
    SEL_DATE = datetime(2023,12,25) # 
    ### Now add SEL_DATE to mpl time stamps, this is needed to get the matching orbits 
    for i in mpl_sites_list:
         i[-1]= datetime.combine(SEL_DATE, i[-1])
         i[-2]= datetime.combine(SEL_DATE, i[-2])
     
    #sys.exit()
    #### now need to get a list of matching orbits, that is a list of orbits that contain a 
    ####  MPL site.
    matching_orbits_list=[\
    [32123,'TROPOMI-Sentinel-5P_L2-TROPOMAER_2023m1225t132102-o32123_v01-2023m1231t105113.nc',
  'Santa_Cruz_Tenerife',28.472,-16.247,0.052,datetime(2023, 12, 25, 13, 15),datetime(2023, 12, 25, 15, 15)],
    [32125,'TROPOMI-Sentinel-5P_L2-TROPOMAER_2023m1225t164402-o32125_v01-2023m1230t212055.nc',
  'GSFC',38.993,-76.84,0.05,datetime(2023, 12, 25, 17, 5),datetime(2023, 12, 25, 19, 5)],
    [32126,'TROPOMI-Sentinel-5P_L2-TROPOMAER_2023m1225t182532-o32126_v01-2023m1230t205012.nc',
  'GSFC',38.993,-76.84,0.05,datetime(2023, 12, 25, 17, 5),datetime(2023, 12, 25, 19, 5)]\
                         ]
                         
    ##### Now find the locations of the MPL sites inside the selected TROPomi orbits
    output_array = process_matching_orbits(matching_orbits_list, current_working_directory, 
                                         current_pth, SEL_DATE, REF_DATE,mpl_sites_list)
    
    ##### Next save outfile one per site
    ### Inputs: array_to_save, year_string,user_string
    create_output_files(output_array, "2023", "test")

# To do, finish editing read_mpl_site_list module, remove input SEL_DATE ,read data as time only and add 
# SEL_DATE later. 