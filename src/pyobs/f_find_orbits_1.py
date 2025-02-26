### created Feb/12/2025 : these are the modules used by test_create_list_orbits_with_sites_2.py


import os,sys,platform
import glob
import h5py as h5
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta


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


def process_hdf5_file(filename, use_indices=False, i_range=None, j_range=None, use_spacing=False, spacing=(1,1)):
    #### Process individual HDF5 file with index-based subsetting
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



#### Read list of MPL sites
def read_mpl_site_list(filename):
#Expected format: comma separated
#line 1 :Name_MPL_sites , name,lat, lon, altiude(km),time_start, time_end
#line 2 : GSFC           ,38.9930,-76.8400,0.050,    hh:mm , hh:mm
#.............

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


def find_matching_orbits(orbits_list, sites_list):
#     """
#     Find orbits that have matching sites based on time overlaps.
#     
#     Parameters:
#     orbits_list: list of tuples (filename, start_time, end_time)
#     sites_list: list of tuples (name, lat, lon, altitude, start_time, end_time)
#     
#     Returns:
#     Dictionary with orbit filenames as keys and lists of matching site names as values
#     """
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

def get_orbit_site_matches(list_mpl_sites,yyyy, mm=None, dd=None, julian=None, verb=True):
#     """
#     Module to match orbits with MPL sites based on temporal overlap.
#     
#     Parameters:
#     -----------
#     yyyy : int
#         Year
#     mm : int, optional
#         Month (1-12)
#     dd : int, optional
#         Day (1-31)
#     julian : int, optional
#         Julian day (alternative to mm/dd)
#     verb : bool, optional
#         Verbose output flag (default: True)
#         
#     Returns:
#     --------
#     list
#         List containing orbit-site matches with format:
#         [orbit_number, orbit_path, site_name, lat, lon, altitude, start_time, end_time]
#     """
    
    if not yyyy:
        raise ValueError('Input year must be provided')
    if not mm and not julian:
        raise ValueError('Either Month and day or Julian day must be provided')

    # Set correct paths
    current_working_directory = os.getcwd()
    current_os = platform.system()
    computer_name = platform.node()
    current_pth, pth_fig = get_path(current_os.lower(), computer_name)
    
    # Convert mm/dd to julian if needed
    if mm:
        date_obj = datetime(yyyy, mm, dd)
        julian = (date_obj - datetime(yyyy, 1, 1)).days + 1
    julian_str = f"{julian:03d}"
    
    # if verb:
        # print('     running test_f_find_orbits_1.py')
        # print(f'    Date is {yyyy},{mm},{dd}')
    if verb:print(f'    Year is {yyyy}, julian {julian_str}')
    # Get list of files for selected date
    path_2_folder_with_orbits = current_pth+str(yyyy)+'/'+julian_str +'/'
    file_list = glob.glob(path_2_folder_with_orbits+'*.nc')
    # Process orbit information
    orbits_list = []
    for full_pathname in file_list:
        Nlines, Ncols, t_start, t_end = get_array_dimensions(full_pathname)
        orbits_list.append((full_pathname[-80:], t_start, t_end))
    
    if verb:
        print("    ------------------------------")
        print(f"    {'Location':<20} {'Latitude':>8} {'Longitude':>9} {'Altitude':>7} {'Start Time':>15}   {'End Time':>15}")
        for entry in list_mpl_sites:
            print(f"    {entry[0]:<20} {entry[1]:>8.3f} {entry[2]:>9.3f} {entry[3]:>7.3f} {entry[4].strftime('%Y-%m-%d %H:%M'):>20} {entry[5].strftime('%Y-%m-%d %H:%M'):>20}")
        print("    ------------------------------")
    
    # Find matching orbits
    matching_orbits = find_matching_orbits(orbits_list, list_mpl_sites)
    
    if verb:
        print("   Matching sites for each orbit:")
        for orbit, sites in matching_orbits.items():
            print(f"   Orbit: {orbit}")
            print("     Matching sites:")
            for site in sites:
                print(f"     - {site}")
            # print("")
#     breakpoint()
    # Create final output list
    orbit_sites = []
    for orbit, site_names in matching_orbits.items():
        orbit_num = int(orbit[51:56])
        for site_name in site_names:
            site_info = next(site for site in list_mpl_sites if site[0] == site_name)
            location_info = [orbit_num, orbit] + site_info
            orbit_sites.append(location_info)
    
    return orbit_sites,path_2_folder_with_orbits

