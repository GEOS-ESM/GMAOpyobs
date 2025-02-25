#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code produces a list of TROPO pixel indexes collocated over MPL sites. 
User must supply data and an input file with basic MPL site information (name,lat,lon, and estimated
 time of tropomi overpass by the site).
 Output is/are text files for each mpl site containing the date, tropomi_filename, scanling and col
 of the pixel over site.
 This output file can be used to get Tropomi L1b radiances over site and MERRA data.

Feb/13/2025 - This code works 
Feb/12/2025 - _2.py : Using _1.py as template, this version is loops over a range of days provided by user

Feb/12/2025 created code. This code is the merge of two existing codes:

test_create_list_orbits_with_sites_2.py: for selected date, this code ingests all orbits 
available for each day (usually about 15) and then it finds which orbits contains MPL sites.
The purpose is to avoide to open and scan for sites in each orbit. The code uses estimate 
time of satellite overpass at the site and compares with the time of start and end of 
sensor scan. This is a fast operation. This operation takes more time 
than just selecting orbits already known to contain the sites.

  
test_find_sites_in_orbit_5.py : For list of orbits found in previous code, it locates in 
the orbit's grid the location of each site. Then it outputs the scanline and column number of 
each orbit. Ouput is a text file with this information for each site. 
 

@author: sgasso
"""
import os,sys,platform
import glob
import h5py as h5
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta

from tropomi_l2_reader import TROPOAER_L2


import test_f_collocate_1   as fc #needed for finding pixels that contain sites 
import test_f_find_orbits_1 as fo # needed for finding orbits that contain sites
 


if __name__ == "__main__":
    # REF_DATE=datetime(2010, 1, 1) # Reference date for Time Array in Tropomi SDS variable
    current_working_directory = os.getcwd()
    current_os    = platform.system()
    computer_name = platform.node()
#     current_pth,pth_fig = fc.get_path(current_os.lower(),computer_name)
    
    #### User Input
    verb = True
#     yyyy0 =2023 ;mm0=12; dd0=24
#     yyyy1 =2023 ;mm1=12; dd1=25
    
    yyyy0,mm0,dd0=[2023,12,1]
    yyyy1,mm1,dd1=[2023,12,1]
    
    #yyyy =2023 ;mm=12; dd=25
    #julian=359
    
    #### # Convert to datetime objects
    start_date = datetime(yyyy0, mm0, dd0)
    end_date   = datetime(yyyy1, mm1, dd1)
    
    # Initialize the combined dictionary to store all results
    combined_results = {}
    # Loop through each day
    current_date = start_date
    start_time = datetime.now() # start timing the loop
    while current_date <= end_date:
        # Extract year, month, day into separate variables
        # to match the previous code
        yyyy = current_date.year
        mm = current_date.month
        dd = current_date.day
        #### previous code starts here.
        #### Some checks 
        if not yyyy : sys.exit('Input year most be provided')
        if not mm and not julian : sys.exit('Either Month and day or Julian day must be provided')
        if verb :
            print('\nThis code is running from directory :', current_working_directory)
            print(f'This code is executed in computer {computer_name} \n')
    
        # if only julian day is provided, then compute mm and dd to get selected day
        if not mm:
            date = datetime(yyyy, 1, 1) + timedelta(days=julian - 1)
            mm = date.month
            dd = date.day
        SEL_DATE = datetime(yyyy,mm,dd)
        if verb: print(f"   Current date: {yyyy}-{mm:02d}-{dd:02d}")
        
        #### Get user defined text file with MPL sites to process
        pth_site_list  = current_working_directory + '/' + 'list_mpl_sites.txt'
        mpl_sites_list = fo.read_mpl_site_list(pth_site_list)
        #### The next list was created by running test_read_level2_geo_2.py
        #### Using a local dictionary now but eventually this will be the output of test_read_level2_geo_2.py
        SEL_DATE = datetime(yyyy,mm,dd) # 
        ### Now add SEL_DATE to mpl time stamps, this is needed to get the matching orbits 
        for i in mpl_sites_list:
             i[-1]= datetime.combine(SEL_DATE, i[-1])
             i[-2]= datetime.combine(SEL_DATE, i[-2])
         
    
        ######  
        matching_orbits_list,current_pth =\
                fo.get_orbit_site_matches(mpl_sites_list,yyyy, mm, dd,verb)
    #     sys.exit()
    #     orbit_matches = get_orbit_site_matches(mpl_sites_list,yyyy=2023, julian=359,verb=True) 
                             
        ##### Now find the locations of the MPL sites inside the selected TROPomi orbits
        output_array = fc.process_matching_orbits(matching_orbits_list, 
                                             current_pth, SEL_DATE,mpl_sites_list,verb)
        
        # Update combined_results with new data
        for key, value in output_array.items():
            if key in combined_results:
                # Extend existing list with new values
                combined_results[key].extend(value)
            else:
                # Create new key-value pair
                combined_results[key] = value.copy()        
        
        
        print(f"   Iteration: {(datetime.now() - start_time).total_seconds():.2f} seconds")

        #### next value for iteration
        current_date += timedelta(days=1)
    
    end_time = datetime.now()
    time_diff = (end_time - start_time).total_seconds()
    print(f"\nLoop Execution time: {time_diff:.2f} seconds")
    ### Save one file per site
    ### Inputs: array_to_save, year_string,user_string
#     fc.create_output_files(output_array, "2023", "test")
    fc.create_output_files(combined_results, "2023", "test2")
    print("\nDone!")

