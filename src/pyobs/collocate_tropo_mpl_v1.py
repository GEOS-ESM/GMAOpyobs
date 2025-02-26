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

--------------------------------------------
Feb/12/2025 Some logic about this code:
There are two major steps in this code:
1)For selected date, this code ingests all orbits 
available for each day (usually about 15) and then it finds which orbits contain MPL sites.
The purpose is to avoid to open and compute distances for sites in EACH orbit. 
This achieved by comparing the time of observation of the firt and last line in the orbit 
and compare with the expected time of overpass in each site. This is a very fast operation
2)  For list of orbits found in previous step, it locates in 
the orbit's grid the location of each site. Then it outputs the scanline and column number of 
each orbit. Ouput is a text file with this information for each site. 
The rational for step 1 is speed as I/O and distance compuation of large arrays takes time, step 1 
is created to make things faster. However if time is not a concern, step 1 can be skipped 
and the proper inputs  need to be created. Both steps are in separate routines so if step 1 
is removed and all orbits are evaluated, step 2 should be easy to implement.

@author: sgasso
"""
import os,sys,platform
import glob
import h5py as h5
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta

from tropomi_l2_reader import TROPOAER_L2


import f_collocate_1   as fc #needed for finding pixels that contain site 
import f_find_orbits_1 as fo # needed for finding orbits that contain sites
 


if __name__ == "__main__":
    # 
    
    #### User Input , Most of the codes can operate with yyyy and Julian day
    #### But I haven''t tested. I do need to change how julian1 and julian2 are set.
    verb = True
    yyyy0,mm0,dd0=[2023,1,1]; yyyy1,mm1,dd1=[2023,12,31]
    # yyyy0,mm0,dd0=[2023,8,7]; yyyy1,mm1,dd1=[2023,8,7]
    user_string = "FullYear" # additional string to add to output file SiteName_yyyy_<user_string>.txt
    #julian=359 

    #### End user input section 
    #--------------------------------------------------
    if not yyyy0 : sys.exit('Input year most be provided')
    # if not mm and not julian : sys.exit('Either Month and day or Julian day must be provided')

    current_working_directory = os.getcwd()
    current_os    = platform.system()
    computer_name = platform.node()
    if verb :
        print('\nThis code is running from directory :', current_working_directory)
        print(f'This code is executed in computer {computer_name} \n')
    
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
        yyyy = current_date.year
        mm   = current_date.month
        dd   = current_date.day
        # if only julian day is provided, then compute mm and dd to get selected day
        if not mm:
            date = datetime(yyyy, 1, 1) + timedelta(days=julian - 1)
            mm = date.month
            dd = date.day
        SEL_DATE = datetime(yyyy,mm,dd)
        if verb: print(f"Current date: {yyyy}-{mm:02d}-{dd:02d}")
        
        #### Get user defined text file with MPL sites
        pth_site_list  = current_working_directory + '/' + 'list_mpl_sites.txt'
        mpl_sites_list = fo.read_mpl_site_list(pth_site_list)
        SEL_DATE = datetime(yyyy,mm,dd) # 
        ### Now add SEL_DATE to mpl time stamps, this is needed to get the matching orbits 
        for i in mpl_sites_list:
             i[-1]= datetime.combine(SEL_DATE, i[-1])
             i[-2]= datetime.combine(SEL_DATE, i[-2])
         
        #### Now for orbits available for this day, find which ones contain a site 
        matching_orbits_list,current_pth =\
                fo.get_orbit_site_matches(mpl_sites_list,yyyy, mm, dd,verb)
        # orbit_matches_sit = fo.get_orbit_site_matches(mpl_sites_list,yyyy=2023, julian=359,verb=True) 
        if  len(matching_orbits_list)==0 : 
            print('No orbits found that contain the sites. Check ', current_pth)
            print('Skip this day')
        else:
            ##### Now find the line and column number of pixels w/MPL sites inside in selected TROPomi orbits
            output_array = fc.process_matching_orbits(matching_orbits_list, 
                                                 current_pth, SEL_DATE,mpl_sites_list,verb)
            print("")
            # Update combined_results with new data
            for key, value in output_array.items():
                if key in combined_results:
                    # Extend existing list with new values
                    combined_results[key].extend(value)
                else:
                    # Create new key-value pair
                    combined_results[key] = value.copy()        
        print(f"   This day time: {(datetime.now() - start_time).total_seconds():.2f} seconds")

        #### next value for iteration
        current_date += timedelta(days=1)
    
    end_time = datetime.now()
    time_diff = (end_time - start_time).total_seconds()
    print(f"\nLoop Execution time: {time_diff:.2f} seconds")
    ### Save one file per site
    ### Inputs: array_to_save, year_string,user_string
    fc.create_output_files(combined_results, "2023", user_string)
    # fc.create_output_files(combined_results, "2023", "test2")
    
    if verb :
        print('\nThis code is running from directory :', current_working_directory)
        print(f'This code is executed in computer {computer_name} \n')

    print("\nDone!")

