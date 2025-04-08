#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mar/4/2024
This code produces a list of TROPO pixel indexes collocated over MPL sites. 
User must provide dates of interest and an input file with basic MPL site information 
(site_name,lat,lon, altitude , stored in a text file with a header then four columns in each line per site).
 Output: text files for each mpl site containing the date, tropomi_filename, scanline and col
 of the pixel over site.
 This output file can be used to get Tropomi L1b radiances over the surface site and MERRA data.

This version is supersedes *_1.py codes. Version _1 worked well but they required the user to provide an 
estimated local satellite overpass time (to make the search of sites faster). This version 2 does not require a local overpass time at the expense to take longer (essentially this code reads all orbits available for selected day and then loops over each one to find for the site, the v1 code did a prescreen to find the sites). V1 is ~6 times faster than V2 in Calculon. 

--------------------------------------------
Some notes about the logic of this code:
1) This code is essentially a wrapper that sets up data to the collocation of sites with the satellite.
2) For selected date or dates, this code loops in each day and ingests all orbits 
available for each day (usually about 15) and then it finds which orbits contain MPL sites.
This operation carried out by the group of subroutines in module f_collocate . 
3) About Input : the output filename format is <SiteName>_YEAR_<user_string>.txt
where SiteName string extracted from the user define input file, YEAR for this dataset , user_string 
is an arbitrary label defined by user. 
4) Code requires the module f_collocate_2 to be located in the same directory
5) Customize this code at will 

Apr/9/2025 - updated script get_orbit_site_matches to deal with missing folder or missing files 
Apr/8/2025 - added loop over years

@author: sgasso@umd.edu
"""
import os,sys,platform
import glob
import h5py as h5
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta

### Load specific modules. 
#from pyobs.tropomi_l2_reader import TROPOAER_L2 # need to read TROPOMI data
import WIP_f_collocate  as fc #needed for finding pixels that contain site 

if __name__ == "__main__":
    # 
    
    #### User Input , Most of the codes can operate with yyyy and Julian day
    #### But I haven''t tested. I do need to change how julian1 and julian2 are set.
    for iyear in [2019,2020,2021,2022,2023,2024]:
    #for iyear in [2020]:    
        verb = 0
        #yyyy0,mm0,dd0=[2020,12,5]; yyyy1,mm1,dd1=[2020,12,6]
        yyyy0,mm0,dd0=[iyear,1,1]; yyyy1,mm1,dd1=[iyear,12,31]
        user_string = "complete" # additional string to add to output file SiteName_yyyy_<user_string>.txt
        
        output_path='/nobackup/2/MPL/Collocations/'  # location of collocation files
        #julian=359 

        #### End user input section 
        #--------------------------------------------------
        if not yyyy0 : sys.exit('Input year most be provided')
        # if not mm and not julian : sys.exit('Either Month and day or Julian day must be provided')

        current_working_directory = os.getcwd()
        current_os    = platform.system()
        computer_name = platform.node()
        if verb<1 :
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
            start_time_day = datetime.now()
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
            if verb>1: print(f"Current date: {yyyy}-{mm:02d}-{dd:02d}")
            
            #### Get user defined text file with MPL sites
            pth_site_list  = current_working_directory + '/' + 'list_mpl_sites.txt'
            mpl_sites_list = fc.read_mpl_site_list(pth_site_list)
            SEL_DATE = datetime(yyyy,mm,dd) # 

             
            #### Now get list of files and the path to corresponding folder
            matching_orbits_list,current_pth =\
                    fc.get_orbit_site_matches(mpl_sites_list,yyyy, mm, dd,verb)

            ##### Now find the line and column number of pixels w/MPL sites inside in selected TROPomi orbits
            output_array = fc.process_matching_orbits(matching_orbits_list, 
                                                 current_pth, SEL_DATE,mpl_sites_list,verb)
            #print("")
            # Update combined_results with new data
            for key, value in output_array.items():
                if key in combined_results:
                    # Extend existing list with new values
                    combined_results[key].extend(value)
                else:
                    # Create new key-value pair
                    combined_results[key] = value.copy()        
            print(f"       This day time: {(datetime.now() - start_time_day).total_seconds():.2f} seconds")

            #### next value for iteration
            current_date += timedelta(days=1)
        
        #### Done with loop over each day
        end_time = datetime.now()
        time_diff = (end_time - start_time).total_seconds()
        print(f"\nLoop Execution time: {time_diff:.2f} seconds")
        
        
        ### Save one file per site
        ### Inputs: array_to_save, year_string,user_string
        fc.create_output_files(output_path,combined_results, str(yyyy0), user_string)
        
        if verb >-1 :
            print('\nThis code is running from directory :', current_working_directory)
            print(f'This code is executed in computer {computer_name} \n')
            print('Files are saved in directory ', output_path)

        print("\nDone!")

