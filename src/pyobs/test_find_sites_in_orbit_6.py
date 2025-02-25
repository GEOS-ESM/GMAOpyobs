#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Feb/12/2026 Same as _v5 but all modules are in an external file

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

import test_f_find_orbits_1 as ff  # # needed for finding pixels that contain sites

if __name__ == "__main__":
    REF_DATE=datetime(2010, 1, 1) # Reference date for Time Array in Tropomi SDS variable
    current_working_directory = os.getcwd()
    print('\nThis code is running from directory :', current_working_directory)
    # Get user inputs
    current_os    = platform.system()
    computer_name = platform.node()
    print(f'This code is executed in computer {computer_name} \n')
    current_pth,pth_fig = ff.get_path(current_os.lower(),computer_name)
    
    #### Get user defined text file with MPL sites to process
    pth_site_list  = current_working_directory + '/' + 'list_mpl_sites.txt'
    mpl_sites_list = ff.read_mpl_site_list(pth_site_list)
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
    output_array = ff.process_matching_orbits(matching_orbits_list, current_working_directory, 
                                         current_pth, SEL_DATE, REF_DATE,mpl_sites_list)
    
    ##### Next save outfile one per site
    ### Inputs: array_to_save, year_string,user_string
    ff.create_output_files(output_array, "2023", "test")

# To do, finish editing read_mpl_site_list module, remove input SEL_DATE ,read data as time only and add 
# SEL_DATE later. 