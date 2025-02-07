#   created 8/28/2024
import os, sys, time,platform
from datetime import date, datetime, timedelta
from tropomi_l2_reader import TROPOAER_L2

###### -------- test area
if __name__ == "__main__":
    current_working_directory = os.getcwd()
    print('\nThis code is running from directory :', current_working_directory)
    # Get user inputs
    now_os   = platform.system().lower()
    now_computer = platform.node()
    print(f'This code is executed in computer {now_computer} \n')

    # now_os=sys.platform
   # if current_os=='win32': #in PC office
      # pthin="D:\Satellite\Tropomi\Level2\/226\/"
   # elif current_os=='darwin': #in Laptop
      # # pthin='/Users/sgasso/'
      # pthin='/Volumes/ExtData1/SatData/Tropomi/Level2/'
   # elif current_os == 'linux': # in Calculon
      # pthin='/nobackup/TROPOMAER/2019/226/'
   # else:
       # print('Current operating system no recognized.')
       # print('Cannot set path to Level1 files. Terminate Run')
       # sys.exit()

    # Your existing get_path function remains the same
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

 
    filename_tro='TROPOMI-Sentinel-5P_L2-TROPOMAER_2023m1225t132102-o32123_v01-2023m1231t105113.nc'
    syn_time = datetime(2023,12,15,13,1,2)
    Files = base_path+filename_tro

    # Extract only GEO data
    TROP = TROPOAER_L2(Files,GEO=True,Verbose=0)
    ## extract defaule list of variables , see list in tropomaer.py
    # TROP2 = TROPOAER_L2(Files,Verbose=1)


####  Check date
    # Both conversions in one-liners
    date1    = (datetime(2010, 1, 1) + timedelta(seconds=int(TROP.time[0])))
    date1_str= date1.strftime('%Y-%m-%dT%H:%M:%S')
    date2     = date1 + timedelta(milliseconds=int(TROP.dtime[0]))
    date2_str = date2.strftime('%Y-%m-%dT%H:%M:%S')
    # Print results
    print(f"Reference day: {date1}")
    print(f"Time first scanline : {date2}")