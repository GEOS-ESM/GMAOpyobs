"""
  Purpose:
    Contains information for speeds altitudes and names for each platform to be used.
    Created from platform.txt in Sam LeBlanc's Moving Lines package.
  
  Dictionary Keys:
    Platform: name of the platform, used internally
    names: [list of string] names to be used to dentify platform when selecting 'new flight path' button
    max_alt: [meters] the maximum altitude of the platform, also used as cruising altitude
    base_speed: [m/s] what the speed would be at Sea level (calculated from offset of slope in altitude vs. speed)
    speed_per_alt:  [m/s] the slope of speed increase, per meter of altitude
    max_speed: [m/s] maximum speed to reach
    max_speed_alt: [meters] altitude at which the max speed is reached
    descent_speed_decrease: [m/s] value of speed decrease (positive) when platform is descending
    climb_vert_speed: [m/s] constant vertical speed for climbs
    descent_vert_speed: [m/s] constant vertical speed for descents (must be negative)
    alt_for_variable_vert_speed: [m/s] altitude above which the vertical speed is variable
    vert_speed_base: [m/s] intercept of vertical speed vs. alt, when vertical speed is variable
    vert_speed_per_alt: [m/s per m] slope of vertical speed vs. alt, when vetical speed is variable
    rate_of_turn: [degree per second] default rate of turn of the platform, if set to None, then uses
                  the turn_bank_angle to calculate turns
    turn_bank_angle: [degrees] the bank angle in turns, used to calculate the rate_of_turn
    warning: [True or False] makes a warning appear to th user that the speeds are not well established.
             Should be set to false when correct speeds are put in. 
  Comments:
    P3 speeds calculated from TRACE-P flight from Steven Howell
    ER-2 speeds calculated from SEAC4RS from Samuel LeBlanc
    C130 speeds from ARISE from Samuel LeBlanc
    DC8 speeds from SEAC4RS from Samuel LeBlanc
    BAE speeds unknown

  Modification history:
    Written: Samuel LeBlanc, NASA Ames, from Santa Cruz, 2016-07-10
    Modified: Converted to a simple Python module defining a master platform dictionary.

"""

platform = dict(
P3 = {'Platform':'p3','names':['p3','P3','P-3','p-3','p 3','P 3'],
    'max_alt':7000.0,'base_speed':110.0,'speed_per_alt':0.00925,
    'max_speed':160.0,'max_speed_alt':5400.0,'descent_speed_decrease':15.0,
    'climb_vert_speed':7.5,'descent_vert_speed':-6.5,'alt_for_variable_vert_speed':5400.0,
    'vert_speed_base':4.5,'vert_speed_per_alt':7e-05,
    'rate_of_turn':None,'turn_bank_angle':15.0,
    'warning':False,'pilot_format':'DD MM'},
ER2 = {'Platform':'er2','names':['er2','ER2','ER-2','er-2','ER 2','er 2','Er','eR'],
    'max_alt':19000.0,'base_speed':70.0,'speed_per_alt':0.0071,
    'max_speed':300.0,'max_speed_alt':30000.0,'descent_speed_decrease':0.0,
    'climb_vert_speed':10.0,'descent_vert_speed':-10.0,'alt_for_variable_vert_speed':0.0,
    'vert_speed_base':24.0,'vert_speed_per_alt':0.0011,
    'rate_of_turn':None,'turn_bank_angle':15.0,
    'warning':False},
DC8 = {'Platform':'dc8','names':['dc8','DC8','DC-8','dc-8','DC 8','dc 8','Dc','dC'],
    'max_alt':13000.0,'base_speed':130.0,'speed_per_alt':0.0075,
    'mean_speed': 136,
    'max_speed':175.0,'max_speed_alt':6000.0,'descent_speed_decrease':15.0,
    'climb_vert_speed':15.0,'descent_vert_speed':-10.0,'alt_for_variable_vert_speed':0.0,
    'vert_speed_base':15.0,'vert_speed_per_alt':0.001,
    'rate_of_turn':None,'turn_bank_angle':15.0,
    'warning':False},
C130 = {'Platform':'c130','names':['c130','C130','C-130','c-130','C 130','c 130'],
    'max_alt':7500.0,'base_speed':130.0,'speed_per_alt':0.0075,
    'max_speed':175.0,'max_speed_alt':6000.0,'descent_speed_decrease':15.0,
    'climb_vert_speed':10.0,'descent_vert_speed':-10.0,'alt_for_variable_vert_speed':0.0,
    'vert_speed_base':10.0,'vert_speed_per_alt':0.001,
    'rate_of_turn':None,'turn_bank_angle':20.0,
    'warning':False},
BAE146 = {'Platform':'bae146','names':['bae','BAE','146','BAe','Bae'],
    'max_alt':8500.0,'base_speed':130.0,'speed_per_alt':0.002,
    'max_speed':150.0,'max_speed_alt':8000.0,'descent_speed_decrease':15.0,
    'climb_vert_speed':5.0,'descent_vert_speed':-5.0,'alt_for_variable_vert_speed':8000.0,
    'vert_speed_base':4.5,'vert_speed_per_alt':7e-05,
    'rate_of_turn':None,'turn_bank_angle':20.0,
    'warning':True},
AJAX = {'Platform':'ajax','names':['ajax','Ajax','AJAX','alphajet','alpha','alpha-jet'],
    'max_alt':9500.0,'base_speed':160.0,'speed_per_alt':0.09,
    'max_speed':250.0,'max_speed_alt':9000.0,'descent_speed_decrease':5.0,
    'climb_vert_speed':5.0,'descent_vert_speed':-5.0,'alt_for_variable_vert_speed':8000.0,
    'vert_speed_base':4.5,'vert_speed_per_alt':7e-05,
    'rate_of_turn':None,'turn_bank_angle':25.0,
    'warning':True},
FALCON = {'Platform':'Falcon','names':['falcon','Falcon','HU25','HU-25'],
    'max_alt':4500.0,'base_speed':115,'speed_per_alt':0.005,
    'max_speed':121,'max_speed_alt':4000.0,'descent_speed_decrease':5.0,
    'climb_vert_speed':5.0,'descent_vert_speed':-5.0,'alt_for_variable_vert_speed':8000.0,
    'vert_speed_base':4.5,'vert_speed_per_alt':7e-05,
    'rate_of_turn':None,'turn_bank_angle':25.0,
    'warning':True},
KINGAIR = {'Platform':'KingAir','names':['KingAir','B200','UC12'],
    'max_alt':8000.0,'base_speed':95,'speed_per_alt':0.00325,
    'max_speed':145,'max_speed_alt':8000.0,'descent_speed_decrease':5.0,
    'climb_vert_speed':5.0,'descent_vert_speed':-5.0,'alt_for_variable_vert_speed':8000.0,
    'vert_speed_base':4.5,'vert_speed_per_alt':7e-05,
    'rate_of_turn':None,'turn_bank_angle':25.0,
    'warning':True},
G3 = {'Platform':'g3','names':['g3','G3','G-3','g-3','g 3','G 3','GIII','gIII'],
    'max_alt':12000.0,'base_speed':150.0,'speed_per_alt':0.017,
    'max_speed':225.0,'max_speed_alt':8500.0,'descent_speed_decrease':15.0,
    'climb_vert_speed':10,'descent_vert_speed':-5.0,'alt_for_variable_vert_speed':8000.0,
    'vert_speed_base':6,'vert_speed_per_alt':7e-05,
    'rate_of_turn':None,'turn_bank_angle':10.0,
    'warning':False,'pilot_format':'DD MM'},
LEARJET = {'Platform':'LearJet','names':['lear','Lear','LJ','LearJet','Lear Jet','learjet','Learjet','Lear jet'],
    'max_alt':12000.0,'base_speed':150.0,'speed_per_alt':0.017,
    'max_speed':235.0,'max_speed_alt':8500.0,'descent_speed_decrease':15.0,
    'climb_vert_speed':10,'descent_vert_speed':-5.0,'alt_for_variable_vert_speed':8000.0,
    'vert_speed_base':7,'vert_speed_per_alt':7e-05,
    'rate_of_turn':None,'turn_bank_angle':10.0,
    'warning':False,'pilot_format':'DD MM'},
TWINOTTER = {'Platform':'TwinOtter','names':['twin','Twin','TWIN','TO','Twin Otter','Twin otter','Twin-Otter','Twin-otter','NPS','T-O','Otter','otter','twinotter','TwinOtter','CIRPAS','cirpas'],
    'max_alt':5400.0,'base_speed':47.0,'speed_per_alt':0.019,
    'max_speed':80.0,'max_speed_alt':2000.0,'descent_speed_decrease':7.0,
    'climb_vert_speed':3,'descent_vert_speed':-3.0,'alt_for_variable_vert_speed':1000.0,
    'vert_speed_base':3,'vert_speed_per_alt':7e-05,
    'rate_of_turn':None,'turn_bank_angle':20.0,
    'warning':False,'pilot_format':'DD MM'},
SHEARWATER = {'Platform':'shearwater','names':['shearwater','R/V','SW','RV','Shear Water'],
    'max_alt':0.0,'base_speed':10.2,'speed_per_alt':0.019,
    'max_speed':13.9,'max_speed_alt':0.0,'descent_speed_decrease':0.0,
    'climb_vert_speed':0,'descent_vert_speed':0.0,'alt_for_variable_vert_speed':0.0,
    'vert_speed_base':0,'vert_speed_per_alt':0.0,
    'rate_of_turn':None,'turn_bank_angle':20.0,
    'warning':False,'pilot_format':'DD MM'}
    )
