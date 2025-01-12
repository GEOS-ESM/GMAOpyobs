"""

   Implements in Python algorithms provided as spreadsheets in the book 

   Fire Science, From Chemistry to Landscape Management
   by Francisco Castro Rego, Penelope Morgan, Paulo Fernandes and Chad Hoffman

"""

import numpy as np

#           (C,H,O)
Species = { (1,0,0)  : 'carbon',
            (0,1,0)  : 'hydrogen',
            (1,4,0)  : 'methane',
            (2,6,0)  : 'ethane', 
            (3,8,0)  : 'propane',
            (8,18,0) : 'octane',
            (1,0,1)  : 'carbon monoxide',
            (6,12,6) : 'glucose',
            (6,9,4)  : 'wood equivalent',
            (1,1.7,0.72)  : 'wood 1',
            (3.4,6.2,2.5) : 'wood 2',
            (10,12,3)     : 'lignin',
          }

# Atomic Weight
Weight = dict ( C=12, H=1, O=16, N=14, O2=16*2, CO=12+16, CO2=12+16*2, H2O=2+16, N2=2*14,CH4=12+4 )

# Low Heat of Combustion (kJ/mol)
#           (C,H,O)
LHC_ = {    (1,0,0)       : 394,
            (0,2,0)       : 242,
            (1,4,0)       : 800,
            (2,6,0)       : 1423,
            (3,8,0)       : 2044, 
            (8,18,0)      : 5104,
            (1,0,1)       : 283,
            (6,12,6)      : 2772,
            (6,9,4)       : 2716,
            (1,1.7,0.72)  : 431,
            (3.4,6.2,2.5) : 1540,
            (10,12,3)     : np.nan,
          }

LHC = dict (   C        = 394,
               H2       = 242,
               CH4      = 800,
               C2H6     = 1423,
               C3H8     = 2044, 
               C8H18    = 5104,
               CO       = 283,
               C6H12O6  = 2772,
               C6H9O4   = 2716,
               WOOD1    = 431,
               WOOD2    = 1540, 
               C10H12O3 = np.nan #  TO DO (Lignin)
         )

# Heat Capacity of Gases (J/mol/K)
#                 298K   1000K
Cp = dict ( N2  = (29.1, 32.7),
            O2  = (29.4, 34.9),
            C   = (8.5,  21.6),
            H2  = (28.8, 30.2),
            CO2 = (40.9, 54.3),
            CO  = (29.1, 33.3),
            CH4 = (35.2, 71.8),
            H2O = (36.0, 41.2),
           )

# Heat Capacity of Solids (J/g/K)
CP = dict ( H2O = 4.2,
            sand = 0.8,
            soil = 1.8,
            wood = 2.0,
          )

# Misc 
L = 2260       # Latent heat of vaporization for water (J/g)
T_i = 598      # Ignition temperature (K)
T_b = 373      # Boiling temperature for water (K)
r_N2_O2 = 3.76 # Molar Nitrogen/Oxygen ratio in air

C, H, O, O2, CO, CO2, CH4, H2O, N, N2 = ('C', 'H', 'O', 'O2', 'CO', 'CO2', 'CH4', 'H2O', 'N', 'N2' )
        
#..........................................................................

class COMBUSTION(object):
 
    def __init__ (self, fuel_C, fuel_H, fuel_O,
                        fuel_lhc, fuel_moist, fuel_temp,
                        air_temp, air_rh, eofr, Verbose=False ):
        """
        Initialize combustion process by specifying fuel and environmental
        conditions.

        TO DO: describe input

        """
        
        # Save inputs
        self.fuel_comp = ( fuel_C, fuel_H, fuel_O)
        self.fuel_lhc = fuel_lhc
        self.fuel_moist = fuel_moist
        self.fuel_temp = fuel_temp
        fuel_mass = fuel_C * Weight[C] + \
                         fuel_H * Weight[H] + \
                         fuel_O * Weight[O]     # molar mass, g/mol
        
        self.air_temp = air_temp
        self.air_rh = air_rh
        self.eofr = eofr

        # Combustion in air
        O2_in_air_H2O = fuel_H / 2
        O2_in_air_CO2 = fuel_C
        O2_in_air_O2  = (fuel_C*2 + fuel_H/2 - fuel_O)/2
        N2_in_air =  O2_in_air_O2 * r_N2_O2
        H2O_in_fuel = fuel_moist * fuel_mass / Weight[H2O]

        if Verbose:
            print('O2_in_air_O2  %5.2f'%O2_in_air_O2)
            print('O2_in_air_CO2 %5.2f'%O2_in_air_CO2)
            print('N2_in_air     %5.2f'%N2_in_air)
            print('O2_in_air_H2O %5.2f'%O2_in_air_H2O)
            print('H2O_in_fuel   %5.2f\n'%H2O_in_fuel)

        # EOFR adjustments
        adj_N2  = N2_in_air    * (eofr - 1)
        adj_O2_ = O2_in_air_O2 * (eofr - 1)
        adj_O2  = max(adj_O2_,0)
        factor  = min(adj_O2_,0)/1.41*2 
        adj_CO  = -0.75 * factor
        adj_CH4 = -0.08 * factor
        adj_C   = -0.17 * factor
        adj_CO2 =         factor
        adj_H2O =  0.16 * factor

        if Verbose:
            print('Adjustment N2       %5.2f'%adj_N2)
            print('Adjustment O2 (raw) %5.2f'%adj_O2_)
            print('Adjustment O2       %5.2f\n'%adj_O2)
            print('Adjustment CO       %5.2f'%adj_CO)
            print('Adjustment CH4      %5.2f'%adj_CH4)
            print('Adjustment C        %5.2f'%adj_C)
            print('Adjustment CO2      %5.2f'%adj_CO2)
            print('Adjustment H2O      %5.2f\n'%adj_H2O)

        # Combustion Efficiency, (C as CO2) / (C in fuel)
        self.ce  = min(1, (O2_in_air_CO2+adj_CO2)/fuel_C )

        # Modified Combustion Efficiency, (C as CO2) / (C in CO2 and CO)
        self.mce = min(1, (O2_in_air_CO2+adj_CO2)/(O2_in_air_CO2+adj_CO2+adj_CO) )

        # Stoichiometric air-to-fuel ratio, (gram / gram)
        self.s_a2f = ( O2_in_air_O2 * Weight[O] * 2   + \
                       N2_in_air    * Weight[N] * 2 ) / fuel_mass

        # Actual air-to-fuel ratio, (gram / gram)
        self.r_a2f = eofr * self.s_a2f 

        if Verbose:
            print('Combustion Efficiency             %5.2f'%self.ce)
            print('Modified Combustion Efficiency    %5.2f'%self.mce)
            print('Stoichiometric air-to-fuel ratio  %5.2f'%self.s_a2f)
            print('Actual air-to-fuel ratio          %5.2f\n'%self.r_a2f)
      
        # Combustion sources, in grams
        self.sources = dict ( FUEL =  fuel_mass,
                              O2  = Weight[O2] * ( O2_in_air_O2 + adj_O2_ ),
                              N2 =  Weight[N2] * ( N2_in_air + adj_N2 ),
                             )

        # Water vapor in the air
        H2O_in_air = (self.sources[O2]+self.sources[N2])*(air_rh/100)*0.00388*np.exp((air_temp-273.15)*0.0656) \
                     / Weight[H2O]
    
        self.sources[H2O] = fuel_moist * fuel_mass + H2O_in_air * Weight[H2O]

        if Verbose:
                print('H2O_in_air              %5.2f\n'%H2O_in_air)
                print('Source FUEL     (grams) %5.2f'%self.sources['FUEL'])
                print('Source OXYGEN   (grams) %5.2f'%self.sources[O2])
                print('Source WATER    (grams) %5.2f'%self.sources[H2O])
                print('Source NITROGEN (grams) %5.2f\n'%self.sources[N2])

        # Combustion products, in grams
        self.products = dict ( CO2  =  Weight[CO2] * (O2_in_air_CO2 + adj_CO2),
                               H2O  =  Weight[H2O] * (O2_in_air_H2O + H2O_in_fuel + H2O_in_air + adj_H2O),
                               N2   =  self.sources[N2],
                               CO   =  Weight[CO]  * adj_CO,
                               CH4  =  Weight[CH4] * adj_CH4,
                               SOOT =  Weight[C]   * adj_C,
                               O2   =  Weight[O2]  * adj_O2,
                               )
        if Verbose:
                print('Combustion Products CARBON DIOXIDE  (grams) %5.2f'%self.products['CO2'])
                print('Combustion Products WATER           (grams) %5.2f'%self.products['H2O'])
                print('Combustion Products NITROGEN        (grams) %5.2f'%self.products['N2'])
                print('Combustion Products CARBON MONOXIDE (grams) %5.2f'%self.products['CO'])
                print('Combustion Products METHANE         (grams) %5.2f'%self.products['CH4'])
                print('Combustion Products SOOT            (grams) %5.2f'%self.products['SOOT'])
                print('Combustion Products OXYGEN          (grams) %5.2f\n'%self.products['O2'])

        # Heat of Combustion (Dry Fuel)
        self.dry_heat = dict ( per_FUEL = fuel_lhc / self.sources['FUEL'],                    # kJ/gram
                               per_O2   = fuel_lhc / self.sources['O2'],                      # kJ/gram
                               per_AIR  = fuel_lhc / (self.sources['O2']+self.sources['N2']), # kJ/gram
                               r_O2f    = self.sources['O2'] / self.sources['FUEL'],          # gram/gram
                              )

        if Verbose:
            print('Combustion Dry Heat per FUEL    %5.2f'%self.dry_heat['per_FUEL'],'kJ/g')
            print('Combustion Dry Heat per OXYGEN  %5.2f'%self.dry_heat['per_O2'],'kJ/g')
            print('Combustion Dry Heat per AIR     %5.2f'%self.dry_heat['per_AIR'],'kJ/g')
            print('Ratio OXYGEN / FUEL             %5.2f\n'%self.dry_heat['r_O2f'])


        # Heat of Pre-ignition
        raise_wood    = (T_i-fuel_temp)*(CP['wood']/1000)                             # kJ/gram
        raise_H2O_liq = (T_b-fuel_temp)*(CP['H2O']*fuel_moist/1000)                   # kJ/gram
        vaporize_H2O  = fuel_moist * L / 1000                                         # kJ/gram
        raise_H2O_gas = fuel_moist * (T_i-373)*(Cp['H2O'][1] / 1000 / Weight['H2O'])  # kJ/gram
 
        self.heat_ignition = raise_wood + raise_H2O_liq + raise_H2O_gas + vaporize_H2O

        if Verbose:
            print('Total Heat of Ignition               %5.2f'%self.heat_ignition,'kJ/g')
            print('To Raise Temperature of wood         %5.2f'%raise_wood,'kJ/g')
            print('To Raise Temperature of Liquid  H2O  %5.2f'%raise_H2O_liq,'kJ/g')
            print('To Vaporize H2O                      %5.2f'%vaporize_H2O,'kJ/g')
            print('To Raise Temperature of Gaseous H2O  %5.2f'%raise_H2O_gas,'kJ/g\n')


        # Heat in Combustible Products (kJ/g)
        heat_combustibles = (adj_CO2 
                                  
if __name__ == "__main__":

    c = COMBUSTION(fuel_C=6, fuel_H=9, fuel_O=4,
                   fuel_lhc=2716, fuel_moist=0.50, fuel_temp=298,
                   air_temp=298, air_rh=90, eofr=0.8,
                   Verbose=True)








    
