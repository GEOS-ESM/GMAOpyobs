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

FUEL, HC, PM25, PM10, PM = ('FUEL', 'HC', 'PM25', 'PM10', 'PM')

#..........................................................................

class COMBUSTION(object):
 
    def __init__ (self, fuel_C, fuel_H, fuel_O,
                        fuel_lhc, fuel_moist, 
                        air_temp, air_rh, eofr,
                        fuel_temp=None, Verbose=False ):
        """
        Initialize combustion process by specifying fuel and environmental
        conditions.

        fuel_C, fuel_H, Fuel_O  ---  fuel composition: carbon, hydrogen, oxygen
        fuel_lhc                ---  low heat combustion of fuel (kJ/mol)
        fuel_moist              ---  fuel moisture content [0-1]
        air_temp                ---  air temperature (K)
        air_rh                  ---  air relative humidity [0-100]
        eofr                    ---  equivalent oxygen to fuel ratio [0-1]
        fuel_temp               ---  fuel temperature (K). If None, same as
                                     air_temp
        Verbose                 ---  prints out intermediate results.
        

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
        self.O2_in_air_H2O = fuel_H / 2
        self.O2_in_air_CO2 = fuel_C
        self.O2_in_air_O2  = (fuel_C*2 + fuel_H/2 - fuel_O)/2
        self.N2_in_air =  self.O2_in_air_O2 * r_N2_O2
        self.H2O_in_fuel = fuel_moist * fuel_mass / Weight[H2O]

        if Verbose:
            print('O2_in_air_O2  %5.2f'%self.O2_in_air_O2)
            print('O2_in_air_CO2 %5.2f'%self.O2_in_air_CO2)
            print('N2_in_air     %5.2f'%self.N2_in_air)
            print('O2_in_air_H2O %5.2f'%self.O2_in_air_H2O)
            print('H2O_in_fuel   %5.2f\n'%self.H2O_in_fuel)

        # EOFR adjustments
        self.adj_N2  = self.N2_in_air    * (eofr - 1)
        self.adj_O2_ = self.O2_in_air_O2 * (eofr - 1)
        self.adj_O2  = max(self.adj_O2_,0)
        factor  = min(self.adj_O2_,0)/1.41*2 
        self.adj_CO  = -0.75 * factor
        self.adj_CH4 = -0.08 * factor
        self.adj_C   = -0.17 * factor
        self.adj_CO2 =         factor
        self.adj_H2O =  0.16 * factor

        if Verbose:
            print('Adjustment N2       %5.2f'%self.adj_N2)
            print('Adjustment O2 (raw) %5.2f'%self.adj_O2_)
            print('Adjustment O2       %5.2f\n'%self.adj_O2)
            print('Adjustment CO       %5.2f'%self.adj_CO)
            print('Adjustment CH4      %5.2f'%self.adj_CH4)
            print('Adjustment C        %5.2f'%self.adj_C)
            print('Adjustment CO2      %5.2f'%self.adj_CO2)
            print('Adjustment H2O      %5.2f\n'%self.adj_H2O)

        # Combustion Efficiency, (C as CO2) / (C in fuel)
        self.ce  = min(1, (self.O2_in_air_CO2+self.adj_CO2)/fuel_C )

        # Modified Combustion Efficiency, (C as CO2) / (C in CO2 and CO)
        self.mce = min(1, (self.O2_in_air_CO2+self.adj_CO2)/(self.O2_in_air_CO2+self.adj_CO2+self.adj_CO) )

        # Stoichiometric air-to-fuel ratio, (gram / gram)
        #s_a2f = ( self.O2_in_air_O2 * Weight[O] * 2   + \
        #               N2_in_air    * Weight[N] * 2 ) / fuel_mass

        # Actual air-to-fuel ratio, (gram / gram)
        #r_a2f = eofr * s_a2f 

        if Verbose:
            print('Combustion Efficiency             %5.2f'%self.ce)
            print('Modified Combustion Efficiency    %5.2f'%self.mce)
            #print('Stoichiometric air-to-fuel ratio  %5.2f'%s_a2f)
            #print('Actual air-to-fuel ratio          %5.2f\n'%r_a2f)
      
        # Combustion sources, in grams
        self.sources = dict ( FUEL =  fuel_mass,
                              O2  = Weight[O2] * ( self.O2_in_air_O2 + self.adj_O2_ ),
                              N2 =  Weight[N2] * ( self.N2_in_air + self.adj_N2 ),
                             )

        # Water vapor in the air
        self.H2O_in_air = (self.sources[O2]+self.sources[N2])*(air_rh/100)*0.00388*np.exp((air_temp-273.15)*0.0656) \
                     / Weight[H2O]
    
        self.sources[H2O] = fuel_moist * fuel_mass + self.H2O_in_air * Weight[H2O]

        if Verbose:
                print('H2O_in_air              %5.2f\n'%self.H2O_in_air)
                print('Source FUEL     (grams) %5.2f'%self.sources[FUEL])
                print('Source OXYGEN   (grams) %5.2f'%self.sources[O2])
                print('Source WATER    (grams) %5.2f'%self.sources[H2O])
                print('Source NITROGEN (grams) %5.2f\n'%self.sources[N2])

        # Combustion products, in grams
        self.products = dict ( CO2  =  Weight[CO2] * (self.O2_in_air_CO2 + self.adj_CO2),
                               H2O  =  Weight[H2O] * (self.O2_in_air_H2O + self.H2O_in_fuel + self.H2O_in_air + self.adj_H2O),
                               N2   =  self.sources[N2],
                               CO   =  Weight[CO]  * self.adj_CO,
                               CH4  =  Weight[CH4] * self.adj_CH4,
                               SOOT =  Weight[C]   * self.adj_C,
                               O2   =  Weight[O2]  * self.adj_O2,
                               )
        if Verbose:
                print('Combustion Products CARBON DIOXIDE  (grams) %5.2f'%self.products['CO2'])
                print('Combustion Products WATER           (grams) %5.2f'%self.products['H2O'])
                print('Combustion Products NITROGEN        (grams) %5.2f'%self.products['N2'])
                print('Combustion Products CARBON MONOXIDE (grams) %5.2f'%self.products['CO'])
                print('Combustion Products METHANE         (grams) %5.2f'%self.products['CH4'])
                print('Combustion Products SOOT            (grams) %5.2f'%self.products['SOOT'])
                print('Combustion Products OXYGEN          (grams) %5.2f\n'%self.products['O2'])


    def getT(self,Verbose=False):
        """
        Calculates Flame temperature.

            T = c.getT()
        
        """

        fuel_lhc   = self.fuel_lhc 
        fuel_temp  = self.fuel_temp 
        fuel_moist = self.fuel_moist
               
        # Heat of Combustion (Dry Fuel)
        dry_heat = dict ( per_FUEL = fuel_lhc / self.sources['FUEL'],                    # kJ/gram
                               per_O2   = fuel_lhc / self.sources['O2'],                      # kJ/gram
                               per_AIR  = fuel_lhc / (self.sources['O2']+self.sources['N2']), # kJ/gram
                               r_O2f    = self.sources['O2'] / self.sources['FUEL'],          # gram/gram
                              )

        if Verbose:
            print('Combustion Dry Heat per FUEL    %5.2f'%dry_heat['per_FUEL'],'kJ/g')
            print('Combustion Dry Heat per OXYGEN  %5.2f'%dry_heat['per_O2'],'kJ/g')
            print('Combustion Dry Heat per AIR     %5.2f'%dry_heat['per_AIR'],'kJ/g')
            print('Ratio OXYGEN / FUEL             %5.2f\n'%dry_heat['r_O2f'])


        # Heat of Pre-ignition
        raise_wood    = (T_i-fuel_temp)*(CP['wood']/1000)                             # kJ/gram
        raise_H2O_liq = (T_b-fuel_temp)*(CP['H2O']*fuel_moist/1000)                   # kJ/gram
        vaporize_H2O  = fuel_moist * L / 1000                                         # kJ/gram
        raise_H2O_gas = fuel_moist * (T_i-373)*(Cp['H2O'][1] / 1000 / Weight['H2O'])  # kJ/gram
 
        heat_ignition = raise_wood + raise_H2O_liq + raise_H2O_gas + vaporize_H2O

        if Verbose:
            print('Total Heat of Ignition               %5.2f'%heat_ignition,'kJ/g')
            print('To Raise Temperature of wood         %5.2f'%raise_wood,'kJ/g')
            print('To Raise Temperature of Liquid  H2O  %5.2f'%raise_H2O_liq,'kJ/g')
            print('To Vaporize H2O                      %5.2f'%vaporize_H2O,'kJ/g')
            print('To Raise Temperature of Gaseous H2O  %5.2f'%raise_H2O_gas,'kJ/g\n')


        # Heat in Combustible Products (kJ/g)
        heat_combustibles = (self.adj_CO * LHC[CO] + self.adj_CH4 * LHC[CH4] + self.adj_C * LHC[C] ) / self.sources['FUEL']

        # ADJUSTED HEAT YIELD (kJ/gram)
        heat_yield = dry_heat['per_FUEL'] - heat_ignition - heat_combustibles
        heat_yield_per_mol = heat_yield * self.sources['FUEL']
        
        if Verbose:
            print('Heat in CO, CH4 and C  %7.2f'%heat_combustibles,'kJ/g')
            print('Adjusted Heat Yield    %7.2f'%heat_yield,'kJ/g')
            print('Adjusted Heat Yield    %7.2f'%heat_yield_per_mol,'kJ/mol\n')


        # Heat capacity of the mixture of products (per mol of fuel)
        heat_capacity =   self.O2_in_air_CO2 * Cp[CO2][1] + \
                         (self.O2_in_air_H2O+self.H2O_in_air+self.H2O_in_fuel) * Cp[H2O][1] + \
                         (self.N2_in_air+self.adj_N2) * Cp[N2][1] + \
                         self.adj_O2 * Cp[O2][1]
        delta_T = heat_yield_per_mol / heat_capacity * 1000
        flame_T = fuel_temp + delta_T
                         
        if Verbose:
            print('Estimated temperature of gases in the flame  %5.2f'%flame_T,'K')
            print('Heat capacity of the mixture of products     %5.2f'%heat_capacity,'J/K')
            print('Estimated increase in temperature            %5.2f'%delta_T,'K\n')

        return flame_T
            
    def getEFs(self,Verbose=False):
        """
        Return emission factors for key species.

        EF = c.getEFs()
        
        """
            
        # Emission Factors
        EF = dict ( H2O  = 1000 * self.products[H2O] / self.sources[FUEL],
                    CO2  = 1000 * self.products[CO2] / self.sources[FUEL],
                    CO   = 1000 * self.products[CO] / self.sources[FUEL],
                    CH4  = (1000-408) * self.products[CH4] / self.sources[FUEL],
                    HC   =       408  * self.products[CH4] / self.sources[FUEL],
                    PM   = 1000 * self.products['SOOT'] / self.sources[FUEL]
                   )
        EF['PM25'] = 0.692 * EF['PM']
        EF['PM10'] = 0.816 * EF['PM']

        if Verbose:
            print('Emission Factor H2O    %7.2f'%EF[H2O], 'g/Kg')
            print('Emission Factor CO2    %7.2f'%EF[CO2], 'g/Kg')
            print('Emission Factor CO     %7.2f'%EF[CO],  'g/Kg')
            print('Emission Factor CH4    %7.2f'%EF[CH4], 'g/Kg')
            print('Emission Factor HC     %7.2f'%EF[HC],  'g/Kg')
            print('Emission Factor PM2.5  %7.2f'%EF[PM25],'g/Kg')
            print('Emission Factor PM10   %7.2f'%EF[PM10],'g/Kg')
            print('Emission Factor PM     %7.2f'%EF[PM],  'g/Kg')
  
        return EF
            
if __name__ == "__main__":
 
    c = COMBUSTION(fuel_C=6, fuel_H=9, fuel_O=4,
                   fuel_lhc=2716, fuel_moist=0.50, fuel_temp=298,
                   air_temp=298, air_rh=90, eofr=0.8,
                   Verbose=True)
    
    T = c.getT(Verbose=True)
    
    c.getEFs(Verbose=True)








    
