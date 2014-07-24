# this file contains constants and info on atmospheric gases

#solar constant:
sol_const = 1361e0 # [W/m^2]
# solar radius
sol_radius = 6.955e8 # [m]
# 1 AU in meter
AU = 1.49597870700e11 # [m]
#proton mass
mu = 1.67262158e-27 # [kg]
# gas constant
R_gas = 8.314e0 # J/K/mol
R_Earth = 6371000e0 # m
sigma = 5.67e-8 # J m^-2 s^-1 K^-4
grav_const = 6.67e-11 # N m^2 / kg^2
M_Earth = 5.97e24 # kg


# heat capacity, Cp in [J/mol/K] of atmospheric gases
# data is from the NIST-JANAF thermochemical tables
# http://kinetics.nist.gov/janaf/
# the data tables can also be found in codepath/thermodata/

TCp1  = [100.  , 200.  , 298.15, 300.  , 400.  , 500.  , 600.  , 700.  , 800.  , 900.  , 1000. ]
TCp2  = [100.  , 200.  , 250.  , 298.15, 300.  , 350.  , 400.  , 450.  , 500.  , 600.  , 700.  , 800.  , 900.  , 1000. ]
CpH2O = [33.299, 33.349, 33.590, 33.596, 34.262, 35.226, 36.325, 37.495, 38.721, 39.987, 41.268]
CpCO2 = [29.208, 32.359, 37.129, 37.221, 41.325, 44.627, 47.321, 49.564, 51.434, 52.999, 54.308]
CpCH4 = [33.258, 33.473, 34.216, 35.639, 35.708, 37.874, 40.500, 43.374, 46.342, 52.227, 57.794, 62.932, 67.601, 71.795]
CpNH3 = [33.284, 33.757, 35.652, 35.701, 38.716, 42.048, 45.293, 48.354, 51.235, 53.948, 56.491]
CpO3  = [33.292, 35.058, 39.238, 39.330, 43.744, 47.262, 49.857, 51.752, 53.154, 54.208, 55.024]
CpN2  = [29.104, 29.107, 29.111, 29.124, 29.125, 29.165, 29.249, 29.387, 29.580, 30.110, 30.754, 31.433, 32.090, 32.697]
CpO2  = [29.106, 29.126, 29.201, 29.376, 29.385, 29.694, 30.106, 30.584, 31.091, 32.090, 32.981, 33.733, 34.355, 34.870]
CpH2  = [28.154, 27.447, 28.344, 28.836, 28.849, 29.081, 29.181, 29.229, 29.260, 29.327, 29.441, 29.624, 29.881, 30.205]
CpHe  = [20.786, 20.786, 20.786, 20.786, 20.786, 20.786, 20.786, 20.786, 20.786, 20.786, 20.786, 20.786, 20.786, 20.786]

# moldata contains info on molecules.
# the order of data is:
# mol. name, molecule number in HITRAN if available (0 if not), molecular weight,
# partition fn (linear or non-linear molecule), adiabatic index (gamma),
# array of Cp values as defined on a temperatuer grid
# todo: the adiabatic index somewhat varies with temperature especially for complex molecules.

mol={}
mol['H2O'] = ['H2O', 1  , 18., 'NonLin', 1.33 , TCp1, CpH2O ]
mol['CO2'] = ['CO2', 2  , 44., 'Lin'   , 1.3  , TCp1, CpCO2 ]
mol['CH4'] = ['CH4', 6  , 16., 'NonLin', 1.3  , TCp2, CpCH4 ]
mol['NH3'] = ['NH3', 11 , 17., 'NonLin', 1.31 , TCp1, CpNH3 ]
mol['O3' ] = ['O3' , 3  , 48., 'NonLin', 1.3  , TCp1, CpO3  ]
mol['N2' ] = ['N2' , 0  , 28., 'Lin'   , 1.4  , TCp2, CpN2  ]
mol['O2' ] = ['O2' , 0  , 32., 'Lin'   , 1.4  , TCp2, CpO2  ]
mol['H2' ] = ['H2' , 0  , 2. , ''      , 1.38 , TCp2, CpH2  ]
mol['He' ] = ['He' , 0  , 4. , ''      , 1.66 , TCp2, CpHe  ]

# refracive indices and King correction factors of gases. 

# from Sneep et al, 2005
def refind_N2(wavenr):
   import numpy as np
   refind = np.zeros(len(wavenr))
   refind[(wavenr < 21360e0)] = 1e0 + (6498.2e0 + 307.43305e12 / (14.4e9 - wavenr[(wavenr < 21360e0)]**2e0))*1e-8
   refind[(wavenr >= 21360e0)] = 1e0 + (5677.465e0 + 318.81874e12 / (14.4e9 - wavenr[(wavenr >= 21360e0)]**2e0))*1e-8
   return refind
   
def King_N2(wavenr):
   return 1.034 + 3.17e-12*(wavenr)**2e0

# def refind_N2_b(wavenr):
#    import numpy as np
#    return np.zeros(len(wavenr)) + 1.000298e0
# 
# def King_N2_b(wavenr):
#    import numpy as np
#    return np.ones(len(wavenr))

#from Sneep et al, 2005
def refind_CO2(wavenr):
   return 1e0 + (5799.25 / ((128908.9)**2e0 - wavenr**2e0) + 120.05 / ((89223.8)**2e0 - wavenr**2e0) + 5.3334 / ((75037.5)**2e0 - wavenr**2e0) + 4.3244 / ((67837.7)**2e0 - wavenr**2e0) + 0.1218145e0 / ((2418.136)**2e0 - wavenr**2e0)) * 1.1427e3
   
def King_CO2(wavenr):
   return 1.1364 + 25.3e-12*(wavenr)**2e0
   
#from Sneep et al, 2005
def refind_NH3(wavenr):
   return 1e0 + (46662e0 + 4.02e-6*(wavenr)**2e0) * 1e-8

def King_NH3(wavenr):
   import numpy as np
   return np.ones(len(wavenr))
   
# from Sneep et al, 2005
def refind_O2(wavenr):
   return 1e0 + (20564.8e0 + 24.809 / (4.09e9 - (wavenr)**2e0))*1e-8

def King_O2(wavenr):
   return 1.096e0 + 1.385e-11*(wavenr)**2e0 + 1.448e-20*(wavenr)**2e0
   
#from http://www.kayelaby.npl.co.uk/general_physics/2_5/2_5_7.html
def refind_H2O(wavenr):
   import numpy as np
   return np.zeros(len(wavenr)) + 1.000256e0

def King_H2O(wavenr):
   import numpy as np
   return np.ones(len(wavenr))

#http://refractiveindex.info/?group=GASES&material=Methane
def refind_CH4(wavenr):
   import numpy as np
   return np.zeros(len(wavenr)) + 1.000444e0

def King_CH4(wavenr):
   import numpy as np
   return np.ones(len(wavenr))

# I did not find data
def refind_O3(wavenr):
   import numpy as np
   return np.zeros(len(wavenr)) + 1e0

def King_O3(wavenr):
   import numpy as np
   return np.ones(len(wavenr))

#http://refractiveindex.info/?group=GASES&material=Hydrogen
# fitting Cauchy equation to tabulated data provided on the webpage
#http://en.wikipedia.org/wiki/Cauchy's_equation
def refind_H2(wavenr):
   import numpy as np
   wavel = 1e0/wavenr * 1e4
   return np.zeros(len(wavenr)) + 1.00013526e+00 + 1.24454843e-06/wavel**2e0

def King_H2(wavenr):
   import numpy as np
   return np.ones(len(wavenr))

#http://refractiveindex.info/?group=GASES&material=Helium
def refind_He(wavenr):
   import numpy as np
   # convert wavenumber to wavelength
   wavel = 1e0/wavenr * 1e4
   return np.zeros(len(wavel)) + 1e0 + 0.01470091 / (423.98 - wavel**(-2e0))

def King_He(wavenr):
   import numpy as np
   return np.ones(len(wavenr))

# saturation vapor pressures

# data from http://www.ddbst.com/en/online/Online_Calc_vap_Form.php
# and http://cires.colorado.edu/~voemel/vp.html
def Psat_H2O(T):
   import numpy as np
   #below 0 C. the equation is valid down to -100 C, the accuracy of the values below that T is uncertain
   if (T <=273.15):
      return 10e0**(-9.09718 * (273.16/T - 1e0) - 3.56654 * np.log10(273.16/T) + \
         0.876793 * (1e0 - T/273.16) + np.log10(6.1071) ) * 100e0
   #between 0 C and 100 C
   if (T > 273.15) and (T < 373.15):
      return 10e0**(8.07131 - 1730.63 / ( 233.426 + T - 273.15)) * 101.325/760e0*1e3
   #between 100 C and the triple point, 374 C 
   if (T >= 373.15) and (T < 647.15):
      return 10e0**(8.14019 - 1810.94	/ (244.485+ T - 273.15)) * 101.325/760e0*1e3
   else:
      return 0e0

# source: http://jcp.aip.org/resource/1/jcpsa6/v5/i1/p45_s1
# and http://www.ddbst.com/en/online/Online_Calc_vap_Form.php
# http://encyclopedia.airliquide.com/images_encyclopedie/VaporPressureGraph/Carbon_dioxide_Vapor_Pressure.GIF
def Psat_CO2(T):
   import numpy as np
   if (T < 210e0):
      return 10e0**(-(1354.210/T) + 8.69903 + 0.0015880 * T - 4.5107e-6 * T**2e0) * 1333.2239
   if (T >= 210e0) and (T < 304.15):
      return 10e0**(7.5322 - 835.06 / (268.223 + T - 273.15)) * 101.325/760e0*1e3
   else:
      return 0e0


# functions to calculate saturation vapor pressure of water from Pierrehumbert
#Saturation vapor pressure over ice (Smithsonian formula)
#    Input: Kelvin. Output: Pascal
def satvpi(T):
   import numpy as np
   #
   #  Compute es over ice (valid between -153 c and 0 c)
   #  see smithsonian meteorological tables page 350
   #
   #  Original source: GFDL climate model, circa 1995
   esbasi =    6107.1
   tbasi =     273.16
   #
   aa  = -9.09718 *(tbasi/T-1.0)
   b   = -3.56654 *np.log10(tbasi/T)
   c   =  0.876793*(1.0-T/tbasi)
   e   = np.log10(esbasi)
   esice = 10.**(aa+b+c+e)
   return .1*esice  #Convert to Pascals

#Saturation vapor pressure over liquid water (Smithsonian formula)
#    Input: Kelvin. Output: Pascal
def satvpw(T):
   import numpy as np
   #  compute es over liquid water between -20c and freezing.
   #  see smithsonian meteorological tables page 350.
   #
   #  Original source: GFDL climate model, circa 1995
   esbasw = 1013246.0
   tbasw =     373.16
   #
   aa  = -7.90298*(tbasw/T-1)
   b   =  5.02808*np.log10(tbasw/T)
   c   = -1.3816e-07*(  10.**( ((1-T/tbasw)*11.344)-1 )  )
   d   =  8.1328e-03*(  10.**( ((tbasw/T-1)*(-3.49149))-1)  )
   e   = np.log10(esbasw)
   esh2O  = 10.**(aa+b+c+d+e)
   return .1*esh2O  #Convert to Pascals

def satvpg(T):
   import numpy as np
#This is the saturation vapor pressure computation used in the
#GFDL climate model.  It blends over from water saturation to
#ice saturation as the temperature falls below 0C.
   
   T = np.array(T)
   satP = np.zeros(T.size)
   
   #print 'T:',T
   
   if T.size == 1:
      if ((T-273.16) <  -20.):
         satP = satvpi(T)
      if ( ((T-273.16) >= -20.)&((T-273.16)<=0.)):
         satP = 0.05*(273.16-T)*satvpi(T) + 0.05*(T-253.16)*satvpw(T)
      if ((T-273.16)>0.):
         satP = satvpw(T)
   else:
      satP[(T-273.16) < -20e0] = satvpi(T[(T-273.16) < -20e0])
      
      indx = ((T-273.16) >= -20e0) & ((T-273.16)<=0.)
      #print 'indx:',indx
      satP[indx] = 0.05*(273.16-T[indx])*satvpi(T[indx]) + 0.05*(T[indx]-253.16)*satvpw(T[indx])
      
      satP[(T-273.16)>0.] = satvpw(T[(T-273.16)>0.])
   
   #print 'satP:',satP
   
   return satP
