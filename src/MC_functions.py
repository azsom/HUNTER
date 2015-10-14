## functions called in MC_int

import numpy as np
import pickle
from scipy import stats
import random as ran
import scipy.interpolate 
from numerics import interp5
from numerics import interp6
import data as d
import matplotlib.pyplot as plt 
import matplotlib
from scipy.interpolate import interp1d

# stellar mass
#@profile
def stellar_mass():
   # IMF of Kroupa (2001)
   # IMF(Ms < 0.5) ~ m^(-1.3)
   # IMF(Ms > 0.5) ~ m^(-2.3)   
   
   # integrals of the two joint power laws. these values are used to decide which power law we sample from.
   int_l = 1e0/(-1.3+1.0)*2e0*0.5**(-1.3+1.0) - 1e0/(-1.3+1.0)*2e0*0.09**(-1.3+1.0)
   int_h = 1e0/(-2.3+1.0)*1.1**(-2.3+1.0) - 1e0/(-2.3+1.0)*0.5**(-2.3+1.0)
   
   int_low = int_l / (int_l+int_h)
   int_high = int_h / (int_l+int_h)
   
   mass = 0e0
   
   if ran.random() <= int_low:
      # we sample low mass part of the distribution
      while mass < 0.09 or mass > 0.5:
         mass = ((0.5**(-1.3+1.0) - 0.09**(-1.3+1.0))*ran.random() + 0.09**(-1.3+1.0))**(1.0/(-1.3+1.0))
   else:
      # we sample from the high mass end
      while mass < 0.5 or mass > 1.1:
         mass = ((1.1**(-2.3+1.0) - 0.5**(-2.3+1.0))*ran.random() + 0.5**(-2.3+1.0))**(1.0/(-2.3+1.0))
   return mass
   
# planet mass
#@profile
def Mplanet_gauss_opt_sigma3(RP,M_water,M_silicate,M_iron,R_silicate,R_water):
   # own recipe
   
   RP_crit = 1.79
   
   MP = 0e0
   if RP < RP_crit:
      # small rocky planets. the mass is normally distributed around the mean
      # and the sigma is the fraction of the mean
      M_w = M_water(RP)
      M_s = M_silicate(RP)
      M_i = M_iron(RP) 
      while MP > M_i or MP < M_w:
         mean = (M_w+M_s)/2e0
         sigma = mean / 3e0
         MP = np.random.normal(mean,sigma)
   if RP >= RP_crit:
      # a potentially gaseous planet:
      M_w = M_water(RP_crit)
      while MP < M_w:
         # mean from Weiss&Marcy 2014, Eq. 3
         mean = 2.69 * RP**0.93
         sigma = mean / 3e0
         MP = np.random.normal(mean,sigma)
   
   # critical density and mass above which this planet is considered rocky:
   # if it is pure water and below 50 earth masses
   dens_crit = 1e5
   if MP < 50 and RP < R_water(1000e0):
      M_crit = M_water(RP)
      dens_crit = M_crit*d.M_Earth / (4e0/3e0*np.pi*(RP*d.R_Earth)**3e0)

   return MP, dens_crit
   
def Mplanet_gauss_norm_sigma3(RP,M_water,M_silicate,M_iron,R_silicate,R_water):
   # own recipe
   
   RP_crit = 1.54
   
   MP = 0e0
   if RP < RP_crit:
      # small rocky planets. the mass is normally distributed around the mean
      # and the sigma is the fraction of the mean
      M_w = M_water(RP)
      M_s = M_silicate(RP)
      M_i = M_iron(RP) 
      while MP > M_i or MP < (M_w+M_s)/2e0:
         mean = M_s
         sigma = mean / 3e0
         MP = np.random.normal(mean,sigma)
   if RP >= RP_crit:
      # a potentially gaseous planet:
      M_w = M_water(RP_crit)
      while MP < M_w:
         # mean from Weiss&Marcy 2014, Eq. 3
         mean = 2.69 * RP**0.93
         sigma = mean / 3e0
         MP = np.random.normal(mean,sigma)
   
   # critical mass above which this planet is considered rocky:
   # if it is water-silicate mix and below 50 earth masses
   dens_crit = 1e5
   if MP < 50 and RP < R_silicate(1000e0):
      M_w = M_water(RP)
      M_s = M_silicate(RP)
      M_crit = (M_w + M_s) / 2e0      
      dens_crit = M_crit*d.M_Earth / (4e0/3e0*np.pi*(RP*d.R_Earth)**3e0)

   return MP, dens_crit   

def Mplanet_gauss_pes_sigma3(RP,M_water,M_silicate,M_iron,R_silicate,R_water):
   # own recipe
   
   RP_crit = 1.22
   
   MP = 0e0
   if RP < RP_crit:
      # small rocky planets. the mass is normally distributed around the mean
      # and the sigma is the fraction of the mean
      M_w = M_water(RP)
      M_s = M_silicate(RP)
      M_i = M_iron(RP) 
      while MP > M_i or MP < M_s:
         mean = (M_s+M_i)/2e0
         sigma = mean / 3e0
         MP = np.random.normal(mean,sigma)
   if RP >= RP_crit:
      # a potentially gaseous planet:
      M_w = M_water(RP_crit)
      while MP < M_w:
         # mean from Weiss&Marcy 2014, Eq. 3
         mean = 2.69 * RP**0.93
         sigma = mean / 3e0
         MP = np.random.normal(mean,sigma)
   
   # critical mass above which this planet is considered rocky:
   # if it is pure silicate and below 50 earth masses
   dens_crit = 1e5
   if MP < 50 and RP < R_silicate(1000e0):
      M_crit = M_silicate(RP)
      dens_crit = M_crit*d.M_Earth / (4e0/3e0*np.pi*(RP*d.R_Earth)**3e0)
   
   return MP, dens_crit


#@profile
def Mplanet_gauss_opt_sigma6(RP,M_water,M_silicate,M_iron,R_silicate,R_water):
   # own recipe
   
   RP_crit = 1.79
   
   MP = 0e0
   if RP < RP_crit:
      # small rocky planets. the mass is normally distributed around the mean
      # and the sigma is the fraction of the mean
      
      M_w = M_water(RP)
      M_s = M_silicate(RP)
      M_i = M_iron(RP) 
      while MP > M_i or MP < M_w:
         mean = (M_w+M_s)/2e0
         sigma = mean / 6e0
         MP = np.random.normal(mean,sigma)
   if RP >= RP_crit:
      # a potentially gaseous planet:
      M_w = M_water(RP_crit)
      while MP < M_w:
         # mean from Weiss&Marcy 2014, Eq. 3
         mean = 2.69 * RP**0.93
         sigma = mean / 6e0
         MP = np.random.normal(mean,sigma)
   
   # critical density and mass above which this planet is considered rocky:
   # if it is pure water and below 50 earth masses
   dens_crit = 1e5
   if MP < 50 and RP < R_water(1000e0):
      M_crit = M_water(RP)
      dens_crit = M_crit*d.M_Earth / (4e0/3e0*np.pi*(RP*d.R_Earth)**3e0)

   return MP, dens_crit
   
def Mplanet_gauss_norm_sigma6(RP,M_water,M_silicate,M_iron,R_silicate,R_water):
   # own recipe
   
   RP_crit = 1.54
   
   MP = 0e0
   if RP < RP_crit:
      # small rocky planets. the mass is normally distributed around the mean
      # and the sigma is the fraction of the mean
      M_w = M_water(RP)
      M_s = M_silicate(RP)
      M_i = M_iron(RP) 
      while MP > M_i or MP < (M_w+M_s)/2e0:
         mean = M_s
         sigma = mean / 6e0
         MP = np.random.normal(mean,sigma)
   if RP >= RP_crit:
      # a potentially gaseous planet:
      M_w = M_water(RP_crit)
      while MP < M_w:
         # mean from Weiss&Marcy 2014, Eq. 3
         mean = 2.69 * RP**0.93
         sigma = mean / 6e0
         MP = np.random.normal(mean,sigma)
   
   # critical mass above which this planet is considered rocky:
   # if it is water-silicate mix and below 50 earth masses
   dens_crit = 1e5
   if MP < 50 and RP < R_silicate(1000e0):
      M_w = M_water(RP)
      M_s = M_silicate(RP)
      M_crit = (M_w + M_s) / 2e0      
      dens_crit = M_crit*d.M_Earth / (4e0/3e0*np.pi*(RP*d.R_Earth)**3e0)

   return MP, dens_crit   

def Mplanet_gauss_pes_sigma6(RP,M_water,M_silicate,M_iron,R_silicate,R_water):
   # own recipe
   
   RP_crit = 1.22
   
   MP = 0e0
   if RP < RP_crit:
      # small rocky planets. the mass is normally distributed around the mean
      # and the sigma is the fraction of the mean
      M_w = M_water(RP)
      M_s = M_silicate(RP)
      M_i = M_iron(RP) 
      while MP > M_i or MP < M_s:
         mean = (M_s+M_i)/2e0
         sigma = mean / 6e0
         MP = np.random.normal(mean,sigma)
   if RP >= RP_crit:
      # a potentially gaseous planet:
      M_w = M_water(RP_crit)
      while MP < M_w:
         # mean from Weiss&Marcy 2014, Eq. 3
         mean = 2.69 * RP**0.93
         sigma = mean / 6e0
         MP = np.random.normal(mean,sigma)
   
   # critical mass above which this planet is considered rocky:
   # if it is pure silicate and below 50 earth masses
   dens_crit = 1e5
   if MP < 50 and RP < R_silicate(1000e0):
      M_crit = M_silicate(RP)
      dens_crit = M_crit*d.M_Earth / (4e0/3e0*np.pi*(RP*d.R_Earth)**3e0)
   
   return MP, dens_crit



# atmosphere type
#@profile
def atm_type_all():
   x = ran.random()
   if x <= 1e0/3e0:
     atm = 'H2'
   elif x > 1e0/3e0 and x <= 2e0/3e0:
     atm = 'CO2'
   else:
     atm = 'N2'
   return atm
   
def atm_type_N2_CO2():
   x = ran.random()
   if x <= 0.5:
      atm = 'CO2'
   else:
      atm = 'N2'
   return atm      
   
def atm_type_H2():
   atm = 'H2'
   return atm

def atm_type_N2():
   atm = 'N2'
   return atm

def atm_type_CO2():
   atm = 'CO2'
   return atm

# surface pressure
#@profile
def Psurf_log_min1e4(Mp,g_surf,Rp,atm):
   # uniformly distributed in log space
   Pmin = 1e4
   # max surface pressure
   # F = m*g = P*A 
   # mg = P*4*pi*r^2 
   # P = (mg) / (4*pi*r^2)
   if atm != 'H2':
      Pmax = 1e-2 * (Mp*d.M_Earth) * g_surf / (4e0*np.pi*(Rp*d.R_Earth)**2e0)
   else:
      Pmax = 1e-5 * (Mp*d.M_Earth) * g_surf / (4e0*np.pi*(Rp*d.R_Earth)**2e0)      
   a = np.log(Pmin)
   b = np.log(Pmax)
   PS = np.exp(a + (b-a) * ran.random())
   return PS

def Psurf_log_min1e2(Mp,g_surf,Rp,atm):
   # uniformly distributed in log space
   Pmin = 1e2
   # max surface pressure
   # F = m*g = P*A 
   # mg = P*4*pi*r^2 
   # P = (mg) / (4*pi*r^2)
   if atm != 'H2':
      Pmax = 1e-2 * (Mp*d.M_Earth) * g_surf / (4e0*np.pi*(Rp*d.R_Earth)**2e0)
   else:
      Pmax = 1e-5 * (Mp*d.M_Earth) * g_surf / (4e0*np.pi*(Rp*d.R_Earth)**2e0)      
   a = np.log(Pmin)
   b = np.log(Pmax)
   PS = np.exp(a + (b-a) * ran.random())
   return PS

def Psurf_log_min1e0(Mp,g_surf,Rp,atm):
   # uniformly distributed in log space
   Pmin = 1e0
   # max surface pressure
   # F = m*g = P*A 
   # mg = P*4*pi*r^2 
   # P = (mg) / (4*pi*r^2)
   if atm != 'H2':
      Pmax = 1e-2 * (Mp*d.M_Earth) * g_surf / (4e0*np.pi*(Rp*d.R_Earth)**2e0)
   else:
      Pmax = 1e-5 * (Mp*d.M_Earth) * g_surf / (4e0*np.pi*(Rp*d.R_Earth)**2e0)      
   a = np.log(Pmin)
   b = np.log(Pmax)
   PS = np.exp(a + (b-a) * ran.random())
   return PS

def Psurf_lognorm_1e4(Mp,g_surf,Rp,atm):
   Pmin = 0e0
   if atm != 'H2':
      Pmax = 1e-2 * (Mp*d.M_Earth) * g_surf / (4e0*np.pi*(Rp*d.R_Earth)**2e0)
   else:
      Pmax = 1e-5 * (Mp*d.M_Earth) * g_surf / (4e0*np.pi*(Rp*d.R_Earth)**2e0) 
   PS = -1e0
   while PS < Pmin or PS > Pmax:
      PS = np.exp(np.random.normal(0e0,1e0))
      mode = np.exp(0e0-1e0**2e0)
      PS = PS / mode * 1e4
   return PS

def Psurf_lognorm_1e5(Mp,g_surf,Rp,atm):
   Pmin = 0e0
   if atm != 'H2':
      Pmax = 1e-2 * (Mp*d.M_Earth) * g_surf / (4e0*np.pi*(Rp*d.R_Earth)**2e0)
   else:
      Pmax = 1e-5 * (Mp*d.M_Earth) * g_surf / (4e0*np.pi*(Rp*d.R_Earth)**2e0) 
   PS = -1e0
   while PS < Pmin or PS > Pmax:
      PS = np.exp(np.random.normal(0e0,1e0))
      mode = np.exp(0e0-1e0**2e0)
      PS = PS / mode * 1e5
   return PS

def Psurf_lognorm_1e6(Mp,g_surf,Rp,atm):
   Pmin = 0e0
   if atm != 'H2':
      Pmax = 1e-2 * (Mp*d.M_Earth) * g_surf / (4e0*np.pi*(Rp*d.R_Earth)**2e0)
   else:
      Pmax = 1e-5 * (Mp*d.M_Earth) * g_surf / (4e0*np.pi*(Rp*d.R_Earth)**2e0) 
   PS = -1e0
   while PS < Pmin or PS > Pmax:
      PS = np.exp(np.random.normal(0e0,1e0))
      mode = np.exp(0e0-1e0**2e0)
      PS = PS / mode * 1e6
   return PS

def Psurf_lognorm_1e7(Mp,g_surf,Rp,atm):
   Pmin = 0e0
   if atm != 'H2':
      Pmax = 1e-2 * (Mp*d.M_Earth) * g_surf / (4e0*np.pi*(Rp*d.R_Earth)**2e0)
   else:
      Pmax = 1e-4 * (Mp*d.M_Earth) * g_surf / (4e0*np.pi*(Rp*d.R_Earth)**2e0) 
   PS = -1e0
   while PS < Pmin or PS > Pmax:
      #print Pmax/1e5, Mp,Rp,g_surf,PS/1e5
      PS = np.exp(np.random.normal(0e0,1e0))
      mode = np.exp(0e0-1e0**2e0)
      PS = PS / mode * 1e7
   return PS

# surface albedo
albmin = 0e0
albmax = 1e0

#@profile
def alb_surf_uniform():
   # uniformly distributed
   alb = ran.uniform(albmin,albmax)
   return alb
def alb_surf_normal():
   alb = -1e0
   while alb < albmin or alb > albmax:
      alb = np.random.normal(0.2,0.1)
   if alb < albmin or alb > albmax:
      print 'albedo aout of range.',alb
      raise ValueError
      
   return alb
def alb_surf_lognormal():
   alb = -1e0
   while alb < albmin or alb > albmax:
      alb = np.random.lognormal(0.0,1e0)
      mode = np.exp(0e0-1e0**2e0)
      alb = alb / mode * 0.1
   if alb < albmin or alb > albmax:
      print 'albedo out of range.'
      raise ValueError
   return alb


# relative humidity   
#@profile
def relhum_log():
   # uniformly distributed in log space
   relhummin = 1e-2
   relhummax = 1e0
   a = np.log(relhummin)
   b = np.log(relhummax)
   relh = np.exp(a + (b-a) * ran.random())
   return relh
   
def relhum_lin():
   relhummin = 1e-2
   relhummax = 1e0
   relh = ran.uniform(relhummin,relhummax)
   return relh
   
def relhum_lognorm():
   # the mode is 0.5 (average tropospheric relhum needed to reproduce surface temperature on Earth)
   relhummin = 1e-2
   relhummax = 1e0   
   relh = 0e0
   while relh < relhummin or relh > relhummax:
      relh = np.exp(np.random.normal(0e0,1e0))
      mode = np.exp(0e0-1e0**2e0)
      relh = relh / mode * 0.5
   return relh
   
# CO2 mixing ratio
#@profile
def atm_N2_CO2_log(atm):
   # uniformly distributed in log space between 1e-5 and 0.5
   if atm == 'N2':
      a = np.log(1e-5)
      b = np.log(0.5)
      CO2 = np.exp(a + (b-a) * ran.random())
      N2 = 1e0 - CO2
   if atm == 'CO2':
      a = np.log(1e-5)
      b = np.log(0.5)
      N2 = np.exp(a + (b-a) * ran.random())
      CO2 = 1e0 - N2
   
   return N2,CO2

def atm_N2_CO2_lin(atm):
   # uniformly distributed between 0 and 0.5
   if atm == 'N2':
      CO2 = ran.uniform(0e0,0.5)
      N2 = 1e0 - CO2
   if atm == 'CO2':
      N2 = ran.uniform(0e0,0.5)
      CO2 = 1e0 - N2
   return N2,CO2

def atm_N2_CO2_lognorm(atm):
   # lognorm distribution, where CO2 mode is 4e-2 as on Earth, and
   # N2 mode is 3.5e-2 as on Venus.
   if atm == 'N2':
      CO2min = 1e-5
      CO2max = 0.5
      CO2 = 0e0
      while CO2 < CO2min or CO2 > CO2max:
         CO2 = np.exp(np.random.normal(0e0,1e0))
         mode = np.exp(0e0-1e0**2e0)
         CO2 = CO2 / mode * 4e-2
      N2 = 1e0 - CO2
   if atm == 'CO2':
      N2min = 1e-5
      N2max = 0.5
      N2 = 0e0
      while N2 < N2min or N2 > N2max:
         N2 = np.exp(np.random.normal(0e0,1e0))
         mode = np.exp(0e0-1e0**2e0)
         N2 = N2 / mode * 3.5e-2
      CO2 = 1e0 - N2
   return N2,CO2

# this function calculates new occurrence rates on a non-uniform grid from the Dressing&Charbonneau data
def occurrence_rate():

   print '* calculating gaussian KDE of the Dressing&Charbonneau data *'

   # to fit CDPP
   def func(x, a, b):
      return np.sqrt(a/(b*x))

   # constants in SI units
   AU = 1.496e11
   RSun = 6.955e8
   MSun = 1.9891e30
   Grav = 5.67e-8
   REarth = 6.371e6
   day = 86400e0
   solconst = 1361e0

   ## read in data

   # read in from table 2
   # kepler ID, orbital period, rp/rs

   f = open('../data/apj464844t2_mrt.txt','r')
   lines = f.readlines()[45:]
   f.close()

   kepIDp = []
   orbp = []
   RpRs = []
   flux_test = []
   radius_test = []
   asemi_test = []

   for line in lines:
      line_split = line.split(' ')
      kepIDp = np.append(kepIDp,line_split[-15])
      orbp = np.append(orbp,float(line_split[-13]))
      RpRs = np.append(RpRs,float(line_split[-11]))
      flux_test = np.append(flux_test,float(line_split[-6]))
      radius_test = np.append(radius_test,float(line_split[-9]))
      asemi_test = np.append(asemi_test,float(line_split[-2])*float(line_split[-12])*RSun/AU)

   # read in from table 1
   # kepler ID, effective temperature and error bar, stellar radius and error bar, stellar mass and error bar

   f = open('../data/apj464844t1_mrt.txt','r')
   lines = f.readlines()[28:]
   f.close()

   kepID = []
   Teff = []
   Teffp = []
   Teffn = []
   Rs = []
   Rsp = []
   Rsn = []
   Ms = []
   Msp = []
   Msn = []

   for line in lines:
      line_split = line.split(' ')
      kepID = np.append(kepID,line_split[-20])
      Teff = np.append(Teff,float(line_split[-19]))
      Teffp = np.append(Teffp,float(line_split[-18]))
      Teffn = np.append(Teffn,float(line_split[-17]))
      Rs = np.append(Rs,float(line_split[-16]))
      Rsp = np.append(Rsp,float(line_split[-15]))
      Rsn = np.append(Rsn,float(line_split[-14]))
      Ms = np.append(Ms,float(line_split[-13]))
      Msp = np.append(Msp,float(line_split[-12]))
      Msn = np.append(Msn,float(line_split[-11]))


   # CDPP values and LIVETIMEs
   file = open('../data/CDPP.dat','r')
   [kepIDCDPP,LIVETIME,CDPP] = pickle.load(file)
   file.close()

   ## calculate incoming fluxes and planet radii with error bars

   # stellar luminosities and errors (most of the error is in Rs)
   Ls = 4e0 * np.pi * (Rs*RSun)**2e0 * Teff**4e0 * Grav # [W]
   Lsp = 4e0 * np.pi * ((Rs+Rsp)*RSun)**2e0 * (Teff)**4e0 * Grav # [W]
   Lsn = 4e0 * np.pi * ((Rs-Rsn)*RSun)**2e0 * (Teff)**4e0 * Grav # [W]

   # converting orbital period to semi major axis
   # indices of kepIDs with planets
   indx = [np.where(kepID == i)[0][0] for i in kepIDp]

   mu = 6.67e-11 * Ms[indx] * MSun # [m^3/s^2]
   mup = 6.67e-11 * (Ms[indx] + Msp[indx]) * MSun
   mun = 6.67e-11 * (Ms[indx] - Msn[indx]) * MSun

   a_semi = ((orbp*day/(2e0*np.pi))**2e0*mu)**(1e0/3e0) / AU # [AU]
   a_semip = ((orbp*day/(2e0*np.pi))**2e0*mup)**(1e0/3e0) / AU
   a_semin = ((orbp*day/(2e0*np.pi))**2e0*mun)**(1e0/3e0) / AU

   # calculating the flux received at semi manor axis with error bars (the error is mostly due to Ls)
   flux = Ls[indx] / (4e0 * np.pi * ((a_semi * AU))**2e0) / solconst # [solar constant]
   fluxp = Lsp[indx] / (4e0 * np.pi * ((a_semi * AU))**2e0) / solconst
   fluxn = Lsn[indx] / (4e0 * np.pi * ((a_semi * AU))**2e0) / solconst

   # planet radii and errors
   Rp = RpRs * Rs[indx] * RSun / REarth # [Earth radius]
   Rpp = RpRs * Rsp[indx] * RSun / REarth
   Rpn = RpRs * Rsn[indx] * RSun / REarth

   ## calculate occurrence rate for each point
   # calculate Ns: around how many stars in the sample we can observe the planet with period p and radius Rp
   # calculate for all stars SNR = (Rp/Rs)^2 * sqrt(n) / CDPP where n is the number of transits within 
   #   the observation time (LIVETIME), CDPP is the sigma CDPP interpolated to the transit duration.
   # Ns is the number of stars where SNR > 7.1
   #
   # the occurrence rate around each planet is 
   # f(Rp,F) = a / Rs / Ns
   #

   # loop through each planet and each star.
   Ns = np.zeros(len(kepIDp))
   p_transit = np.zeros(len(kepIDp))
   f_occur = np.zeros(len(kepIDp))

   flux_points = []
   Rp_points = []

   for i in range(len(kepIDp)):
      count = 0
      for ii in range(len(kepID)):
         # calculate SNR
         # the orbital period of planet[i] receiving flux[i] around star[ii]
         asem = np.sqrt(Ls[ii] / (flux[i]*solconst*4e0*np.pi)) # [m]
         orb = 2e0 * np.pi * np.sqrt(asem**3e0 / (6.67e-11 * Ms[ii] * MSun)) # [s]
         
         # number of transits observed in Q1-6
         n = LIVETIME[ii]/(orb/day)
         
         # transit time
         b = 0e0
         t_trans = orb/np.pi * np.arcsin((Rs[ii]*RSun)/asem*np.sqrt((1e0 + (Rp[i]*REarth)/(Rs[ii]*RSun))**2e0 - b**2e0)) / 3600e0 # [hours]
         
         # estimate CDPP based on Christiansen et al 2012 (above eq. 4)
         t_CDPP = np.array([3e0,6e0,12e0])
         # find t_CDPP nearest to t_trans
         ind = np.where(np.abs(t_CDPP - t_trans) == np.min(np.abs(t_CDPP - t_trans)))
         # estimate CDPP:
         CDPP_p = np.sqrt(t_CDPP[ind[0][0]]/t_trans) * CDPP[ii,ind[0][0],1]
         
         if CDPP_p > 0e0:
            SNR = ((Rp[i]*REarth)/(Rs[ii]*RSun))**2e0 * np.sqrt(n) / (CDPP_p/1e6)
         else:
            SNR = 0e0
         
         if SNR > 7.1:
            count = count+1
      
      Ns[i] = count  
      # probability of transt (1/p_transit is the number in upper left corner in brackets)
      p_transit[i] = (Rs[indx[i]]*RSun) / (a_semi[i]*AU)
      f_occur[i] = 1e0 / (Ns[i] * p_transit[i])
      
      flux_points = np.append(flux_points, np.zeros(int(round(f_occur[i]*10000e0)))+flux[i])
      Rp_points = np.append(Rp_points, np.zeros(int(round(f_occur[i]*10000e0)))+Rp[i])
      
   
   values = np.vstack([np.log10(flux_points), np.log10(Rp_points)])
   kernel = stats.gaussian_kde(values, bw_method = len(flux)**(-1./(2.0+4.0)))

# ## generate gaussianKDE.pdf:
#    sample = kernel.resample(size=10000)
# 
#    x_sample = 10e0**sample[0,:]
#    y_sample = 10e0**sample[1,:]
# 
# 
#    ## if x and y ranges are changes, check the tick locs and labels below!
#    xmin = -1e0
#    xmax = 2e0
#    ymin = -0.5
#    ymax = 1.5
#    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
#    positions = np.vstack([X.ravel(), Y.ravel()])
#    Z = np.reshape(kernel(positions).T, X.shape)
# 
#    matplotlib.rcParams.update({'font.size': 18})
# 
#    plt.figure(figsize=(8,6))
# 
#    plt.axis([xmin,xmax,ymin,ymax])
#    plt.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,extent=[xmin, xmax, ymin, ymax],aspect='auto')
#    plt.scatter(np.log10(flux),np.log10(Rp),s=80,marker='+',color='k',linewidth=2)
#    #CB=plt.colorbar(CS)
#    #CB.set_label('occurrence rate')
#    tick_locs = [-1,0,1,2]
#    tick_lbls = 10e0**np.array(tick_locs)
#    plt.xticks(tick_locs, tick_lbls)
#    tick_locs = [-0.5,0,0.5,1,1.5]
#    tick_lbls = [0.32,1e0,3.16,10,31.6]
#    plt.yticks(tick_locs,tick_lbls)
#    plt.tick_params(length=5,width=2,which='major')
#    plt.tick_params(length=4,width=2,which='minor')
#    plt.xlabel('incoming stellar flux [solar constant]')
#    plt.ylabel('planet radius [$R_\oplus$]')
#    plt.savefig('../figs/gaussianKDE.pdf', bbox_inches='tight')
#    plt.close()
#    
#   stop

   return kernel,np.sum(f_occur)
   

# monte carlo integration
#@profile
def MC_int(n,Rp_F_kernel, sum_occur,Mplanet,atm_type,Psurf,alb_surf,relhum,atm_N2_CO2):

   ## data and functions
   # stuff u need to calculate the boiling point
   T_grid = np.linspace(273,647,num = 647-273+1,endpoint=True)
   Psat = d.satvpg(T_grid)

   # stellar parameters from Allard:
   Mstar_B = [0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.05,1.10] 
   Teff_B = [2641,2812,3292,3436,3522,3649,3958,4418,4910,5370,5814,5987,6070]
   LogL_B = [-3.27,-3.07,-2.31,-1.96,-1.71,-1.46,-1.15,-0.84,-0.55,-0.27,0.01,0.15,0.27]

   # reading the atmosphere grids
   file = open('../grids/CO2atm_grid.out','r')
   [g_surf_grid,T_boil,T_num,Psurf_grid,relhum_grid,N2_grid,semi_major_CO2,fluxfac_CO2,T_top_CO2] = pickle.load(file)
   file.close()

   file = open('../grids/N2atm_grid.out','r')
   [g_surf_grid,T_boil,T_num,Psurf_grid,relhum_grid,CO2_grid,semi_major_N2,fluxfac_N2,T_top_N2] = pickle.load(file)
   file.close()

   file = open('../grids/H2atm_grid.out','r')
   [g_surf_grid,T_boil,T_num,Psurf_grid,relhum_grid,semi_major_H2,fluxfac_H2,T_top_H2] = pickle.load(file)
   file.close()

   alb_grid = np.linspace(albmin,albmax,num=5,endpoint=True)

   # functions to calculate M_iron, M_silicate, and M_water based on Seager et al 2007:
   pmass = 10e0**np.linspace(np.log10(0.01),np.log10(1000),num=1000)
   R_water = 10e0**(-0.209396 + np.log10(pmass/5.52)/3e0 - 0.0807*(pmass/5.52)**(0.375)) * 4.43
   R_sil = 10e0**(-0.209594 + np.log10(pmass/10.55)/3e0 - 0.0799*(pmass/10.55)**(0.413)) * 3.90
   R_iron = 10e0**(-0.209490 + np.log10(pmass/5.8)/3e0 - 0.0804*(pmass/5.8)**(0.394)) * 2.52   

   M_water = interp1d(R_water,pmass)
   M_silicate = interp1d(R_sil,pmass)
   M_iron = interp1d(R_iron,pmass)

   R_water = interp1d(pmass,R_water)
   R_silicate = interp1d(pmass,R_sil)
   R_iron = interp1d(pmass,R_iron)


   ## arrays to save the results
   planet = 0e0
   all_MS = np.zeros(n)
   all_flux = np.zeros(n)
   all_Rp = np.zeros(n)
   all_Mp = np.zeros(n)
   all_atmtype = np.chararray(n) + '    '
   all_Psurf = np.zeros(n)
   all_Tsurf = np.zeros(n)
   all_alb = np.zeros(n)
   all_relhum = np.zeros(n)
   all_CO2 = np.zeros(n)
   all_N2 = np.zeros(n)
   all_habi = np.chararray(n) + ' '

   # number of rocky planets drawn from the distribution
   rocky = 0e0

   # number of habitable planets drawn from the distribution
   habitable = 0e0

   ## start the integration
   i = 0
   while i < n:
      # if i % (n/10) == 0:
      #    print i,rocky,habitable

      # draw stellar properties
      # Teff should be between 3100 and 4000 K - Dressing&Charbonneau
      TS = 0e0
      while TS < 5500e0 or TS > 6000e0:
         MS = stellar_mass()   
         spl = scipy.interpolate.splrep(Mstar_B, Teff_B)
         TS = scipy.interpolate.splev(MS, spl, ext=1)
      
#       # draw planet radius and orbital period values from the distribution of Dressing&Charbonneau 2013
#       flux = 0e0
#       Rp = 0e0
#       while flux < 0e0 or Rp < 0.5:
#          sample = Rp_F_kernel.resample(size=1)
#          flux = 10e0**sample[0,0]
#          Rp = 10e0**sample[1,0]
      
      # planet flux and radius are uniformly distributed
      logF_min = -2e0
      logF_max = 2e0
      logRp_min = np.log10(0.5)
      logRp_max = np.log10(5e0)
      flux = 10e0**np.random.uniform(logF_min,logF_max)
      Rp = 10e0**np.random.uniform(logRp_min,logRp_max)
      
      planet = planet + 1e0
      Mp, dens_crit = Mplanet(Rp,M_water,M_silicate,M_iron,R_silicate,R_water)
      
      g_surf = d.grav_const * (Mp*d.M_Earth) / (Rp*d.R_Earth)**2e0
      
      dens = Mp*d.M_Earth / (4e0/3e0*np.pi*(Rp*d.R_Earth)**3e0)
      
      all_MS[i] = MS
      all_Rp[i] = Rp
      all_Mp[i] = Mp
      all_flux[i] = flux
      all_habi[i] = 'no'
      all_atmtype[i] = 'SuperE'
      
      # is the planet rocky? 
      # check whether its density is largen than a critical value
      if dens > dens_crit:
         rocky = rocky + 1e0
         
         # atmosphere type: N2, CO2 or H2 dominated
         atm = atm_type()
         # surface pressure, albedo, and relative humidity values
         PS = Psurf(Mp,g_surf,Rp,atm)
         
         if Rp < 0.5 or PS < 1e3:
            all_atmtype[i] = 'rock'
         else:
            # rocky planet with atmosphere
            all_atmtype[i] = 'p_w_atm'
            alb = alb_surf()
            RH = relhum()
            
            # if the atmosphere is not H2 dominated, we draw a CO2 or N2 mixing ratio
            if atm != 'H2':
               N2, CO2 = atm_N2_CO2(atm)
            else:
               CO2 = 0e0
               N2 = 0e0
            
            # stupidity check: are parameters within the grid?
            if g_surf < g_surf_grid[0] or g_surf > g_surf_grid[-1]:
               #print 'surface gravity outside range.'
               #print Mp,Rp,g_surf
               if g_surf > g_surf_grid[-1]:
                  g_surf = g_surf_grid[-1]
               else:
                  g_surf = g_surf_grid[0]
            if alb < alb_grid[0] or alb > alb_grid[-1]:
               print 'albedo is outside grid.'
               raise ValueError
            if RH < relhum_grid[0] or RH > relhum_grid[-1]:
               print 'relative humidity is outside grid.'
               raise ValueError
            
            all_atmtype[i] = atm
            all_Psurf[i] = PS
            all_alb[i] = alb
            all_relhum[i] = RH
            all_CO2[i] = CO2
            all_N2[i] = N2
            
            # given PS and RH, calculate the boiling point of water.  
            PS_b = Psat * RH + PS
            if np.any(Psat > PS_b):
               indx = np.min(np.where(Psat > PS_b))
               T_B = T_grid[indx]
            else:
               T_B = 647e0
            
            # T_gr where the flux array should be interpolated to
            T_gr = np.linspace(273e0,T_B,num = 10,endpoint=True)         
            # normalized T_surf grid, 0e0 = 273 K, 1e0 = boiling temperature
            T_surf_norm = np.linspace(0e0,1e0,num=T_num,endpoint = True)
            
            T_surf = 0e0
            T_top = 0e0
            
            # atmosphere is set up. let's check whether the surface climate is habitable :)
            if atm == 'H2' and PS < 1e7:     
               if flux >= np.min(fluxfac_H2) and flux <= np.max(fluxfac_H2):
                  
                  # flux values at T_gr
                  flux_grid = interp5(g_surf_grid,alb_grid, T_surf_norm,Psurf_grid,relhum_grid,fluxfac_H2,np.zeros(10)+g_surf, np.zeros(10)+alb, (T_gr - 273e0)/(T_B-273e0),np.zeros(10)+PS, np.zeros(10)+RH)
                  
                  # top of atmosphere temperature at T_gr
                  T_top_grid = interp5(g_surf_grid,alb_grid, T_surf_norm,Psurf_grid,relhum_grid,T_top_H2,np.zeros(10)+g_surf, np.zeros(10)+alb, (T_gr - 273e0)/(T_B-273e0),np.zeros(10)+PS, np.zeros(10)+RH)
                  
                  # check whether liuqid water can be on the surface at flux value
                  # also check whether T_top isn't too large!
                  T_surf = np.interp(flux,flux_grid,T_gr,left = 0e0, right = 1000e0)
                  T_top = np.interp(flux,flux_grid,T_top_grid,left = 0e0, right = 1000e0)
                  
            if atm == 'N2' and PS < 1e7:
               if flux >= np.min(fluxfac_N2) and flux <= np.max(fluxfac_N2):
                  
                  # flux values at T_gr
                  flux_grid = interp6(g_surf_grid,alb_grid, T_surf_norm,Psurf_grid,relhum_grid, CO2_grid,fluxfac_N2,np.zeros(10)+g_surf, np.zeros(10)+alb, (T_gr - 273e0)/(T_B-273e0),np.zeros(10)+PS, np.zeros(10)+RH, np.zeros(10)+CO2)
                  
                  # top of atmosphere temperature at T_gr
                  T_top_grid = interp6(g_surf_grid,alb_grid, T_surf_norm,Psurf_grid,relhum_grid, CO2_grid,T_top_N2,np.zeros(10)+g_surf, np.zeros(10)+alb, (T_gr - 273e0)/(T_B-273e0),np.zeros(10)+PS, np.zeros(10)+RH, np.zeros(10)+CO2)
                  
                  # check whether liuqid water can be on the surface at flux value
                  # also check whether T_top isn't too large!
                  T_surf = np.interp(flux,flux_grid,T_gr,left = 0e0, right = 1000e0)
                  T_top = np.interp(flux,flux_grid,T_top_grid,left = 0e0, right = 1000e0)
                  
            if atm == 'CO2' and PS < 1e7:
               if flux >= np.min(fluxfac_CO2) and flux <= np.max(fluxfac_CO2):
                  
                  # flux values at T_gr
                  flux_grid = interp6(g_surf_grid,alb_grid, T_surf_norm,Psurf_grid,relhum_grid, N2_grid,fluxfac_CO2,np.zeros(10)+g_surf, np.zeros(10)+alb, (T_gr - 273e0)/(T_B-273e0),np.zeros(10)+PS, np.zeros(10)+RH, np.zeros(10)+N2)
                  
                  # top of atmosphere temperature at T_gr
                  T_top_grid = interp6(g_surf_grid,alb_grid, T_surf_norm,Psurf_grid,relhum_grid, N2_grid,T_top_CO2,np.zeros(10)+g_surf, np.zeros(10)+alb, (T_gr - 273e0)/(T_B-273e0),np.zeros(10)+PS, np.zeros(10)+RH, np.zeros(10)+N2)
                  
                  # check whether liuqid water can be on the surface at flux value
                  # also check whether T_top isn't too large!
                  T_surf = np.interp(flux,flux_grid,T_gr,left = 0e0, right = 1000e0)
                  T_top = np.interp(flux,flux_grid,T_top_grid,left = 0e0, right = 1000e0)
            
            # let's check whether the atmosphere is stable 
            # following Heng & Kopparla, 2012, I compare the advective and radiative time scales
            if T_surf >= 273e0 and T_surf <= 647e0 and T_top <= 100e0:
               # number of degrees of freedom of freedom of the gas, we adapt 5 for simplicity
               ndof = 5e0
               # adiabatic gas index
               gamma = 1e0 + 2e0/ndof
               # irradiation temperature
               Tirr = (flux/d.sigma)**(1e0/4e0)
               # mean molecular mass
               if atm == 'H2':
                  mean_mm = d.mu*2e0
                  cp = (2e0+ndof)*d.R_gas/(2e0*2e0)
               if atm == 'N2':
                  mean_mm = d.mu*(28e0*N2 + 44e0*CO2)
                  cp = (2e0+ndof)*d.R_gas/(2e0*(28e0*N2 + 44e0*CO2))
               if atm == 'CO2':
                  mean_mm = d.mu*(28e0*N2 + 44e0*CO2)
                  cp = (2e0+ndof)*d.R_gas/(2e0*(28e0*N2 + 44e0*CO2))
               # sound speed
               cs = np.sqrt(gamma*d.kb*T_surf/mean_mm)
               
               # advective time scale:
               t_adv = Rp*d.R_Earth / cs               
               # radiative time scale
               t_rad = cp * T_surf * PS / (g_surf*flux*1370e0)
               
               #print t_adv, t_rad
               
               if t_adv < t_rad:
                  habitable = habitable + 1e0
                  all_habi[i] = 'yes'
                  all_Tsurf[i] = T_surf
               else:
                  all_habi[i] = 'no'
                  all_Tsurf[i] = T_surf
            else:
               all_habi[i] = 'no'
               all_Tsurf[i] = T_surf
            
            
            
            
      i = i + 1

   print 'fraction of rocky planets: ', rocky/n * sum_occur
   print 'fraction of habitable planets: ', habitable/n * sum_occur

   results = [all_MS,all_flux,all_Rp,all_Mp,all_atmtype,all_Psurf,all_Tsurf,all_alb,all_relhum,all_CO2,all_N2,all_habi,planet,rocky,habitable]   
   
   return results
   
