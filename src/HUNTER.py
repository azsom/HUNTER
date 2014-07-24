# this code performs an MC sampling to estimate what fraction of rocky planets could 
# harbor liquid water in the data set of Dressing&Charbonneau 2013.
# it uses precalcualted all-convective T-P profiles from Clima_expsum to speed up the calculation
# currently there are tables for H2, N2, and CO2 dominated atmospheres

import numpy as np
import data as d
import MC_functions as MC
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec
import pickle
from time import time


print '** start MC **'

# read in masses and radii of exoplanets
data = open('../data/catalogue.txt','r')
lines = [line.strip() for line in data]
data.close()
R_p = []
M_p = []
for line in lines:
   if line[0] != '#':
      line_split = line.split(',')
      if (line_split[3] != '') and (line_split[2] != ''):          
         R_p = np.append(R_p, np.float(line_split[3]))
         M_p = np.append(M_p, np.float(line_split[2]))

R_p = R_p * 10.97e0
M_p = M_p * 318e0

## data and functions needed in the code:

# Dressing&Charbonneau: calculating gaussian KDE kernel 
# this kernel is used to sample from the planet radius and flux distribution.
Rp_F_kernel, sum_occur = MC.occurrence_rate()
# Dressing & Charbonneau planet radius histrogram values
RP_data = np.array([0.5,0.7,1.0,1.4,2.0,2.8,4.0,5.7,8.0,11.3,16.0,22.6])

## integrating

print '* integrating *'

# number of MC shoots:
n = 1e5

# distributions to use
planet_mass = ['Mplanet_gauss_opt_sigma6']#['Mplanet_gauss_opt_sigma6','Mplanet_gauss_norm_sigma6','Mplanet_gauss_pes_sigma6','Mplanet_gauss_opt_sigma3','Mplanet_gauss_norm_sigma3','Mplanet_gauss_pes_sigma3']
atmosphere = ['atm_type_H2']#['atm_type_all','atm_type_N2_CO2','atm_type_H2','atm_type_N2','atm_type_CO2']
surf_pressure = ['Psurf_lognorm_1e6']#['Psurf_log_min1e0','Psurf_log_min1e2','Psurf_log_min1e4','Psurf_lognorm_1e4','Psurf_lognorm_1e5','Psurf_lognorm_1e6','Psurf_lognorm_1e7']
surf_albedo = ['alb_surf_uniform']#['alb_surf_uniform','alb_surf_normal','alb_surf_lognormal']
rel_humidity = ['relhum_lognorm']#['relhum_log','relhum_lin','relhum_lognorm']
N2_CO2_mix_rat = ['atm_N2_CO2_log']#['atm_N2_CO2_log','atm_N2_CO2_lin','atm_N2_CO2_lognorm']

habi = []

i = 0
# looping over all possible combinations
for iter1 in range(len(planet_mass)):
   for iter2 in range(len(atmosphere)):
      for iter3 in range(len(surf_pressure)):
         for iter4 in range(len(surf_albedo)):
            for iter5 in range(len(rel_humidity)):
               for iter6 in range(len(N2_CO2_mix_rat)):
                  
                  print '* ',planet_mass[iter1],',',atmosphere[iter2],',',surf_pressure[iter3],',',surf_albedo[iter4],',',rel_humidity[iter5],',',N2_CO2_mix_rat[iter6],' *'
                  
                  # distributions to call:
                  Mplanet = getattr(MC,planet_mass[iter1])
                  atm_type = getattr(MC,atmosphere[iter2])
                  Psurf = getattr(MC,surf_pressure[iter3])
                  alb_surf = getattr(MC,surf_albedo[iter4])
                  relhum = getattr(MC,rel_humidity[iter5])
                  atm_N2_CO2 = getattr(MC,N2_CO2_mix_rat[iter6])
                  # integrating
                  results = MC.MC_sample(n,Rp_F_kernel, sum_occur,Mplanet,atm_type,Psurf,alb_surf,relhum,atm_N2_CO2)
                  
                  # save results                  
                  [all_MS,all_flux,all_Rp,all_Mp,all_atmtype,all_Psurf,all_Tsurf,all_alb,all_relhum,all_CO2,all_N2,all_habi,planet,rocky,habitable] = results
                  
                  file = open('../results/30057.dat','w')
                  #file = open('../results/'+"%05d" % i + '.dat','w')
                  string = planet_mass[iter1]+','+atmosphere[iter2]+','+surf_pressure[iter3]+','+surf_albedo[iter4]+','+rel_humidity[iter5]+','+N2_CO2_mix_rat[iter6]+'\n'
                  file.write(string)
                  file.write(str(n)+'\n')
                  file.write(str(sum_occur)+'\n')
                  file.write(str(planet)+'\n')
                  file.write(str(rocky)+'\n')
                  file.write(str(habitable)+'\n')
                  for ii in range(len(all_MS)):
                     string = str(all_MS[ii])+' '+str(all_flux[ii])+' '+str(all_Rp[ii])+' '+str(all_Mp[ii])+' '+all_atmtype[ii]+' '+str(all_Psurf[ii])+' '+str(all_Tsurf[ii])+' '+str(all_alb[ii])+' '+str(all_relhum[ii])+' '+str(all_CO2[ii])+' '+str(all_N2[ii])+' '+all_habi[ii]+'\n'
                     file.write(string)
                  file.close()
                  
                  habi = np.append(habi,habitable/n * sum_occur)
                  
                  i = i + 1

print '* finished *'
