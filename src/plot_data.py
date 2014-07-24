# to read and plot the output of HUNTE

import matplotlib.pyplot as plt 
import numpy as np
import data as d
import os
from scipy.interpolate import interp1d
from bisect import bisect_left
import matplotlib
import matplotlib.gridspec as gridspec
import pickle

#matplotlib.rcParams.update({'font.size': 18})

## plot HZ occurrence rate histogram

# read data files
i = 0
path = '../results'
for (path, dir, files) in os.walk(path):
   i += 1
   if i == 1:
      break
      
files = np.array(files)
mask = [(i[0] != '0') and (i[0] != '.') for i in files]

files = files[np.where(mask)]

Rp_lim = [0.5,5e0]
flux_lim = [0.1,10e0]
nr_bin = 18
Rp_grid = 10e0**np.linspace(np.log10(Rp_lim[0]),np.log10(Rp_lim[1]),num = nr_bin+1)
flux_grid = 10e0**np.linspace(np.log10(flux_lim[0]),np.log10(flux_lim[1]),num = nr_bin+1)


## use this if files are read for the first time:

# arrays to store stuff
Mp_dist = []
atm_dist = []
Psurf_dist = []
alb_dist = []
relhum_dist = []
N2_CO2_dist = []

HZ_occur = []
rock_occur = []
SE_occur = []
non_HZ_occur = []


# for plotting the Rp and flux distributions
planet_total = np.zeros([len(files),nr_bin,nr_bin])
habi_total = np.zeros([len(files),nr_bin,nr_bin])
Rp_dist = np.zeros([len(files),nr_bin])
Rp_habi = np.zeros([len(files),nr_bin])
flux_dist = np.zeros([len(files),nr_bin])
flux_habi = np.zeros([len(files),nr_bin])

jj = 0

for i in files:
   print i
   file = open('../results/'+i,'r')
   lines = file.readlines()
   
   line = lines[0]
   line_split = line.split(',')
   Mp_dist = np.append(Mp_dist,line_split[0])
   atm_dist = np.append(atm_dist,line_split[1])
   Psurf_dist = np.append(Psurf_dist,line_split[2])
   alb_dist = np.append(alb_dist,line_split[3])
   relhum_dist = np.append(relhum_dist,line_split[4])
   N2_CO2_dist = np.append(N2_CO2_dist,line_split[5][:-1])
   
   f_occur = float(lines[2])
   n_total = float(lines[3])
   n_rocky = float(lines[4])
   n_HZ = float(lines[5])
   
   HZ_occur = np.append(HZ_occur,n_HZ/n_total * f_occur)
   
   for j in lines[6:]:
      line_split = j.split(' ')
      Rp = float(line_split[2])
      flux = float(line_split[1])
      habi = line_split[11]
      indx_Rp = bisect_left(Rp_grid,Rp) - 1
      indx_flux = bisect_left(flux_grid,flux) - 1
      if indx_Rp >= 0 and indx_Rp <= nr_bin-1:
         if indx_flux >= 0 and indx_flux <= nr_bin-1:
            planet_total[jj,indx_Rp,indx_flux] = planet_total[jj,indx_Rp,indx_flux] + 1
            Rp_dist[jj,indx_Rp] = Rp_dist[jj,indx_Rp] + 1 
            flux_dist[jj,indx_flux] = flux_dist[jj,indx_flux] + 1 
            if habi[0:2] == 'ye':
               habi_total[jj,indx_Rp,indx_flux] = habi_total[jj,indx_Rp,indx_flux] + 1
               Rp_habi[jj,indx_Rp] = Rp_habi[jj,indx_Rp] + 1 
               flux_habi[jj,indx_flux] = flux_habi[jj,indx_flux] + 1 
      
   file.close()
   jj = jj + 1
   
   
# dump all arrays so I don't need to read them again in the future
to_save = [Mp_dist,atm_dist,atm_dist,Psurf_dist,alb_dist,relhum_dist,N2_CO2_dist,HZ_occur,planet_total,Rp_dist,flux_dist,habi_total,Rp_habi,flux_habi,n_total,f_occur]
file = open('../data/save_data.dat','w')
pickle.dump(to_save,file)
file.close()



## use this, if files were read once already.
# 
# # read pickle file:
# file = open('../data/save_data.dat','r')
# [Mp_dist,atm_dist,atm_dist,Psurf_dist,alb_dist,relhum_dist,N2_CO2_dist,HZ_occur,planet_total,Rp_dist,flux_dist,habi_total,Rp_habi,flux_habi,n_total,f_occur] = pickle.load(file)
# file.close()
# 

## plot occurrence rate as a function of Rp and flux
fig = plt.figure()
gs = gridspec.GridSpec(2, 2, width_ratios=[1, 3],height_ratios=[3,1])

# Rp distribution
to_plot = Rp_habi / Rp_dist
to_plot[Rp_dist == 0e0] = 0e0

ax0 = plt.subplot(gs[0,0])
plt.ylabel('planet radius [R$_\oplus$]')
plt.xlabel('fr. of HZ planets')
ax0.semilogy()
plt.axis([0e0,1e0,Rp_lim[0],Rp_lim[1]])

xerrmin = np.min(to_plot,axis=0)
xerrmax = np.max(to_plot,axis=0)
plt.barh(Rp_grid[:-1],xerrmax,height = Rp_grid[1:]-Rp_grid[:-1],left = 0,color='#FFF5EE',edgecolor='#FFF5EE')
# min and max values:
x = np.ravel(zip(xerrmax,xerrmax))
y = np.ravel(zip(Rp_grid[:-1],Rp_grid[1:]))
ax0.plot(x,y,'k')

# 25 and 75 percentiles
xerr25 = np.percentile(to_plot,25,axis=0)
xerr75 = np.percentile(to_plot,75,axis=0)
ax0.errorbar(np.mean(to_plot,axis=0), (Rp_grid[:-1]+Rp_grid[1:])/2e0, xerr=[np.mean(to_plot,axis=0)-xerr25,xerr75-np.mean(to_plot,axis=0)],fmt='ko',capthick=2)

plt.tick_params(length=5,width=1.5,which='major')
plt.tick_params(length=4,width=1.5,which='minor')


# flux-Rp distribution
ax1 = plt.subplot(gs[0,1])
to_plot = np.mean(habi_total,axis=0) / np.mean(planet_total,axis=0) 
to_plot[np.mean(planet_total,axis=0) == 0e0] = 0e0

#to_plot = np.min(habi_total,axis=0) / np.min(planet_total,axis=0)
#to_plot[np.mean(planet_total,axis=0) == 0e0] = 0e0
#print np.max(habi_total,axis=0)

ax1 = plt.loglog()
plt.hot()
plt.yticks(visible = False)
plt.xticks(visible = False)
plt.tick_params(axis='both',which='both',direction='out')
plt.tick_params(length=5,width=2,which='major')
plt.tick_params(length=4,width=2,which='minor')

ax1 = plt.axis([flux_grid.min(),flux_grid.max(),Rp_grid.min(),Rp_grid.max()])
ax1 = plt.pcolormesh(flux_grid,Rp_grid,to_plot,vmin=0,vmax=1) 
cbaxes = fig.add_axes([0.93, 0.42, 0.028, 0.53]) 
CB = plt.colorbar(ax1,cax=cbaxes)
CB.ax.yaxis.set_ticks_position('left')
CB.set_label('average fraction of HZ planets',color='w',labelpad=-62)
cbytick_obj = plt.getp(CB.ax.axes, 'yticklabels')
plt.setp(cbytick_obj, color='w')


# flux distribution
ax3 = plt.subplot(gs[1,1])
to_plot = flux_habi / flux_dist
to_plot[flux_dist == 0e0] = 0e0

plt.xlabel('incoming stellar flux [sol. constant]')
plt.ylabel('fr. of HZ planets')
plt.axis([flux_lim[0],flux_lim[1],0e0,1e0])
ax3.semilogx()
# min and max errors
yerrmin = np.min(to_plot,axis=0)
yerrmax = np.max(to_plot,axis=0)
ax3.bar(flux_grid[:-1],yerrmax,width = flux_grid[1:]-flux_grid[:-1],bottom = 0e0,color='#FFF5EE',edgecolor='#FFF5EE')
x = np.ravel(zip(flux_grid[:-1],flux_grid[1:]))
y = np.ravel(zip(yerrmax,yerrmax))
plt.plot(x,y,'k')

# 25 and 75 percentiles
yerr25 = np.percentile(to_plot,25,axis=0)
yerr75 = np.percentile(to_plot,75,axis=0)

ax3.errorbar(10e0**((np.log10(flux_grid[:-1])+np.log10(flux_grid[1:]))/2e0),np.mean(to_plot,axis=0), yerr=[np.mean(to_plot,axis=0)-yerr25,yerr75-np.mean(to_plot,axis=0)],fmt='ko',capthick=2)
plt.vlines([1.5,0.85,0.23, 0.2],[0e0,0e0,0e0,0e0],[1e0,1e0,1e0,1e0],colors = ['r','g','g','b'],linewidth=3,zorder=20)
plt.tight_layout(w_pad=0, h_pad=0.2)
plt.tick_params(length=5,width=1.5,which='major')
plt.tick_params(length=4,width=1.5,which='minor')

plt.savefig('../figs/habi_fraction_Rp_flux.pdf')
plt.close()


## flux and Rp distributions for the most optimistic scenario
print files[HZ_occur == np.max(HZ_occur)][0]

matplotlib.rcParams.update({'font.size': 18})


# read file:
file = open('../results/'+files[HZ_occur == np.max(HZ_occur)][0],'r')
lines = file.readlines()

n_tot = n_total * 10

# arrays for the data:
all_MS = np.zeros(n_tot)
all_flux = np.zeros(n_tot)
all_Rp = np.zeros(n_tot)
all_Mp = np.zeros(n_tot)
all_atmtype = np.chararray(n_tot) + '    '
all_Psurf = np.zeros(n_tot)
all_Tsurf = np.zeros(n_tot)
all_alb = np.zeros(n_tot)
all_relhum = np.zeros(n_tot)
all_CO2 = np.zeros(n_tot)
all_N2 = np.zeros(n_tot)
all_habi = np.chararray(n_tot) + ' '

# arrays to plot
planet_total = np.zeros([nr_bin,nr_bin])
habi_total = np.zeros([nr_bin,nr_bin])
Rp_dist = np.zeros([nr_bin])
Rp_habi = np.zeros([nr_bin])
flux_dist = np.zeros([nr_bin])
flux_habi = np.zeros([nr_bin])

j = 0
for i in lines[6:]:
   line_split = i.split(' ')
   
   all_MS[j] = float(line_split[0])
   all_flux[j] = float(line_split[1])
   all_Rp[j] = float(line_split[2])
   all_Mp[j] = float(line_split[3])
   all_atmtype[j] = line_split[4]
   all_Psurf[j] = float(line_split[5])
   all_Tsurf[j] = float(line_split[6])
   all_alb[j] = float(line_split[7])
   all_relhum[j] = float(line_split[8])
   all_CO2[j] = float(line_split[9])
   all_N2[j] = float(line_split[10])
   all_habi[j] = line_split[11]
   
   indx_Rp = bisect_left(Rp_grid,all_Rp[j]) - 1
   indx_flux = bisect_left(flux_grid,all_flux[j]) - 1
   if indx_Rp >= 0 and indx_Rp <= nr_bin-1:
      if indx_flux >= 0 and indx_flux <= nr_bin-1:
         planet_total[indx_Rp,indx_flux] = planet_total[indx_Rp,indx_flux] + 1
         Rp_dist[indx_Rp] = Rp_dist[indx_Rp] + 1 
         flux_dist[indx_flux] = flux_dist[indx_flux] + 1 
         if all_habi[j] == 'ye':
            habi_total[indx_Rp,indx_flux] = habi_total[indx_Rp,indx_flux] + 1
            Rp_habi[indx_Rp] = Rp_habi[indx_Rp] + 1 
            flux_habi[indx_flux] = flux_habi[indx_flux] + 1 
   j = j + 1


# flux-Rp distribution
to_plot = habi_total/ planet_total
to_plot[planet_total == 0e0] = 0e0

fig = plt.figure()
plt.loglog()
plt.hot()
plt.tick_params(axis='both',which='both',direction='out')
plt.tick_params(length=5,width=2,which='major')
plt.tick_params(length=4,width=2,which='minor')
plt.axis([flux_grid.min(),flux_grid.max(),Rp_grid.min(),Rp_grid.max()])
plt.xlabel('incoming stellar flux [sol. constant]')
plt.ylabel('planet radius [R$_\oplus$]')
CS = plt.pcolormesh(flux_grid,Rp_grid,to_plot,vmin=0,vmax=1) 
plt.vlines([1.5,0.85,0.23, 0.2],[Rp_grid.min(),Rp_grid.min(),Rp_grid.min(),Rp_grid.min()],[Rp_grid.max(),Rp_grid.max(),Rp_grid.max(),Rp_grid.max()],colors = ['r','g','g','b'],linewidth=4,zorder=20)
cbaxes = fig.add_axes([0.895, 0.2, 0.04, 0.7]) 
CB = plt.colorbar(CS,cax=cbaxes)
CB.ax.yaxis.set_ticks_position('left')
CB.set_label('fraction of HZ planets',color='w',labelpad=-74)
cbytick_obj = plt.getp(CB.ax.axes, 'yticklabels')
plt.setp(cbytick_obj, color='w')
plt.tight_layout(w_pad=0, h_pad=0)
plt.savefig('../figs/habi_opt_Rp_flux.pdf')
plt.close()


## flux and Rp distribution of the most pessimistic scenario
print files[HZ_occur == np.min(HZ_occur)][0]


# read file:
file = open('../results/'+files[HZ_occur == np.min(HZ_occur)][0],'r')
lines = file.readlines()

n_tot = n_total * 10

# arrays for the data:
all_MS = np.zeros(n_tot)
all_flux = np.zeros(n_tot)
all_Rp = np.zeros(n_tot)
all_Mp = np.zeros(n_tot)
all_atmtype = np.chararray(n_tot) + '    '
all_Psurf = np.zeros(n_tot)
all_Tsurf = np.zeros(n_tot)
all_alb = np.zeros(n_tot)
all_relhum = np.zeros(n_tot)
all_CO2 = np.zeros(n_tot)
all_N2 = np.zeros(n_tot)
all_habi = np.chararray(n_tot) + ' '

# arrays to plot
planet_total = np.zeros([nr_bin,nr_bin])
habi_total = np.zeros([nr_bin,nr_bin])
Rp_dist = np.zeros([nr_bin])
Rp_habi = np.zeros([nr_bin])
flux_dist = np.zeros([nr_bin])
flux_habi = np.zeros([nr_bin])

j = 0
for i in lines[6:]:
   line_split = i.split(' ')
   
   all_MS[j] = float(line_split[0])
   all_flux[j] = float(line_split[1])
   all_Rp[j] = float(line_split[2])
   all_Mp[j] = float(line_split[3])
   all_atmtype[j] = line_split[4]
   all_Psurf[j] = float(line_split[5])
   all_Tsurf[j] = float(line_split[6])
   all_alb[j] = float(line_split[7])
   all_relhum[j] = float(line_split[8])
   all_CO2[j] = float(line_split[9])
   all_N2[j] = float(line_split[10])
   all_habi[j] = line_split[11]
   
   indx_Rp = bisect_left(Rp_grid,all_Rp[j]) - 1
   indx_flux = bisect_left(flux_grid,all_flux[j]) - 1
   if indx_Rp >= 0 and indx_Rp <= nr_bin-1:
      if indx_flux >= 0 and indx_flux <= nr_bin-1:
         planet_total[indx_Rp,indx_flux] = planet_total[indx_Rp,indx_flux] + 1
         Rp_dist[indx_Rp] = Rp_dist[indx_Rp] + 1 
         flux_dist[indx_flux] = flux_dist[indx_flux] + 1 
         if all_habi[j] == 'ye':
            habi_total[indx_Rp,indx_flux] = habi_total[indx_Rp,indx_flux] + 1
            Rp_habi[indx_Rp] = Rp_habi[indx_Rp] + 1 
            flux_habi[indx_flux] = flux_habi[indx_flux] + 1
   j = j + 1


# flux-Rp distribution
to_plot = habi_total/ planet_total
to_plot[planet_total == 0e0] = 0e0

fig = plt.figure()
plt.loglog()
plt.hot()
plt.tick_params(axis='both',which='both',direction='out')
plt.tick_params(length=5,width=2,which='major')
plt.tick_params(length=4,width=2,which='minor')
plt.axis([flux_grid.min(),flux_grid.max(),Rp_grid.min(),Rp_grid.max()])
plt.xlabel('incoming stellar flux [sol. constant]')
plt.ylabel('planet radius [R$_\oplus$]')
CS = plt.pcolormesh(flux_grid,Rp_grid,to_plot,vmin=0,vmax=1) 
plt.vlines([1.5,0.85,0.23, 0.2],[Rp_grid.min(),Rp_grid.min(),Rp_grid.min(),Rp_grid.min()],[Rp_grid.max(),Rp_grid.max(),Rp_grid.max(),Rp_grid.max()],colors = ['r','g','g','b'],linewidth=4,zorder=20)
cbaxes = fig.add_axes([0.895, 0.2, 0.04, 0.7]) 
CB = plt.colorbar(CS,cax=cbaxes)
CB.ax.yaxis.set_ticks_position('left')
CB.set_label('fraction of HZ planets',color='w',labelpad=-74)
cbytick_obj = plt.getp(CB.ax.axes, 'yticklabels')
plt.setp(cbytick_obj, color='w')
plt.tight_layout(w_pad=0, h_pad=0)
plt.savefig('../figs/habi_pes_Rp_flux.pdf')
plt.close()




## study the influence of parameters on the HZ
# todo: update distributions if necessary based on MC_int.py
planet_mass = ['Mplanet_gauss_opt_sigma6','Mplanet_gauss_norm_sigma6','Mplanet_gauss_pes_sigma6','Mplanet_gauss_opt_sigma3','Mplanet_gauss_norm_sigma3','Mplanet_gauss_pes_sigma3']
atmosphere = ['atm_type_all','atm_type_N2_CO2','atm_type_H2','atm_type_N2','atm_type_CO2']
surf_pressure = ['Psurf_log_min1e0','Psurf_log_min1e2','Psurf_log_min1e4','Psurf_lognorm_1e4','Psurf_lognorm_1e5','Psurf_lognorm_1e6','Psurf_lognorm_1e7']
surf_albedo = ['alb_surf_uniform','alb_surf_normal','alb_surf_lognormal']
rel_humidity = ['relhum_log','relhum_lin','relhum_lognorm']
N2_CO2_mix_rat = ['atm_N2_CO2_log','atm_N2_CO2_lin','atm_N2_CO2_lognorm']


dbar = 0.03
dwidth = 0.03
dgroup = 0.02

fig = plt.figure()
ax = fig.add_subplot(111)
heights = [np.average(HZ_occur[Mp_dist == planet_mass[0]]),np.average(HZ_occur[Mp_dist == planet_mass[1]]),np.average(HZ_occur[Mp_dist == planet_mass[2]]),np.average(HZ_occur[Mp_dist == planet_mass[3]]),np.average(HZ_occur[Mp_dist == planet_mass[4]]),np.average(HZ_occur[Mp_dist == planet_mass[5]])]
plt.bar([0e0, dbar,2e0*dbar,3e0*dbar,4e0*dbar,5e0*dbar],heights,[dwidth,dwidth,dwidth,dwidth,dwidth,dwidth],color='#FFF5EE')

heights = [np.average(HZ_occur[atm_dist == atmosphere[0]]),np.average(HZ_occur[atm_dist == atmosphere[1]]),np.average(HZ_occur[atm_dist == atmosphere[2]]),np.average(HZ_occur[atm_dist == atmosphere[3]]),np.average(HZ_occur[atm_dist == atmosphere[4]])]
plt.bar([dwidth+5e0*dbar+dgroup,dwidth+6e0*dbar+dgroup,dwidth+7e0*dbar+dgroup,dwidth+8e0*dbar+dgroup,dwidth+9e0*dbar+dgroup],heights,[dwidth,dwidth,dwidth,dwidth,dwidth],color='#FFF5EE')

heights = [np.average(HZ_occur[Psurf_dist == surf_pressure[0]]),np.average(HZ_occur[Psurf_dist == surf_pressure[1]]),np.average(HZ_occur[Psurf_dist == surf_pressure[2]]),np.average(HZ_occur[Psurf_dist == surf_pressure[3]]),np.average(HZ_occur[Psurf_dist == surf_pressure[4]]),np.average(HZ_occur[Psurf_dist == surf_pressure[5]]),np.average(HZ_occur[Psurf_dist == surf_pressure[6]])]
plt.bar([2e0*dwidth + 9e0*dbar + 2e0*dgroup,2e0*dwidth + 10e0*dbar + 2e0*dgroup,2e0*dwidth + 11e0*dbar + 2e0*dgroup,2e0*dwidth + 12e0*dbar + 2e0*dgroup,2e0*dwidth + 13e0*dbar + 2e0*dgroup,2e0*dwidth + 14e0*dbar + 2e0*dgroup,2e0*dwidth + 15e0*dbar + 2e0*dgroup],heights,[dwidth,dwidth,dwidth,dwidth,dwidth,dwidth,dwidth],color='#FFF5EE')

heights = [np.average(HZ_occur[alb_dist == surf_albedo[0]]),np.average(HZ_occur[alb_dist == surf_albedo[1]]),np.average(HZ_occur[alb_dist == surf_albedo[2]])]
plt.bar([3e0*dwidth + 15e0*dbar + 3e0*dgroup,3e0*dwidth + 16e0*dbar + 3e0*dgroup,3e0*dwidth + 17e0*dbar + 3e0*dgroup],heights,[dwidth,dwidth,dwidth],color='#FFF5EE')

heights = [np.average(HZ_occur[relhum_dist == rel_humidity[0]]),np.average(HZ_occur[relhum_dist == rel_humidity[1]]),np.average(HZ_occur[relhum_dist == rel_humidity[2]])]
plt.bar([4e0*dwidth + 17e0*dbar + 4e0*dgroup,4e0*dwidth + 18e0*dbar + 4e0*dgroup,4e0*dwidth + 19e0*dbar + 4e0*dgroup],heights,[dwidth,dwidth,dwidth],color='#FFF5EE')

heights = [np.average(HZ_occur[N2_CO2_dist == N2_CO2_mix_rat[0]]),np.average(HZ_occur[N2_CO2_dist == N2_CO2_mix_rat[1]]),np.average(HZ_occur[N2_CO2_dist == N2_CO2_mix_rat[2]])]
plt.bar([5e0*dwidth + 19e0*dbar + 5e0*dgroup,5e0*dwidth + 20e0*dbar + 5e0*dgroup,5e0*dwidth + 21e0*dbar + 5e0*dgroup],heights,[dwidth,dwidth,dwidth],color='#FFF5EE')

plt.xticks((3e0*dbar,dwidth+7e0*dbar+dgroup+dwidth/2e0,2e0*dwidth + 12e0*dbar + 2e0*dgroup+dwidth/2e0,3e0*dwidth + 16e0*dbar + 3e0*dgroup+dwidth/2e0,4e0*dwidth + 18e0*dbar + 4e0*dgroup+dwidth/2e0,5e0*dwidth + 20e0*dbar + 5e0*dgroup+dwidth/2e0),('planet mass','atm. type','surface pressure','surf. albedo','rel. hum.','N$_2$/CO$_2$'),rotation = 20)
plt.tight_layout(w_pad=0, h_pad=0)


plt.text(0e0+dwidth/2,0.005,'PM1',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(dbar+dwidth/2,0.005,'PM2',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(2e0*dbar+dwidth/2,0.005,'PM3',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(3e0*dbar+dwidth/2,0.005,'PM4',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(4e0*dbar+dwidth/2,0.005,'PM5',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(5e0*dbar+dwidth/2,0.005,'PM6',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')

plt.text(dwidth+5e0*dbar+dgroup+dwidth/2,0.005,'AT1',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(dwidth+6e0*dbar+dgroup+dwidth/2,0.005,'AT2',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(dwidth+7e0*dbar+dgroup+dwidth/2,0.005,'AT3',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(dwidth+8e0*dbar+dgroup+dwidth/2,0.005,'AT4',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(dwidth+9e0*dbar+dgroup+dwidth/2,0.005,'AT5',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')


plt.text(2e0*dwidth+9e0*dbar+2e0*dgroup+dwidth/2,0.005,'SP1',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(2e0*dwidth+10e0*dbar+2e0*dgroup+dwidth/2,0.005,'SP2',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(2e0*dwidth+11e0*dbar+2e0*dgroup+dwidth/2,0.005,'SP3',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(2e0*dwidth+12e0*dbar+2e0*dgroup+dwidth/2,0.005,'SP4',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(2e0*dwidth+13e0*dbar+2e0*dgroup+dwidth/2,0.005,'SP5',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(2e0*dwidth+14e0*dbar+2e0*dgroup+dwidth/2,0.005,'SP6',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(2e0*dwidth+15e0*dbar+2e0*dgroup+dwidth/2,0.005,'SP7',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')

plt.text(3*dwidth+15e0*dbar+3*dgroup+dwidth/2,0.005,'SA1',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(3*dwidth+16e0*dbar+3*dgroup+dwidth/2,0.005,'SA2',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(3*dwidth+17e0*dbar+3*dgroup+dwidth/2,0.005,'SA3',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')

plt.text(4*dwidth+17e0*dbar+4*dgroup+dwidth/2,0.005,'RH1',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(4*dwidth+18e0*dbar+4*dgroup+dwidth/2,0.005,'RH2',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(4*dwidth+19e0*dbar+4*dgroup+dwidth/2,0.005,'RH3',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')

plt.text(5*dwidth+19e0*dbar+5*dgroup+dwidth/2,0.005,'MR1',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(5*dwidth+20e0*dbar+5*dgroup+dwidth/2,0.005,'MR2',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')
plt.text(5*dwidth+21e0*dbar+5*dgroup+dwidth/2,0.005,'MR3',rotation = 90,fontsize = 13,horizontalalignment='center',verticalalignment='bottom')


plt.hlines(np.average(HZ_occur),0e0,5e0*dwidth + 21e0*dbar + 5e0*dgroup+dwidth,linewidth=3)
plt.ylabel('average occurrence rate')
plt.axis([0e0,5e0*dwidth + 21e0*dbar + 5e0*dgroup+dwidth,0e0,0.2])
plt.tick_params(length=5,width=2,which='major')
plt.tick_params(length=4,width=2,which='minor')
plt.tight_layout(w_pad=0, h_pad=0)

plt.savefig('../figs/para_imp.pdf')
plt.close()

print np.average(HZ_occur)

print 'planet mass'
for i in planet_mass:
   print np.average(HZ_occur[Mp_dist == i])


print 'atmosphere'
for i in atmosphere:
   print np.average(HZ_occur[atm_dist == i])


print 'surface pressure'
for i in surf_pressure:
   print np.average(HZ_occur[Psurf_dist == i])
   
print 'surface albedo'
for i in surf_albedo:
   print np.average(HZ_occur[alb_dist == i])

print 'relative humidity'
for i in rel_humidity:
   print np.average(HZ_occur[relhum_dist == i])

print 'N2 and CO2 mixing ratios'
for i in N2_CO2_mix_rat:
   print np.average(HZ_occur[N2_CO2_dist == i])




## plot the planet mass PDFs


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

R_w = R_water(pmass)
R_s = R_silicate(pmass)
R_i = R_iron(pmass) 

Rp = 10e0**np.linspace(np.log10(0.5),np.log10(5e0),num=2000)


mean_mp_opt = np.zeros(len(Rp))
mean_mp_nom = np.zeros(len(Rp))
mean_mp_pes = np.zeros(len(Rp))

Rp_crit = 1.54

mean_mp_nom[Rp < Rp_crit] = M_silicate(Rp[Rp < Rp_crit])
mean_mp_nom[Rp >= Rp_crit] = 2.69 * Rp[Rp >= Rp_crit]**0.93
sigma_nom = mean_mp_nom / 6e0

Rp_crit = 1.22
mean_mp_pes[Rp < Rp_crit] = (M_silicate(Rp[Rp < Rp_crit]) + M_iron(Rp[Rp < Rp_crit]))/2e0
mean_mp_pes[Rp >= Rp_crit] = 2.69 * Rp[Rp >= Rp_crit]**0.93
sigma_pes = mean_mp_pes / 6e0

Rp_crit = 1.79
mean_mp_opt[Rp < Rp_crit] = (M_silicate(Rp[Rp < Rp_crit]) + M_water(Rp[Rp < Rp_crit]))/2e0
mean_mp_opt[Rp >= Rp_crit] = 2.69 * Rp[Rp >= Rp_crit]**0.93
sigma_opt = mean_mp_opt / 6e0

plt.semilogx()
plt.xlabel(r'Planet radius [R$_\oplus$]')
plt.ylabel(r'Planet mass [M$_\oplus$]')
plt.axis([0.5,5,0.1,15])

pes, = plt.plot(Rp,mean_mp_pes,'g',linewidth=2)
plt.fill_between(Rp,mean_mp_pes - sigma_pes,mean_mp_pes + sigma_pes,facecolor = 'g',alpha = 0.2)

nom, = plt.plot(Rp,mean_mp_nom,'r',linewidth=2)
plt.fill_between(Rp,mean_mp_nom - sigma_nom,mean_mp_nom + sigma_nom,facecolor = 'r',alpha = 0.2)

opt, = plt.plot(Rp,mean_mp_opt,'c',linewidth=2)
plt.fill_between(Rp,mean_mp_opt - sigma_opt,mean_mp_opt + sigma_opt,facecolor = 'c',alpha = 0.2)

water, = plt.plot(Rp[M_water(Rp[Rp < 3e0]) < 13e0],M_water(Rp[M_water(Rp[Rp < 3e0]) < 13e0]),'--b',linewidth=1.3)
Si, = plt.plot(Rp[M_silicate(Rp[Rp < 3e0]) < 13e0],M_silicate(Rp[M_silicate(Rp[Rp < 3e0]) < 13e0]),'--y',linewidth=1.3)
iron, = plt.plot(Rp[M_iron(Rp[Rp < 2e0]) < 13e0],M_iron(Rp[M_iron(Rp[Rp < 2e0]) < 13e0]),'--k',linewidth=1.3)

plt.legend([water,Si,iron,pes,nom,opt],['water planet','silicate planet','iron planet','pessimistic PDF','nominal PDF','optimistic PDF'],prop={'size':16},loc='upper left',bbox_to_anchor=(0, 0, 1, 13e0/15e0))

plt.text(0.55,14e0,'rocky planets',fontsize=16)
plt.text(1.3e0,14e0,'super-Earths',fontsize=16)
plt.text(2.5e0,14e0,'mini-Neptunes',fontsize=16)

plt.tick_params(axis='both',which='both',direction='out')
plt.tick_params(length=5,width=2,which='major')
plt.tick_params(length=4,width=2,which='minor')
plt.xticks([0.5,1,5])
plt.tight_layout(w_pad=0, h_pad=0)
plt.savefig('../figs/planet_mass_PDF.pdf')
plt.close()

