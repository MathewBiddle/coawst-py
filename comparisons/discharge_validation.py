import netCDF4
import pandas as pd
import matplotlib.pyplot as plt
import coawstpy
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
import numpy as np


# bring in the data
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final_noveg'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final_post_lee'
if dir.split("_")[-1] == 'noveg':
    run = "noveg"
elif dir.split("_")[-1] == 'lee':
    run = 'post-lee'
else:
    run = "veg"

## river data
river_frc = dir+'/river_frc.nc'
f_river = netCDF4.Dataset(river_frc, 'r')
river_time = f_river.variables['river_time'][:]
river_transport = f_river.variables['river_transport'][:, 3] # middle of water column
river_datetime_list=[]
for sec in river_time:
    river_datetime_list.append(
        netCDF4.num2date(sec, units=f_river.variables['river_time'].units, calendar='standard'))
#river_transport = f_river.variables['river_transport'][:] # (river_time, river)


inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')
ocean_time = f.variables['ocean_time'][:]
lat = f.variables['lat_rho'][:][:]
lon = f.variables['lon_rho'][:][:]
datetime_list=[]
for sec in ocean_time:
    datetime_list.append(netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))

#xlim = (datetime_list.min(), datetime_list.max())
## point selection
x = np.arange(20,30)
y = np.array([1]*len(x))
#x = 26
#y = 1
#plant_height = f.variables['plant_height'][0, 0, :, :]
#plant_height = np.ma.masked_greater(plant_height, 1)
#plant_height[x, y] = 1
#plt.figure(1)
#plt.pcolor(f.variables['lon_rho'][:][:], f.variables['lat_rho'][:][:], plant_height)


# calculate water column volume for points x and y
h = f.variables['h'][:]
zeta = f.variables['zeta'][:]
height = h[x, y]+zeta[:, x, y]
pm = f.variables['pm'][:] # width - x
pn = f.variables['pn'][:]
#coord = pm[x, y]*pn[x, y]
#vol = height/coord
sa_pt = height/pn[x, y] # surface area

#ubar = f.variables['ubar_eastward'][:]
vbar = f.variables['vbar_northward'][:] # only care about whats coming in
#ubar_pt = ubar[:, x, y]
vbar_pt = vbar[:, x, y]

#sim_q = sa*np.sqrt(vbar_pt**2+ubar_pt**2) # velocity times surface area = discharge
sim_q_pt = -1*sa_pt*vbar_pt
sim_q_mean = np.mean(sim_q_pt,axis=1)
#input_q = river_transport[:,3]
#plt.figure()
fig, ax = plt.subplots()
myFmt = DateFormatter("%m/%d")
dayint=10
ax.plot_date(river_datetime_list, river_transport, label='Input', xdate=True, linestyle='', linewidth=0.5,
                marker='.', markersize=1)
ax.plot_date(datetime_list, sim_q_mean, label='COAWST', xdate=True, linestyle='-', linewidth=0.5,
                marker='o', markersize=1)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
ax.xaxis.set_major_formatter(myFmt)
ax.set_ylabel('Discharge (m3/s)')
ax.set_xlim([min(datetime_list),max(datetime_list)])
ax.legend()

#plt.xlim(xlim)
