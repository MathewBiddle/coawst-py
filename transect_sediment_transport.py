import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
#import warnings
import numpy as np
import datetime
import time
import netCDF4
import coawstpy
#import pandas as pd

# bring in the data
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')

ocean_time = f.variables['ocean_time'][:]
lat = f.variables['lat_rho'][:][:]
lon = f.variables['lon_rho'][:][:]

## Do some date conversions ##
datetime_list=[]
for sec in ocean_time:
    datetime_list.append(
        netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))


# figure out the transect:

# identify the points to extract data
# Susquehanna River mouth
trans_name_sr = 'T1'
x_sr = np.array([29,30,31,32])
y_sr = np.array([13,12,11,10])

# Turkey Point to Sandy Point
#trans_name = 'Turkey Pt to Sandy Pt'
#x = np.array(list(range(42,67)))  #
#y = np.array([58]*len(x))

# Verify point location
plant_height = f.variables['plant_height'][0, 0, :, :]
plant_height = np.ma.masked_greater(plant_height, 1)
plant_height[x_sr, y_sr] = 1
# plt.figure(1)
# plt.pcolor(f.variables['lon_rho'][:][:], f.variables['lat_rho'][:][:], plant_height)
# plt.title('%s Transect' % trans_name)

# calculate the total concentration at each time step for the transect (x,y) and depth.
# (h+zeta)/(pm*pn) * mud_01_depth_sum(time,x,y)
# sum over depth
mud_01 = f.variables['mud_01'][:]#, :, x, y]
mud_01_trans_sr = mud_01[:, :, x_sr, y_sr]
mud_sr = np.sum(mud_01_trans_sr, axis=1)

sand_01 = f.variables['sand_01'][:]
sand_01_trans_sr = sand_01[:, :, x_sr, y_sr]
sand_sr = np.sum(sand_01_trans_sr, axis=1)

# calculate water column volume for points x and y
h = f.variables['h'][:]
zeta = f.variables['zeta'][:]
height_sr = h[x_sr, y_sr]+zeta[:, x_sr, y_sr]

#sys.exit()
pm = f.variables['pm'][:]
pn = f.variables['pn'][:]
coord_sr = pm[x_sr, y_sr]*pn[x_sr, y_sr]
vol_sr = height_sr/coord_sr
#sys.exit()
mud_trans_conc_sr = vol_sr*mud_sr
mud_tot_conc_sr = np.sum(mud_trans_conc_sr, axis=1)/1000  # metric ton

sand_trans_conc_sr = vol_sr*sand_sr
sand_tot_conc_sr = np.sum(sand_trans_conc_sr, axis=1)/1000  # metric ton

total_sediment_sr = np.float(np.sum(sand_tot_conc_sr)+np.sum(mud_tot_conc_sr))  # metric ton
print("Transect at the mouth of S.R. total suspended sediment mass = %e tons" % total_sediment_sr)
# totals to 4.62 x 10^6 metric tons entering the bay
#
# From Palinkas 2014:
# Various estimates have been generated for the amount of sediment delivered to the Chesapeake Bay associated with
# TS Lee, ranging from 6.7 to 19.0 X 10^6 t(Cheng et al., 2013; Hirsch, 2012),
#
# Based upon the model simulation, we estimated that approximately 18 days after the flood (28 September 2011), nearly
# 6.73 X 10^6 t of fluvial sediment deposited inside the Bay, 0.15?106

#sys.exit()
#vel_mag = np.sqrt(ubar**2 + vbar**2)

# initialize the figure
myFmt = DateFormatter("%m/%d")
dayint = 10
fig, (ax) = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(12, 8))
fig.subplots_adjust(hspace=0.25)

# start plotting
ax[0].plot_date(datetime_list,mud_tot_conc_sr, label='mud_01',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1, color='r')
ax[0].plot_date(datetime_list,sand_tot_conc_sr,label='sand_01',
              xdate=True, linestyle='--', linewidth=0.5,
              marker='', markersize=1, color='r')
ax[0].set_ylabel('%s sed [t]'%trans_name_sr)
ax[0].xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
ax[0].xaxis.set_major_formatter(myFmt)
#ax[0].grid(True)
ax[0].legend(loc='upper left')
ax[2].plot_date(datetime_list,(mud_tot_conc_sr+sand_tot_conc_sr), label=trans_name_sr,
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1, color='r')


# Turkey Point to Sandy Point
trans_name_t2s = 'T2'
x_t2s = np.array(list(range(42,67)))  #
y_t2s = np.array([58]*len(x_t2s))

# verify location
plant_height[x_t2s, y_t2s] = 2


# calculate the total concentration at each time step for the transect (x,y) and depth.
# (h+zeta)/(pm*pn) * mud_01_depth_sum(time,x,y)
# sum over depth

mud_01_trans_t2s = mud_01[:, :, x_t2s, y_t2s]
mud_t2s = np.sum(mud_01_trans_t2s, axis=1)

sand_01_trans_t2s = sand_01[:, :, x_t2s, y_t2s]
sand_t2s = np.sum(sand_01_trans_t2s, axis=1)

# calculate water column volume for points x and y
height_t2s = h[x_t2s, y_t2s]+zeta[:, x_t2s, y_t2s]

#sys.exit()
coord_t2s = pm[x_t2s, y_t2s]*pn[x_t2s, y_t2s]
vol_t2s = height_t2s/coord_t2s
#sys.exit()
mud_trans_conc_t2s = vol_t2s*mud_t2s
mud_tot_conc_t2s = np.sum(mud_trans_conc_t2s, axis=1)/1000  # metric ton

sand_trans_conc_t2s = vol_t2s*sand_t2s
sand_tot_conc_t2s = np.sum(sand_trans_conc_t2s, axis=1)/1000  # metric ton

total_sediment_t2s = (np.sum(sand_tot_conc_t2s)+np.sum(mud_tot_conc_t2s))  # metric ton
print("Transect Turkey Point to Sandy Point total suspended sediment mass = %e tons" % total_sediment_t2s)
# totals to 4.62 x 10^6 metric tons exiting the bay
#
# From Palinkas 2014:
# Various estimates have been generated for the amount of sediment delivered to the Chesapeake Bay associated with
# TS Lee, ranging from 6.7 to 19.0 X 10^6 t(Cheng et al., 2013; Hirsch, 2012),
#
# Based upon the model simulation, we estimated that approximately 18 days after the flood (28 September 2011), nearly
# 6.73 X 10^6 t of fluvial sediment deposited inside the Bay, 0.15?106

#sys.exit()
#vel_mag = np.sqrt(ubar**2 + vbar**2)

myFmt = DateFormatter("%m/%d")
dayint = 10

#fig, (ax) = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(12, 8))
#fig.subplots_adjust(hspace=0.05)
ax[1].plot_date(datetime_list,mud_tot_conc_t2s, label='mud_01',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1, color='b')
ax[1].plot_date(datetime_list,sand_tot_conc_t2s,label='sand_01',
              xdate=True, linestyle='--', linewidth=0.5,
              marker='', markersize=1, color='b')
ax[1].set_ylabel('%s sed [t]'%trans_name_t2s)
ax[1].xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
ax[1].xaxis.set_major_formatter(myFmt)
#ax[1].grid(True)
ax[1].legend(loc='upper left')
#ylim = ax[1].get_ylim()

ax[2].plot_date(datetime_list,(mud_tot_conc_t2s+sand_tot_conc_t2s), label=trans_name_t2s,
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1, color='b')
ax[2].legend(loc='upper left')
ax[2].set_ylabel('sand+mud [t]')
plt.suptitle('Total Suspended Sediment Mass')
plt.figure()
plt.pcolor(f.variables['lon_rho'][:][:], f.variables['lat_rho'][:][:], plant_height)
plt.title('Transect Locations')