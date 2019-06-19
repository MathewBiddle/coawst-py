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
trans_name = 'S.R. Mouth'
x = np.array([29,30,31,32])
y = np.array([13,12,11,10])

# Turkey Point to Sandy Point
#trans_name = 'Turkey Pt to Sandy Pt'
#x = np.array(list(range(42,67)))  #
#y = np.array([58]*len(x))

try:
    lat_pt
    lon_pt
except NameError:
    print("Using indexes x, y = (%s, %s)" % (x, y))
else:
    print("Using geo-coords lat, lon = (%f, %f)" % (lat_pt, lon_pt))
    x = np.abs(f.variables['lon_rho'][:, 1]-lon_pt).argmin()
    y = np.abs(f.variables['lat_rho'][1, :]-lat_pt).argmin()

# Verify point location
plant_height = f.variables['plant_height'][0, 0, :, :]
plant_height = np.ma.masked_greater(plant_height, 1)
plant_height[x, y] = 1
# plt.figure(1)
# plt.pcolor(f.variables['lon_rho'][:][:], f.variables['lat_rho'][:][:], plant_height)
# plt.title('%s Transect' % trans_name)

# calculate the total concentration at each time step for the transect (x,y) and depth.
# (h+zeta)/(pm*pn) * mud_01_depth_sum(time,x,y)
# sum over depth
mud = np.sum(f.variables['mud_01'][:, :, x, y], axis=1)
sand = np.sum(f.variables['sand_01'][:, :, x, y], axis=1)
# calculate water column volume for points x and y
height = f.variables['h'][x, y]+f.variables['zeta'][:, x, y]
coord = f.variables['pm'][x, y]*f.variables['pn'][x, y]
vol = height/coord

mud_trans_conc = vol*mud
mud_tot_conc = np.sum(np.sum(mud_trans_conc, axis=1), axis=1)  # kg

sand_trans_conc = vol*sand
sand_tot_conc = np.sum(np.sum(sand_trans_conc, axis=1), axis=1)  # kg

total_sediment = (np.sum(sand_tot_conc)+np.sum(mud_tot_conc))/1000  # metric ton
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

myFmt = DateFormatter("%m/%d")
dayint = 10

fig, (ax) = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(12, 8))
fig.subplots_adjust(hspace=0.5)
ax[0].plot_date(datetime_list,mud_tot_conc, label='mud_01',
              xdate=True, linestyle=':', linewidth=0.5,
              marker='', markersize=1, color='r')
ax[0].plot_date(datetime_list,sand_tot_conc,label='sand_01',
              xdate=True, linestyle='--', linewidth=0.5,
              marker='', markersize=1, color='r')
ax[0].set_ylabel('%s sed [kg]'%trans_name)
ax[0].xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
ax[0].xaxis.set_major_formatter(myFmt)
ax[0].grid(True)
ax[0].legend(loc='upper left')
ax[2].plot_date(datetime_list,(mud_tot_conc+sand_tot_conc)/1000, label='S.R.',
              xdate=True, linestyle='--', linewidth=0.5,
              marker='', markersize=1, color='r')
#ylim = ax[0].get_ylim()
#sys.exit()
#ax[2].plot_date(datetime_list,vel_mag,label=trans_name,
#                xdate=True, linestyle='-', linewidth=0.5,
#                marker='', markersize=1, color='r')
# q = coawstpy.stick_plot(datetime_list, ubar, vbar, ax=ax[0])
# ax[0].set_ylabel(trans_name)
# ax[0].xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
# ax[0].xaxis.set_major_formatter(myFmt)
# ax[0].grid(True)
# ax[0].quiverkey(q, X=0.1, Y=0.85, U=1, label='%s 1 m/s' % trans_name, labelpos='N', coordinates='axes')

# Turkey Point to Sandy Point
trans_name = 'Turkey2Sandy'
x = np.array(list(range(42,67)))  #
y = np.array([58]*len(x))

# verify location
plant_height[x, y] = 1

# calculate the total concentration at each time step for the transect (x,y) and depth.
# (h+zeta)/(pm*pn) * mud_01_depth_sum(time,x,y)
# sum over depth
mud = np.sum(f.variables['mud_01'][:, :, x, y], axis=1)
sand = np.sum(f.variables['sand_01'][:, :, x, y], axis=1)
# calculate water column volume for points x and y
height = f.variables['h'][x, y]+f.variables['zeta'][:, x, y]
coord = f.variables['pm'][x, y]*f.variables['pn'][x, y]
vol = height/coord

mud_trans_conc = vol*mud
mud_tot_conc = np.sum(np.sum(mud_trans_conc, axis=1), axis=1)  # kg

sand_trans_conc = vol*sand
sand_tot_conc = np.sum(np.sum(sand_trans_conc, axis=1), axis=1)  # kg

total_sediment = (np.sum(sand_tot_conc)+np.sum(mud_tot_conc))/1000  # metric ton
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
ax[1].plot_date(datetime_list,mud_tot_conc, label='mud_01',
              xdate=True, linestyle=':', linewidth=0.5,
              marker='', markersize=1, color='b')
ax[1].plot_date(datetime_list,sand_tot_conc,label='sand_01',
              xdate=True, linestyle='--', linewidth=0.5,
              marker='', markersize=1, color='b')
ax[1].set_ylabel('%s sed [kg]'%trans_name)
ax[1].xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
ax[1].xaxis.set_major_formatter(myFmt)
ax[1].grid(True)
ax[1].legend(loc='upper left')
#ylim = ax[1].get_ylim()

ax[2].plot_date(datetime_list,(mud_tot_conc+sand_tot_conc)/1000, label='Turkey Pt',
              xdate=True, linestyle='--', linewidth=0.5,
              marker='', markersize=1, color='b')
ax[2].legend(loc='upper left')
ax[2].set_ylabel('Total sediment [t]')

plt.figure()
plt.pcolor(f.variables['lon_rho'][:][:], f.variables['lat_rho'][:][:], plant_height)
plt.title('Transect Locations')