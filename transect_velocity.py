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
#dir='/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110702_20111101'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110906_20110926'
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
#dir = '/Volumes/Documents/COAWST_34_UPPER_CHES_FULL'
#inputfile = dir+'/upper_ches_his.nc'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110714201800_20111031231800'
inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')

#roms_in_dir = '/Volumes/Documents/ROMS_INPUT_FILES'
river_frc = dir+'/river_frc.nc'

f_river = netCDF4.Dataset(river_frc, 'r')
river_time = f_river.variables['river_time'][:]
river_transport = f_river.variables['river_transport'][:, 0]

## Do some date conversions ##
#epoch_date = '%s %s'%(f_river.variables['river_time'].units.split(' ')[-2], f_river.variables['river_time'].units.split(' ')[-1])
#dt_obj = datetime.datetime.strptime(epoch_date, '%Y-%m-%d %H:%M:%S')
#time_diff = time.mktime(dt_obj.timetuple())-time.mktime(datetime.datetime(1970, 1, 1, 0, 0, 0).timetuple())
river_datetime_list=[]
for sec in river_time:
    #ts = sec/(12*3600)
    river_datetime_list.append(
        netCDF4.num2date(sec, units=f_river.variables['river_time'].units))
    #river_datetime_list.append(datetime.datetime.fromtimestamp(sec+time_diff))

#sys.exit()
ocean_time = f.variables['ocean_time'][:]
lat = f.variables['lat_rho'][:][:]
lon = f.variables['lon_rho'][:][:]

## Do some date conversions ##
#epoch_date = '%s %s'%(f.variables['ocean_time'].units.split(' ')[-2], f.variables['ocean_time'].units.split(' ')[-1])
#dt_obj = datetime.datetime.strptime(epoch_date, '%Y-%m-%d %H:%M:%S')
#time_diff = time.mktime(dt_obj.timetuple())-time.mktime(datetime.datetime(1970, 1, 1, 0, 0, 0).timetuple())
datetime_list=[]
for sec in ocean_time:
   # ts = sec/(12*3600)
    datetime_list.append(
        netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))
    #datetime_list.append(datetime.datetime.fromtimestamp(sec+time_diff))


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

# get the data for the transect and take the average across the transect
# identify axis=1 to get the average across the x's, then y's.
ubar = np.mean(np.mean(f.variables['ubar_eastward'][:, x, y], axis=1), axis=1)
vbar = np.mean(np.mean(f.variables['vbar_northward'][:, x, y], axis=1), axis=1)
vel_mag = np.sqrt(ubar**2 + vbar**2)

myFmt = DateFormatter("%m/%d")
dayint = 10

fig, (ax) = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(12, 8))
fig.subplots_adjust(hspace=0.05)
ax[0].plot_date(datetime_list,ubar, label='East',
              xdate=True, linestyle=':', linewidth=0.5,
              marker='', markersize=1, color='r')
ax[0].plot_date(datetime_list,vbar,label='North',
              xdate=True, linestyle='--', linewidth=0.5,
              marker='', markersize=1, color='r')
ax[0].set_ylabel('%s velocity [m/s]'%trans_name)
ax[0].xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
ax[0].xaxis.set_major_formatter(myFmt)
ax[0].grid(True)
ax[0].legend(loc='lower left')
ylim = ax[0].get_ylim()

ax[2].plot_date(datetime_list,vel_mag,label=trans_name,
                xdate=True, linestyle='-', linewidth=0.5,
                marker='', markersize=1, color='r')
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

ubar = np.mean(np.mean(f.variables['ubar_eastward'][:, x, y], axis=1), axis=1)
vbar = np.mean(np.mean(f.variables['vbar_northward'][:, x, y], axis=1), axis=1)
vel_mag = np.sqrt(ubar**2 + vbar**2)

ax[1].plot_date(datetime_list,ubar,label='East',
              xdate=True, linestyle=':', linewidth=0.5,
              marker='', markersize=1, color='b')
ax[1].plot_date(datetime_list,vbar,label='North',
              xdate=True, linestyle='--', linewidth=0.5,
              marker='', markersize=1, color='b')
ax[1].set_ylabel('%s velocity [m/s]' %trans_name)
ax[1].xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
ax[1].xaxis.set_major_formatter(myFmt)
ax[1].set_ylim(ylim)
ax[1].grid(True)
ax[1].legend(loc='lower left')

ax[2].plot_date(datetime_list,vel_mag,label=trans_name,
                xdate=True, linestyle='-', linewidth=0.5,
                marker='', markersize=1, color='b')
ax[2].legend(loc='upper left')
ax[2].set_ylabel('velocity magnitude [m/s]')
ax[2].grid(True)
# q = coawstpy.stick_plot(datetime_list, ubar, vbar, ax=ax[1])
# ax[1].set_ylabel(trans_name)
# ax[1].xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
# ax[1].xaxis.set_major_formatter(myFmt)
# ax[1].grid(True)
# ax[1].quiverkey(q, X=0.1, Y=0.85, U=1, label='%s 1 m/s' % trans_name, labelpos='N', coordinates='axes')

plt.figure()
plt.pcolor(f.variables['lon_rho'][:][:], f.variables['lat_rho'][:][:], plant_height)
plt.title('Transect Locations')