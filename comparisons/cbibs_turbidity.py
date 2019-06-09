#import os
#os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
#from mpl_toolkits.basemap import Basemap
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

CBIBS_file = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/Initialization_data/CBIBS_insitu_obs/NCEI_copy/S_2011.nc'
fcbibs = netCDF4.Dataset(CBIBS_file, 'r')

ocean_time = f.variables['ocean_time'][:]
lat = f.variables['lat_rho'][:][:]
lon = f.variables['lon_rho'][:][:]

## Do some date conversions ##
#epoch_date = '%s %s'%(f.variables['ocean_time'].units.split(' ')[-2], f.variables['ocean_time'].units.split(' ')[-1])
#dt_obj = datetime.datetime.strptime(epoch_date, '%Y-%m-%d %H:%M:%S')
#time_diff = time.mktime(dt_obj.timetuple())-time.mktime(datetime.datetime(1970, 1, 1, 0, 0, 0).timetuple())
datetime_list=[]
for sec in ocean_time:
    datetime_list.append(
        netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))
    #ts = sec/(12*3600)
    #datetime_list.append(datetime.datetime.fromtimestamp(sec+time_diff))



cbibs_turb = fcbibs.variables['turbidity'][:]
cbibs_turb = np.ma.masked_equal(cbibs_turb, cbibs_turb.min())
cbibs_time = fcbibs.variables['time2'][:]
#Geolocation is 39.5404, -76.0736
cbibs_lon = -76.0736
cbibs_lat = 39.5404
#cbibs_lat = fcbibs.variables['latitude'][:]
#cbibs_lon = fcbibs.variables['longitude'][:]

cbibs_date=[]
for i in range(len(cbibs_time)):
    cbibs_date.append(datetime.datetime.fromordinal(int(cbibs_time[i])) + datetime.timedelta(days=cbibs_time[i]%1) - datetime.timedelta(days = 366))

# epoch_date = '%s %s'%(fcbibs.variables['time2'].units.split(' ')[-3], fcbibs.variables['time2'].units.split(' ')[-2])
# dt_obj = datetime.datetime.strptime(epoch_date, '%Y-%m-%d %H:%M:%S')
# time_diff = time.mktime(dt_obj.timetuple())-time.mktime(datetime.datetime(1970, 1, 1, 0, 0, 0).timetuple())
# cbibs_date=[]
# for sec in cbibs_time:
#     ts = sec/(12*3600)
#     cbibs_date.append(datetime.datetime.fromtimestamp(sec+time_diff))

lat_pt = cbibs_lat
lon_pt = cbibs_lon

try:
    lat_pt
    lon_pt
except NameError:
    print("Using indexes x, y = (%i, %i)" % (x, y))
else:
    print("Using geo-coords lat, lon = (%f, %f)" % (lat_pt, lon_pt))
    x = np.abs(f.variables['lon_rho'][:, 1]-lon_pt).argmin()
    y = np.abs(f.variables['lat_rho'][1, :]-lat_pt).argmin()


#plant_height = f.variables['plant_height'][0, 0, :, :]
#plant_height = np.ma.masked_greater(plant_height, 1)
#plant_height[x, y] = 1
#plt.figure(1)
#plt.pcolor(f.variables['lon_rho'][:][:], f.variables['lat_rho'][:][:], plant_height)


mud_tot = f.variables['mud_01'][:, :, x, y]
#mud = mud_tot
mud = np.mean(mud_tot, axis=1)
sand_tot = f.variables['sand_01'][:, :, x, y]
#sand = sand_tot
sand = np.mean(sand_tot, axis=1)

SSC = mud + sand
fig, (ax) = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(12, 8))
ax.plot_date(datetime_list,SSC,label='COAWST SSC',
              xdate=True, linestyle='-', linewidth=1,
              marker='', markersize=1, color='r')
#ax.spines['left'].set_color('red')
ax.tick_params(axis='y', colors='red')
ax.set_ylabel('COAWST SSC [kg/m3]')
ax.yaxis.label.set_color('red')
xlim = ax.get_xlim()
#ax.grid(True)
ax2v = ax.twinx()
ax2v.plot_date(cbibs_date,cbibs_turb, label='CBIBS',
              xdate=True, linestyle='-', linewidth=1,
              marker='', markersize=1, color='b')
ax2v.set_ylabel('CBIBS Turbidity [NTU]')
ax2v.tick_params(axis='y', colors='blue')
ax2v.yaxis.label.set_color('blue')
#ax2v.grid(True)
ax2v.set_xlim(xlim)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=10))
ax.xaxis.set_major_formatter(DateFormatter("%m/%d"))
#fig.legend(loc='upper left')

