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
import pandas as pd

# bring in the data
#dir='/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110702_20111101'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110906_20110926'
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
#dir = '/Volumes/Documents/COAWST_34_UPPER_CHES_FULL'
#inputfile = dir+'/upper_ches_his.nc'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110714201800_20111031231800'
inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')

eotb_file = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/Initialization_data/Eyes_on_the_bay/EOTBData_Susquehanna_01Jul11_TO_01Nov11.csv'
feotb = pd.read_csv(eotb_file, header=0, parse_dates=['DateTime'], infer_datetime_format=True)

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


#sys.exit()
eotb_turb = feotb['Turb_NTU'].values
#eotb_turb = np.ma.masked_equal(eotb_turb, eotb_turb.min())
feotb['DateTimeUTC'] = feotb['DateTime'] #+pd.DateOffset(hours=4)
eotb_date = feotb['DateTimeUTC'].values
#Geolocation is 39.5404, -76.0736
eotb_lon = -76.0848
eotb_lat = 39.5478
#eotb_lat = feotb.variables['latitude'][:]
#eotb_lon = feotb.variables['longitude'][:]

# epoch_date = '%s %s'%(feotb.variables['time2'].units.split(' ')[-3], feotb.variables['time2'].units.split(' ')[-2])
# dt_obj = datetime.datetime.strptime(epoch_date, '%Y-%m-%d %H:%M:%S')
# time_diff = time.mktime(dt_obj.timetuple())-time.mktime(datetime.datetime(1970, 1, 1, 0, 0, 0).timetuple())
#eotb_date=feotb['DateTime'].values
# for sec in eotb_time:
#     ts = sec/(12*3600)
#     eotb_date.append(datetime.datetime.fromtimestamp(sec+time_diff))

lat_pt = eotb_lat
lon_pt = eotb_lon

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
mud = np.mean(mud_tot, axis=1)  # integrate w/ depth
sand_tot = f.variables['sand_01'][:, :, x, y]
#sand = sand_tot
sand = np.mean(sand_tot, axis=1)  # integrate w/ depth

SSC = mud + sand
fig, (ax) = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(12, 8))
ax.plot_date(datetime_list,sand,label='COAWST sand',
              xdate=True, linestyle='-', linewidth=1,
              marker='', markersize=1, color='r')
ax.plot_date(datetime_list,mud,label='COAWST mud',
              xdate=True, linestyle='-', linewidth=1,
              marker='', markersize=1, color='g')
#ax.spines['left'].set_color('red')
ax.tick_params(axis='y', colors='red')
ax.set_ylabel('COAWST SSC [kg/m3]')
ax.yaxis.label.set_color('red')
xlim = ax.get_xlim()
#ax.grid(True)
ax2v = ax.twinx()
ax2v.plot_date(eotb_date,eotb_turb, label='Havre de Grace',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1, color='b')
ax2v.set_ylabel('Havre de Grace Turbidity [NTU]')
ax2v.tick_params(axis='y', colors='blue')
ax2v.yaxis.label.set_color('blue')
#ax2v.grid(True)
ax2v.set_xlim(xlim)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=10))
ax.xaxis.set_major_formatter(DateFormatter("%m/%d"))
ax.legend(loc='upper left')

