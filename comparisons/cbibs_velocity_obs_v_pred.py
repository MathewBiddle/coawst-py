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
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final_noveg'
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final_post_lee'
if dir.split("_")[-1] == 'noveg':
    run = "noveg"
elif dir.split("_")[-1] == 'lee':
    run = 'post-lee'
else:
    run = "veg"

# read CBIBS data
CBIBS_file = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/Initialization_data/CBIBS_insitu_obs/NCEI_copy/S_2011.nc'
fcbibs = netCDF4.Dataset(CBIBS_file, 'r')
# subset indexes for time. Probably don't need this after we do the searching for matching data.
itime_start = 1059
itime_end = 3618

cbibs_u = fcbibs.variables['eastward_water_velocity'][itime_start:itime_end, 2] # time,depth
cbibs_v = fcbibs.variables['northward_water_velocity'][itime_start:itime_end, 2] # time,depth
cbibs_mag = np.sqrt(cbibs_u**2 + cbibs_v**2)
cbibs_time = fcbibs.variables['time3'][itime_start:itime_end]

#Geolocation is 39.5404, -76.0736
cbibs_lon = -76.0736 # using coords from website
cbibs_lat = 39.5404

# To get all coordinates (which varies)
#cbibs_lat = fcbibs.variables['latitude'][:]
#cbibs_lon = fcbibs.variables['longitude'][:]

cbibs_date=[]
for days in cbibs_time:
    #cbibs_date.append(
    #    netCDF4.num2date(days, units=fcbibs.variables['time3'].units))
    cbibs_date.append(datetime.datetime.fromordinal(int(days)) + datetime.timedelta(days=days%1) - datetime.timedelta(days = 366))

lat_pt = cbibs_lat
lon_pt = cbibs_lon

# read ROMS history file
inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')

ocean_time = f.variables['ocean_time'][:]
lat = f.variables['lat_rho'][:]
lon = f.variables['lon_rho'][:]

## Do some date conversions ##
datetime_list=[]
for sec in ocean_time:
    datetime_list.append(netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))
# find closest point in grid
try:
    lat_pt
    lon_pt
except NameError:
    print("Using indexes x, y = (%i, %i)" % (x, y))
else:
    print("Using geo-coords lat, lon = (%f, %f)" % (lat_pt, lon_pt))
    x = coawstpy.nearest_ind(f.variables['lon_rho'][:, 1],lon_pt)
    y = coawstpy.nearest_ind(f.variables['lat_rho'][1, :],lat_pt)

ubar = f.variables['ubar_eastward'][:, x, y]
vbar = f.variables['vbar_northward'][:, x, y]
mag_vel = np.sqrt(ubar**2 + vbar**2)


# TODO index = coawstpy.nearest_ind(datelist,value)
df = pd.DataFrame(columns=['CBIBS_time','CBIBS_U','CBIBS_V','COAWST_time','COAWST_Ubar','COAWST_Vbar'],
                  index=pd.date_range(start=datetime_list[0], end=datetime_list[-1],freq='1H'))



sys.exit()
## Make some plots.
fig, (ax) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(12, 8))

ax[0].plot_date(datetime_list,mag_vel,label='COAWST',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1, color='r')
#ax.spines['left'].set_color('red')
#ax.tick_params(axis='y', colors='red')
#ax.set_ylabel('COAWST vbar [m/s]')
#ax.yaxis.label.set_color('red')

xlim = ax[0].get_xlim()
#ax.grid(True)
#ax2v = ax.twinx()
ax[0].plot_date(cbibs_date,cbibs_mag, label='CBIBS',
              xdate=True, linestyle='', linewidth=0.5,
              marker='.', markersize=1, color='b')
#ax2v.set_ylabel('CBIBS v [m/s]')
#ax2v.tick_params(axis='y', colors='blue')
#ax2v.yaxis.label.set_color('blue')
#ax2v.grid(True)
#ax2v.set_ylim([-2, 0.25])
#ax2v.set_xlim(xlim)
ax[0].xaxis.set_major_locator(mdates.DayLocator(interval=10))
ax[0].xaxis.set_major_formatter(DateFormatter("%m/%d"))

ax[0].set_ylim([0, 2.5])
ax[0].set_xlim(xlim)
ax[0].set_ylabel('Water Velocity [m/s]')
ax[0].legend(loc='upper left')
plt.suptitle('%s'%run)

ax[1].plot(mag_vel,cbibs_mag,marker='.',linestyle='')
ax[1].set_xlabel('Predicted mag [m/s]')
ax[1].set_ylabel('CBIBS mag [m/s]')


#outfile = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/figures/CBIBS_Velocity_comparison.png'

#plt.savefig(outfile, bbox_inches='tight', dpi = 1000)
