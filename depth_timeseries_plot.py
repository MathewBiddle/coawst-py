import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import numpy as np
import netCDF4

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

eotb_lon = -76.0848
eotb_lat = 39.5478

lat_pt = eotb_lat
lon_pt = eotb_lon

x = np.abs(f.variables['lon_rho'][:, 1] - lon_pt).argmin()
y = np.abs(f.variables['lat_rho'][1, :] - lat_pt).argmin()
z_rho = f.variables['z_rho'][:]
zeta = f.variables['zeta'][:]

fig, (ax) = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(12, 8))
for i in range(len(z_rho[0, :, x, y])):
    ax.plot_date(datetime_list, z_rho[:, i, x, y], label='s_rho=%i' % i,
              xdate=True, linestyle='-', linewidth=1,
              marker='', markersize=1)
ax.plot_date(datetime_list, zeta[:, x, y], label='zeta',
              xdate=True, linestyle='-', linewidth=1,
              marker='', markersize=1)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=10))
ax.xaxis.set_major_formatter(DateFormatter("%m/%d"))
ax.legend(loc='upper left', fontsize='x-small')

