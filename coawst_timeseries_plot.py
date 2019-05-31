#%tb
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
#import numpy.ma as ma
import datetime
import time
import netCDF4
import sys
#import os

# point location
x = 57
y = 27

dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110702_20111101'
inputfile = dir+'/upper_ches_his.nc'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110714201800_20111031231800'
#inputfile = dir+'/upper_ches_his_tide.nc'

f = netCDF4.Dataset(inputfile,'r')

ocean_time = f.variables['ocean_time'][:]
#lat = f.variables['lat_rho'][:][:]
#lon = f.variables['lon_rho'][:][:]
h = f.variables['h'][:][:]

## Do some date conversions ##
epoch_date = '%s %s'%(f.variables['ocean_time'].units.split(' ')[-2], f.variables['ocean_time'].units.split(' ')[-1])
dt_obj = datetime.datetime.strptime(epoch_date, '%Y-%m-%d %H:%M:%S')
time_diff = time.mktime(dt_obj.timetuple())-time.mktime(datetime.datetime(1970, 1, 1, 0, 0, 0).timetuple())
datetime_list=[]
for sec in ocean_time:
    ts = sec/(12*3600)
    datetime_list.append(datetime.datetime.fromtimestamp(sec+time_diff))

# Verify point location
h[x,y] = 100
plt.figure(1)
plt.pcolor(f.variables['lon_rho'][:][:],f.variables['lat_rho'][:][:],h)


# plot variables as time series
myFmt = DateFormatter("%m/%d")

var2plot = ['zeta','rho','Hwave','tke','zetaw']

fig, (ax) = plt.subplots(nrows=len(var2plot), ncols=1, sharex=True, figsize=(12, 8))
fig.subplots_adjust(hspace=0.05)

for i, ax in enumerate(fig.axes):
    if (var2plot[i] == 'rho') or (var2plot[i] == 'salt') or (var2plot[i] == 'tke'):
        ax.plot_date(datetime_list, f.variables[var2plot[i]][:, 0, x, y], xdate=True, linestyle='-', linewidth=1,
                     marker='', markersize=1)
    else:
        ax.plot_date(datetime_list, f.variables[var2plot[i]][:, x, y], xdate=True, linestyle='-', linewidth=1,
                     marker='', markersize=1)
    ax.set_ylabel('%s' % (f.variables[var2plot[i]].name))#, f.variables[var2plot[i]].units))
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=10))
    ax.xaxis.set_major_formatter(myFmt)
    ax.grid(True)


def stick_plot(time, u, v, **kw):
    width = kw.pop('width', 0.001)
    headwidth = kw.pop('headwidth', 0)
    headlength = kw.pop('headlength', 0)
    headaxislength = kw.pop('headaxislength', 0)
    angles = kw.pop('angles', 'uv')
    ax = kw.pop('ax', None)

    if angles != 'uv':
        raise AssertionError("Stickplot angles must be 'uv' so that"
                             "if *U*==*V* the angle of the arrow on"
                             "the plot is 45 degrees CCW from the *x*-axis.")

    time, u, v = map(np.asanyarray, (time, u, v))
    if not ax:
        fig, ax = plt.subplots()

    q = ax.quiver(mdates.date2num(time), [[0] * len(time)], u, v,
                  angles='uv', width=width, headwidth=headwidth,
                  headlength=headlength, headaxislength=headaxislength,
                  **kw)

    ax.axes.get_yaxis().set_visible(False)
    ax.xaxis_date()
    return q

q = stick_plot(datetime_list, f.variables['ubar_eastward'][:,x,y], f.variables['vbar_northward'][:,x,y])