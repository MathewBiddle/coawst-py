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
import pandas as pd
#import sys
#import os

## point location geo
#lat_pt = 39.516345
#lon_pt = -76.034109

# Cindy's locations
locs = pd.DataFrame(columns=['Site', 'lat', 'lon'])
locs['Site'] = ['1','2','3','4','5',
                'Lee7','Lee6','Lee5','Lee2.5','Lee2','Lee0','LeeS2',
                'CBIBS','Tripod']
locs['lat'] = [39.527,39.533,39.515,39.505,39.497,
               39.414,39.380,39.346,39.197,39.135,39.061,38.757,
               39.5396,39.4931]
locs['lon'] = [-76.061,-76.061,-76.051,-76.039,-76.036,
               -76.079,-76.088,-76.197,-76.311,-76.328,-76.328,-76.473,
               -76.0741,-76.0341]
locs['comment'] = ['Russ and Palinkas 2018','','Middle of bed','','',
                   '','','','','','','',
                   'CBIBS Susquehanna Flats','Larry tripod site']

# get coords for Site of choice
lat_pt, lon_pt = locs.loc[locs['Site'] == 'CBIBS', ['lat', 'lon']].values[0]

## Using point location idx
#x = 56   # 57
#y = 27      #

#dir='/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110702_20111101'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110906_20110926'
dir = '/Volumes/Documents/COAWST_34_UPPER_CHES_FULL'
#inputfile = dir+'/upper_ches_his.nc'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110714201800_20111031231800'
inputfile = dir+'/upper_ches_his.nc'

f = netCDF4.Dataset(inputfile, 'r')

try:
    lat_pt
    lon_pt
except NameError:
    print("Using indexes x, y = (%i, %i)" % (x, y))
else:
    print("Using geo-coords lat, lon = (%f, %f)" % (lat_pt, lon_pt))
    x = np.abs(f.variables['lon_rho'][:, 1]-lon_pt).argmin()
    y = np.abs(f.variables['lat_rho'][1, :]-lat_pt).argmin()
#sys.exit()
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

var2plot = ['zeta','velocity','Hwave','tke','mud+sand','bed_thickness']

fig, (ax) = plt.subplots(nrows=len(var2plot), ncols=1, sharex=True, figsize=(12, 8))
fig.subplots_adjust(hspace=0.05)
dayint = 10
for i, ax in enumerate(fig.axes):
    if (var2plot[i] == 'rho') or (var2plot[i] == 'salt') or (var2plot[i] == 'tke') or (var2plot[i] == 'mud_01') or \
            (var2plot[i] == 'sand_01') or (var2plot[i] == 'temp'):
        ax.plot_date(datetime_list, f.variables[var2plot[i]][:, 0, x, y], xdate=True, linestyle='-', linewidth=1,
                     marker='.', markersize=1)
        ax.set_ylabel('%s' % (f.variables[var2plot[i]].name))  # , f.variables[var2plot[i]].units))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
        ax.xaxis.set_major_formatter(myFmt)
        ax.grid(True)

    elif var2plot[i] == 'mud+sand':
        ax.plot_date(datetime_list, f.variables['mud_01'][:, 0, x, y], xdate=True, linestyle='-',linewidth=1,
                     marker='',markersize=1,color='b')
        ax.set_ylabel('%s' % (f.variables['mud_01'].name),color='b')  #
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
        ax.xaxis.set_major_formatter(myFmt)

        ylims = ax.get_ylim()
        ax2v = ax.twinx()
        ax2v.plot_date(datetime_list, f.variables['sand_01'][:, 0, x, y], xdate=True, linestyle='-', linewidth=1,
                     marker='',markersize=1,color='r')
        ax2v.set_ylabel('%s' % ('sand_01'),color='r')
        ax2v.set_ylim(ylims)
        ax2v.xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
        ax2v.xaxis.set_major_formatter(myFmt)
        ax.grid(True)
    elif var2plot[i] == 'velocity':
        #ax.plot_date(datetime_list,
        #             np.sqrt(np.add(f.variables['ubar_eastward'][:, x, y]**2, f.variables['vbar_northward'][:, x, y]**2)),
        #             xdate=True, linestyle='-', linewidth=1, marker='', markersize=1)
        ax.plot_date(datetime_list,f.variables['ubar_eastward'][:, x, y],
                     xdate=True, linestyle='-', linewidth=1, marker='', markersize=1, color='b')
        ax.set_ylabel('%s' % (f.variables['ubar_eastward'].name), color='b')  #
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_formatter(myFmt)

        ax2v = ax.twinx()
        ax2v.plot_date(datetime_list,f.variables['vbar_northward'][:, x, y],
                       xdate=True, linestyle='-', linewidth=1, marker='', markersize=1, color='r')
        ax2v.set_ylabel('%s' % (f.variables['vbar_northward'].name), color='r')  # , f.variables[var2plot[i]].units))
        ax2v.xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
        ax2v.set_ylim(-0.7, 0.7)
        #ylims = ax2v.get_ylim()
        ax.set_ylim(-0.7, 0.7)
        ax2v.xaxis.set_major_formatter(myFmt)
        ax.grid(True)
    elif var2plot[i] == 'bed_thickness':
        ax.stackplot(datetime_list,f.variables[var2plot[i]][:, 2, x, y],
                     f.variables[var2plot[i]][:, 1, x, y],
                     f.variables[var2plot[i]][:, 0, x, y])
        ax.set_ylabel('%s' % (f.variables[var2plot[i]].name))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
        ax.xaxis.set_major_formatter(myFmt)
        ax.grid(True)
    else:
        ax.plot_date(datetime_list, f.variables[var2plot[i]][:, x, y], xdate=True, linestyle='-', linewidth=1,
                     marker='', markersize=1)
        ax.set_ylabel('%s' % (f.variables[var2plot[i]].name))#, f.variables[var2plot[i]].units))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
        ax.xaxis.set_major_formatter(myFmt)
        ax.grid(True)
fig.suptitle('Point %fN %fE' % (f.variables['lat_rho'][x, y], f.variables['lon_rho'][x, y]))
#fig.axes[0].title('test')

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

#q = stick_plot(datetime_list, f.variables['ubar_eastward'][:,x,y], f.variables['vbar_northward'][:,x,y])