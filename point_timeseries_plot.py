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
import datetime
import time
import netCDF4
import pandas as pd
import coawstpy
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
site = '3'
lat_pt, lon_pt = locs.loc[locs['Site'] == site, ['lat', 'lon']].values[0]

z_pt = 4 # 0=bottom 4=surface
## Using point location idx
#x = 56   # 57
#y = 27      #

#dir='/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110702_20111101'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110906_20110926'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
dir = '/Volumes/Documents/COAWST_34_UPPER_CHES_FULL'
#inputfile = dir+'/upper_ches_his.nc'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110714201800_20111031231800'
inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')

#roms_in_dir = '/Volumes/Documents/ROMS_INPUT_FILES'
river_frc = dir+'/../ROMS_INPUT_FILES/river_frc.nc'

f_river = netCDF4.Dataset(river_frc, 'r')
river_time = f_river.variables['river_time'][:]
river_transport = f_river.variables['river_transport'][:, 0]

## Do some date conversions ##
#epoch_date = '%s %s'%(f_river.variables['river_time'].units.split(' ')[-2], f_river.variables['river_time'].units.split(' ')[-1])
#dt_obj = datetime.datetime.strptime(epoch_date, '%Y-%m-%d %H:%M:%S')
#time_diff = time.mktime(dt_obj.timetuple())-time.mktime(datetime.datetime(1970, 1, 1, 0, 0, 0).timetuple())
river_datetime_list=[]
for sec in river_time:
    ts = sec/(12*3600)
    river_datetime_list.append(
        netCDF4.num2date(sec, units=f_river.variables['river_time'].units, calendar='standard'))
    #river_datetime_list.append(datetime.datetime.fromtimestamp(sec+time_diff))

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

## Do some date conversions ##
#epoch_date = '%s %s'%(f.variables['ocean_time'].units.split(' ')[-2], f.variables['ocean_time'].units.split(' ')[-1])
#dt_obj = datetime.datetime.strptime(epoch_date, '%Y-%m-%d %H:%M:%S')
#time_diff = time.mktime(dt_obj.timetuple())-time.mktime(datetime.datetime(1970, 1, 1, 0, 0, 0).timetuple())
datetime_list=[]
for sec in ocean_time:
    if sec == 0.0:
        datetime_list.append(
            netCDF4.num2date(sec + 0.0000000000000001, units=f.variables['ocean_time'].units,
                             calendar=f.variables['ocean_time'].calendar))
    else:
        datetime_list.append(
            netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))
    #datetime_list.append(datetime.datetime.fromtimestamp(sec+time_diff))

# Verify point location
plant_height = f.variables['plant_height'][0, 0, :, :]
plant_height = np.ma.masked_greater(plant_height, 1)
plant_height[x, y] = 1
plt.figure(1)
plt.pcolor(f.variables['lon_rho'][:][:], f.variables['lat_rho'][:][:], plant_height)
plt.title('Site %s' % site)
#outfile = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/figures/%s_map.png' % site
#plt.savefig(outfile, bbox_inches='tight', dpi = 1000)

ptsdir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
ptsfile = ptsdir+"/tripod_wave.pts"
#file='/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/SWAN_20130705_20130715_FRICTION_NOVEG_30SEC_JANSSEN_pt4+Bathy/tripod_wave.pts'


ptsdf = pd.read_fwf(ptsfile, header=4)#,widths=[18,16,16,16,16,16,16,16,16])#,skiprows=range(0,5))

ptsdf.drop([0,1],axis=0,inplace=True)
ptsdf.rename(columns={'%       Time':'Time'},inplace=True)
ptsdf['Yp'] = ptsdf['Yp            Hsig'].astype(str).str.split("    ",expand=True)[0].astype(float)
ptsdf['Hsig'] = ptsdf['Yp            Hsig'].astype(str).str.split("    ",expand=True)[1].astype(float)
ptsdf.drop(columns=['Yp            Hsig'],inplace=True)

ptsdf['Time']=pd.to_datetime(ptsdf['Time'],format='%Y%m%d.%H%M%S',utc=True)
ptsdf['Hsig']=ptsdf['Hsig'].astype(float)
ptsdf['X-Windv']=ptsdf['X-Windv'].astype(float)
ptsdf['Y-Windv']=ptsdf['Y-Windv'].astype(float)


# plot variables as time series
myFmt = DateFormatter("%m/%d")

var2plot = ['depth','Hwave','velocity','mud+sand','bed_thickness','river_transport','wind']

fig, (ax) = plt.subplots(nrows=len(var2plot), ncols=1, sharex=True, figsize=(12, 8))
fig.subplots_adjust(hspace=0.05)
dayint = 10
for i, ax in enumerate(fig.axes):
    if (var2plot[i] == 'rho') or (var2plot[i] == 'salt') or (var2plot[i] == 'tke') or (var2plot[i] == 'mud_01') or \
            (var2plot[i] == 'sand_01') or (var2plot[i] == 'temp'):
        ax.plot_date(datetime_list, f.variables[var2plot[i]][:, z_pt, x, y], xdate=True, linestyle='-', linewidth=0.5,
                     marker='', markersize=1)
        ax.set_ylabel('%s_%s' % (f.variables[var2plot[i]].name,(z_pt+1)))  # , f.variables[var2plot[i]].units))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
        ax.xaxis.set_major_formatter(myFmt)
        ax.grid(True)
        xlim = ax.get_xlim()

    elif var2plot[i] == 'mud+sand':
        ax.plot_date(datetime_list, f.variables['mud_01'][:, z_pt, x, y], label='mud', xdate=True, linestyle='-',linewidth=0.5,
                     marker='', markersize=1, color='b')
        ax.plot_date(datetime_list, f.variables['sand_01'][:, z_pt, x, y], label='sand', xdate=True, linestyle='-', linewidth=0.5,
                     marker='', markersize=1, color='r')
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
        ax.xaxis.set_major_formatter(myFmt)

        ax.set_ylim(0, 0.5)
#        ax.yaxis.tick_right()
#        ax.yaxis.set_label_position("right")
        ax.set_ylabel('SSC [kg/m3]')
        ax.grid(True)
        ax.legend(loc="upper left")
        xlim = ax.get_xlim()

    elif var2plot[i] == 'velocity':
        ax.plot_date(datetime_list,f.variables['ubar_eastward'][:, x, y], label='u',
                     xdate=True, linestyle='-', linewidth=0.5, marker='', markersize=1, color='b')
        ax.plot_date(datetime_list,f.variables['vbar_northward'][:, x, y], label='v',
                       xdate=True, linestyle='-', linewidth=0.5, marker='', markersize=1, color='r')

        ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_formatter(myFmt)
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.set_ylabel('Current [m/s]')
        ax.set_ylim(-0.6, 0.6)
        ax.xaxis.set_major_formatter(myFmt)
        ax.grid(True)
        ax.legend(loc="upper left")
        xlim = ax.get_xlim()

    elif var2plot[i] == 'bed_thickness':
        ax.stackplot(datetime_list,f.variables[var2plot[i]][:, 2, x, y],
                     f.variables[var2plot[i]][:, 1, x, y],
                     f.variables[var2plot[i]][:, 0, x, y])
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.set_ylabel('bed [m]')
        ax.set_ylim(0.8, 1.1)
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
        ax.xaxis.set_major_formatter(myFmt)
        ax.grid(True)
        xlim = ax.get_xlim()

    elif var2plot[i] == 'river_transport':
        # river trnasport had 20% reduction and destributed across cells. Need to back calculate
        ax.plot_date(river_datetime_list, (river_transport + (0.2 * river_transport)) * f_river.variables['river_transport'].shape[1],
                     xdate=True, linestyle='-', linewidth=0.5,
                     marker='', markersize=1)
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
        ax.xaxis.set_major_formatter(myFmt)
        ax.grid(True)
        ax.set_xlim(xlim)
        #ax.yaxis.tick_right()
        #ax.yaxis.set_label_position("right")
        ax.set_ylabel('Q [m3/s]')

    elif var2plot[i] == 'depth':
        ax.plot_date(datetime_list, f.variables['zeta'][:, x, y]+f.variables['h'][x,y], xdate=True, linestyle='-', linewidth=0.5,
                     marker='', markersize=1)
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.set_ylabel('Depth [m]')
        ax.set_ylim(0,6)
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
        ax.xaxis.set_major_formatter(myFmt)
        ax.grid(True)
        xlim = ax.get_xlim()

    elif var2plot[i] == 'wind':
        coawstpy.stick_plot(ptsdf['Time'], ptsdf['X-Windv'], ptsdf['Y-Windv'], ax=ax)
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
        ax.xaxis.set_major_formatter(myFmt)
        ax.xaxis.grid(True)
        #ax.get_yaxis().set_ticks([])
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.set_ylabel('Wind')
        xlim = ax.get_xlim()

    elif var2plot[i] == 'Hwave':
        ax.plot_date(datetime_list, f.variables[var2plot[i]][:, x, y], xdate=True, linestyle='-', linewidth=0.5,
                     marker='', markersize=1)
        #ax.yaxis.tick_right()
        #ax.yaxis.set_label_position("right")
        ax.set_ylabel('Sig. Wave H. [m]')
        ax.set_ylim(0,0.4)
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
        ax.xaxis.set_major_formatter(myFmt)
        ax.grid(True)
        xlim = ax.get_xlim()
    else:
        ax.plot_date(datetime_list, f.variables[var2plot[i]][:, x, y], xdate=True, linestyle='-', linewidth=0.5,
                     marker='', markersize=1)
        ax.set_ylabel('%s' % (f.variables[var2plot[i]].name))#, f.variables[var2plot[i]].units))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
        ax.xaxis.set_major_formatter(myFmt)
        ax.grid(True)
        xlim = ax.get_xlim()

fig.suptitle('Site %s @ %fN %fE' % (site, f.variables['lat_rho'][x, y], f.variables['lon_rho'][x, y]))
#fig.axes[0].title('test')

#outfile = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/figures/timeseries/site_%s_timeseries_lowres.png' % site
#print("Saving to %s" % outfile)
#plt.savefig(outfile, bbox_inches='tight')#, dpi=1000)

