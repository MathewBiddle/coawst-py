import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import numpy as np
import netCDF4
import pandas as pd
import coawstpy

# Cindy's locations
locs = pd.DataFrame(columns=['Site', 'lat', 'lon','comment'])

locations = [
    ['1', 39.527, -76.061, 'Russ and Palinkas 2018'],
    ['2', 39.533, -76.061, ''],
    ['3', 39.515, -76.051, 'Middle of bed'],
    ['4', 39.505, -76.039, ''],
    ['5', 39.497, -76.036, ''],
    ['Lee7', 39.414, -76.079, ''],
    ['Lee6', 39.38, -76.088, ''],
    ['Lee5', 39.346, -76.197, ''],
    ['Lee2.5', 39.197, -76.311, ''],
    ['Lee2', 39.135, -76.328, ''],
    ['Lee0', 39.061, -76.328, ''],
    ['LeeS2', 38.757, -76.473, ''],
    ['CBIBS', 39.5396, -76.0741, 'CBIBS Susquehanna Flats'],
    ['Tripod', 39.4931, -76.0341, 'Larry tripod site']
    ]

i = 0
for location in locations:
    locs.loc[i] = location
    i += 1

z_pt = -1 # 0=bottom 4=surface (-1)

#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final_noveg'
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'
if dir.split("_")[-1] == 'noveg':
    run = "noveg"
else:
    run = "veg"
inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')
ocean_time = f.variables['ocean_time'][:]
lat = f.variables['lat_rho'][:]
lon = f.variables['lon_rho'][:]
## Do some date conversions ##
datetime_list=[]
for sec in ocean_time:
    if sec == 0.0:
        datetime_list.append(
            netCDF4.num2date(sec + 0.0000000000000001, units=f.variables['ocean_time'].units,
                             calendar=f.variables['ocean_time'].calendar))
    else:
        datetime_list.append(
            netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))


river_frc = dir+'/river_frc.nc'
f_river = netCDF4.Dataset(river_frc, 'r')
river_time = f_river.variables['river_time'][:]
river_transport = f_river.variables['river_transport'][:, 0]
river_datetime_list=[]
for sec in river_time:
    ts = sec/(12*3600)
    river_datetime_list.append(
        netCDF4.num2date(sec, units=f_river.variables['river_time'].units, calendar='standard'))

ptsfile = dir+"/tripod_wave.pts"
ptsdf = pd.read_fwf(ptsfile, header=4)
ptsdf.drop([0,1],axis=0,inplace=True)
ptsdf.rename(columns={'%       Time':'Time'},inplace=True)
ptsdf['Yp'] = ptsdf['Yp            Hsig'].astype(str).str.split("    ",expand=True)[0].astype(float)
ptsdf['Hsig'] = ptsdf['Yp            Hsig'].astype(str).str.split("    ",expand=True)[1].astype(float)
ptsdf.drop(columns=['Yp            Hsig'],inplace=True)
ptsdf['Time']=pd.to_datetime(ptsdf['Time'],format='%Y%m%d.%H%M%S',utc=True)
ptsdf['Hsig']=ptsdf['Hsig'].astype(float)
ptsdf['X-Windv']=ptsdf['X-Windv'].astype(float)
ptsdf['Y-Windv']=ptsdf['Y-Windv'].astype(float)


# collect data for site of choice
sites = ['1','2','3','4','5','Lee7','Lee6','CBIBS','Tripod']
for site in sites:
    print("plotting site: %s" % site)
    lat_pt, lon_pt = locs.loc[locs['Site'] == site, ['lat', 'lon']].values[0]
    try:
        lat_pt
        lon_pt
    except NameError:
        print("Using indexes x, y = (%i, %i)" % (x, y))
    else:
        print("Using geo-coords lat, lon = (%f, %f)" % (lat_pt, lon_pt))
        x = np.abs(lon[:, 1]-lon_pt).argmin()
        y = np.abs(lat[1, :]-lat_pt).argmin()

    # Verify point location
    plant_height = f.variables['plant_height'][0, 0, :, :]
    plant_height = np.ma.masked_greater(plant_height, 1)
    plant_height[x, y] = 1
    plt.figure(1)
    plt.pcolor(lon, lat, plant_height)
    plt.title('Site %s' % site)

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
            ax2v = ax.twinx()
            ax2v.plot_date(datetime_list, f.variables['Pwave_top'][:, x, y],
                           xdate=True, linestyle='-', linewidth=0.5,
                           marker='', markersize=1, color='r')
            ax2v.set_ylabel('Period [s]')
            ax2v.tick_params(axis='y', colors='red')
            ax2v.yaxis.label.set_color('red')
            ax2v.set_ylim(0, 6)
            ax2v.set_yticks([0, 3, 6])
            # ax2v.grid(True)
            ax2v.set_xlim(xlim)
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

    #outfile = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/figures/timeseries/veg/site_%s_timeseries_lowres.png' % site
    #print("Saving to %s" % outfile)
    #plt.savefig(outfile, bbox_inches='tight')#, dpi=1000)

