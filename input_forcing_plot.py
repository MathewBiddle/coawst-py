import netCDF4
import pandas as pd
import matplotlib.pyplot as plt
import coawstpy
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
import numpy as np


dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'


## river data
river_frc = dir+'/river_frc.nc'
f_river = netCDF4.Dataset(river_frc, 'r')
river_time = f_river.variables['river_time'][:]
river_transport = f_river.variables['river_transport'][:, 0]
river_datetime_list=[]
for sec in river_time:
    river_datetime_list.append(
        netCDF4.num2date(sec, units=f_river.variables['river_time'].units, calendar='standard'))


## wind data
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

## tide data
bry_file = dir + '/upper_ches_bry.nc'
f_bry = netCDF4.Dataset(bry_file, 'r')
bry_time = f_bry.variables['zeta_time'][:]
bry_zeta = f_bry.variables['zeta_south'][:, 50]
bry_datetime_list=[]
for sec in bry_time:
    bry_datetime_list.append(netCDF4.num2date(sec, units=f_bry.variables['zeta_time'].units, calendar='standard'))

## plotting
fig, (ax) = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(12, 8))
fig.subplots_adjust(hspace=0.1)
myFmt = DateFormatter("%m/%d")
dayint=10
xlim= (ptsdf['Time'].min(), ptsdf['Time'].max())

ax[0].plot_date(river_datetime_list,
                (river_transport + (0.2 * river_transport)) * f_river.variables['river_transport'].shape[1],
                xdate=True, linestyle='-', linewidth=0.5,
                marker='', markersize=1)
ax[0].xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
ax[0].xaxis.set_major_formatter(myFmt)
ax[0].set_xlim(xlim)
ax[0].set_ylabel('Discharge (m3/s)')

coawstpy.stick_plot(ptsdf['Time'],ptsdf['X-Windv'],ptsdf['Y-Windv'], ax=ax[1])
ax[1].xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
ax[1].xaxis.set_major_formatter(myFmt)
ax[1].set_xlim(xlim)
ax[1].set_ylabel('Wind')

ax[2].plot_date(ptsdf['Time'], np.sqrt(ptsdf['X-Windv']**2 + ptsdf['Y-Windv']**2),
                xdate=True, linestyle='-', linewidth=0.5, marker='', markersize=1)
ax[2].xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
ax[2].xaxis.set_major_formatter(myFmt)
ax[2].set_xlim(xlim)
ax[2].set_ylabel('Wind Speed [m/s]')


ax[3].plot_date(bry_datetime_list,bry_zeta, xdate=True, linestyle='-', linewidth=0.5,
                     marker='', markersize=1)
ax[3].xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
ax[3].xaxis.set_major_formatter(myFmt)
ax[3].set_xlim(xlim)
ax[3].set_ylabel('Water surface [m]')



