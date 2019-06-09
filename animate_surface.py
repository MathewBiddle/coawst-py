import pandas as pd
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import datetime
import time


#dir='/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110702_20111101'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110906_20110926'
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
#dir = '/Volumes/Documents/COAWST_34_UPPER_CHES_FULL'
#inputfile = dir+'/upper_ches_his.nc'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110714201800_20111031231800'
inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')
lon = f.variables['lon_rho'][:][:]
lat = f.variables['lat_rho'][:][:]
ocean_time = f.variables['ocean_time'][:]
#lat = f.variables['lat_rho'][:][:]
#lon = f.variables['lon_rho'][:][:]

## Do some date conversions ##
#epoch_date = '%s %s'%(f.variables['ocean_time'].units.split(' ')[-2], f.variables['ocean_time'].units.split(' ')[-1])
#dt_obj = datetime.datetime.strptime(epoch_date, '%Y-%m-%d %H:%M:%S')
#time_diff = time.mktime(dt_obj.timetuple())-time.mktime(datetime.datetime(1970, 1, 1, 0, 0, 0).timetuple())
datetime_list=[]
for sec in ocean_time:
    #ts = sec/(12*3600)
    datetime_list.append(
        netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))
    #datetime_list.append(datetime.datetime.fromtimestamp(sec+time_diff))

dvar = 'mud_01'
mud_01 = f.variables[dvar][:, 0, :, :]

fig, ax = plt.subplots(figsize=(8, 6))
ax.set(xlim=(lon.min(), lon.max()), ylim=(lat.min(), lat.max()))

cax = ax.pcolormesh(lon, lat, mud_01[0, :-1, :-1],
                    vmin=0, vmax=1, cmap='BrBG_r')
fig.colorbar(cax)

#Writer = animation.writers['gif']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

def animate(i):
    cax.set_array(mud_01[i, :-1, :-1].flatten())
    plt.title('%s %s'% (dvar, datetime_list[i]))


anim = animation.FuncAnimation(fig, animate, frames=len(datetime_list))
#plt.draw()
#plt.show()

anim.save('out.gif', writer='imagemagick', fps=10)
