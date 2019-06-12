import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import pandas as pd
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

## READ SWAN PTS data (for wind)
ptsdir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
ptsfile = ptsdir+"/tripod_wave.pts"

ptsdf = pd.read_fwf(ptsfile, header=4)

ptsdf.drop([0, 1], axis=0, inplace=True)
ptsdf.rename(columns={'%       Time':'Time'}, inplace=True)
ptsdf['Yp'] = ptsdf['Yp            Hsig'].astype(str).str.split("    ", expand=True)[0].astype(float)
ptsdf['Hsig'] = ptsdf['Yp            Hsig'].astype(str).str.split("    ", expand=True)[1].astype(float)
ptsdf.drop(columns=['Yp            Hsig'], inplace=True)

ptsdf['Time'] = pd.to_datetime(ptsdf['Time'], format='%Y%m%d.%H%M%S', utc=True)
ptsdf['Hsig'] = ptsdf['Hsig'].astype(float)
ptsdf['X-Windv'] = ptsdf['X-Windv'].astype(float)
ptsdf['Y-Windv'] = ptsdf['Y-Windv'].astype(float)

## Read COAWST data
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')
lon = f.variables['lon_rho'][:][:]
lat = f.variables['lat_rho'][:][:]
ocean_time = f.variables['ocean_time'][:]

## Do some date conversions ##
datetime_list=[]
for sec in ocean_time:
    datetime_list.append(
        netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))

# pick data to plot
dvar = 'sand_01'
sand_01 = f.variables[dvar][:, 0, :, :]

# set up figure
fig, ax = plt.subplots(figsize=(8, 6))

#set up map
m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
    resolution='i', projection='merc', ax=ax)
m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax)
m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax)
#m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

# pcolor variable of interest
cax = m.pcolormesh(lon, lat, sand_01[0, :-1, :-1], latlon=True,
                    vmin=0, vmax=1, cmap='jet', ax=ax)
cbar = fig.colorbar(cax)
cbar.set_label('surface sand SSC [kg/m3]')

#init wind vector
x, y = m(-75.979,39.38)
dax = m.quiver(x,y,ptsdf['X-Windv'].iloc[0],ptsdf['Y-Windv'].iloc[0])

def animate(i):
    # pcolor variable
    cax.set_array(sand_01[i, :-1, :-1].flatten())
    # add wind vector
    dax.set_UVC(ptsdf['X-Windv'].iloc[i*2], ptsdf['Y-Windv'].iloc[i*2])
    plt.title('%s %s'% (dvar, datetime_list[i]))
    #plt.cla()



anim = animation.FuncAnimation(fig, animate, frames=len(datetime_list))
#plt.draw()
#plt.show()
print("Saving animation...")
anim.save('sand_01.mp4', fps=10)
#anim.save('out.gif', writer='imagemagick', fps=10)

print("Done.")
