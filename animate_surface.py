import os
os.environ["PROJ_LIB"] = "/Users/mbiddle/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import pandas as pd
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import datetime
import coawstpy

dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'
## READ SWAN PTS data (for wind)
#ptsdir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
ptsfile = dir+"/tripod_wave.pts"

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
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')
lon = f.variables['lon_rho'][:][:]
lat = f.variables['lat_rho'][:][:]
ocean_time = f.variables['ocean_time'][:]
ubar = f.variables['ubar_eastward'][:]
vbar = f.variables['vbar_northward'][:]
mud_01 = f.variables['mud_01'][:, :, :, :]
mud_01 = np.mean(mud_01,axis=1)

plant_height = f.variables['plant_height'][0, 0, :, :]
plant_height_orig = np.ma.masked_outside(plant_height,0.3,1)
## Do some date conversions ##
datetime_list=[]
for sec in ocean_time:
    datetime_list.append(
        netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))

start_date = datetime.datetime(2011,8,15)
end_date = datetime.datetime(2011,9,15)


start_idx = coawstpy.nearest_ind(datetime_list, start_date)
end_idx = coawstpy.nearest_ind(datetime_list, end_date)

datetime_list = datetime_list[start_idx:end_idx]
mud_01 = mud_01[start_idx:end_idx, :, :]
ubar = ubar[start_idx:end_idx, :, :]
vbar = vbar[start_idx:end_idx, :, :]
ptsdf = ptsdf.iloc[start_idx:end_idx]
# set up figure
fig, ax = plt.subplots(figsize=(8, 6))

#set up map
m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
    resolution='i', projection='merc', ax=ax)
#m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax)
#m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax)
#m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

# pcolor variable of interest
cax = m.pcolormesh(lon, lat, mud_01[0, :-1, :-1], latlon=True,
                    vmin=0, vmax=0.6, cmap='terrain', ax=ax)
cbar = fig.colorbar(cax)
cbar.set_label('$\\widebar{SSC}_{f}$ (kg $m^{-3}$)')

#init wind vector
x, y = m(-75.979,39.38)
dax = m.quiver(x,y,ptsdf['X-Windv'].iloc[0],ptsdf['Y-Windv'].iloc[0])

eax = m.quiver(lon[::2], lat[::2], ubar[0,::2], vbar[0,::2], latlon=True,
             pivot='tail', headaxislength=3, width=0.003)

fax = m.contour(lon, lat, plant_height, levels=[0], colors='k', linewidths=.5, latlon=True, corner_mask=True)


#sys.exit()
def animate(i):
    # pcolor variable
    cax.set_array(mud_01[i, :-1, :-1].flatten())
    eax.set_UVC(ubar[i,::2], vbar[i,::2])
    # add wind vector
    dax.set_UVC(ptsdf['X-Windv'].iloc[i], ptsdf['Y-Windv'].iloc[i])
    plt.title('%s'% (datetime_list[i]))
    #plt.cla()



anim = animation.FuncAnimation(fig, animate, frames=len(datetime_list))
#plt.draw()
#plt.show()
print("Saving animation...")
outdir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Presentations/defense/'
anim.save(outdir+'mud_irene_lee.mp4', fps=20)
#anim.save('out.gif', writer='imagemagick', fps=10)

print("Done.")
