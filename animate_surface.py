import os
os.environ["PROJ_LIB"] = "/Users/mbiddle/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import pandas as pd
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import datetime
import coawstpy

dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final_noveg'
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


# to subset time series
start_date = datetime.datetime(2011,10,20,19,0)
#end_date = datetime.datetime(2011,10,20,19,0)
#
start_idx = coawstpy.nearest_ind(datetime_list, start_date)
#end_idx = coawstpy.nearest_ind(datetime_list, end_date)

# to do entire time series
#start_idx = 1 # zero screws up cause velocities are zero
end_idx = -1

datetime_list = datetime_list[start_idx:end_idx]
mud_01 = mud_01[start_idx:end_idx, :, :]
ubar = ubar[start_idx:end_idx, :, :]
vbar = vbar[start_idx:end_idx, :, :]
ptsdf = ptsdf.iloc[start_idx:end_idx]
# set up figure
fig, ax = plt.subplots(figsize=(8, 6))

#set up map
m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
    resolution='i', projection='merc', ax=ax, epsg=3395)
#m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax)
#m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax)
#m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 300)

# pcolor variable of interest
#cax = m.pcolormesh(lon, lat, mud_01[0, :-1, :-1], latlon=True,
#                    vmin=0, vmax=0.6, cmap='gist_ncar', ax=ax)
terrainBig = cm.get_cmap('terrain',512)
new_terrain = ListedColormap(terrainBig(np.linspace(0.05, 1, 256)))
cax = m.pcolormesh(lon, lat, mud_01[0, :-1, :-1], latlon=True,
                    vmin=0, vmax=0.2, cmap=new_terrain, ax=ax)
cbar = fig.colorbar(cax, extend='max')
cbar.set_label('$\\widebar{SSC}_{f}$ (kg $m^{-3}$)')

#init wind vector
x, y = m(-75.979,39.38)
dax = m.quiver(x,y,ptsdf['X-Windv'].iloc[0],ptsdf['Y-Windv'].iloc[0])
#barb_increments = dict({'half':1,'full':5,'flag':15})
#dax = m.barbs(x,y,ptsdf['X-Windv'].iloc[0],ptsdf['Y-Windv'].iloc[0], barb_increments=barb_increments)

eax = m.quiver(lon[::3], lat[::3], ubar[0,::3], vbar[0,::3], latlon=True,
             pivot='tail', headaxislength=3, width=0.003, headwidth=4, scale_units='inches', scale=1.75)

#xr, yr = m(lon[23,1], lat[23,1])
#eax = m.quiver(xr, yr, ubar[0,23,1], vbar[0,23,1])#, pivot='middle', headaxislength=3, width=0.003)#, latlon=True)
#             pivot='middle', headaxislength=3, width=0.003)

fax = m.contour(lon, lat, plant_height, levels=[0], colors='k', linewidths=.5, latlon=True, corner_mask=True)

#xt, yt = m(lon[5,90], lat[5,90])
#gax = m.quiver(xt, yt, ubar[0,5,90], vbar[0,5,90])

#sys.exit()
def animate(i):
    # pcolor variable
    cax.set_array(mud_01[i, :-1, :-1].flatten())
    eax.set_UVC(ubar[i,::3], vbar[i,::3])
    #eax.set_UVC(ubar[i,23,1], vbar[i,23,1])
    #gax.set_UVC(ubar[i,5,90], vbar[i,5,90])
    # add wind vector
    dax.set_UVC(ptsdf['X-Windv'].iloc[i], ptsdf['Y-Windv'].iloc[i])
    plt.title('%s'% (datetime_list[i]))
    #plt.cla()



anim = animation.FuncAnimation(fig, animate, frames=len(datetime_list))
#plt.draw()
#plt.show()
print("Saving animation...")
outdir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/defense/'
anim.save(outdir+'mud_postLee-end_veg.mp4', fps=15, dpi=200)
#anim.save('out.gif', writer='imagemagick', fps=10)

print("Done.")
