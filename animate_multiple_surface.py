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

vegdir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'
novegdir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final_noveg'
## READ SWAN PTS data (for wind)
#ptsdir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
vegptsfile = vegdir+"/tripod_wave.pts"

vegptsdf = pd.read_fwf(vegptsfile, header=4)

vegptsdf.drop([0, 1], axis=0, inplace=True)
vegptsdf.rename(columns={'%       Time':'Time'}, inplace=True)
vegptsdf['Yp'] = vegptsdf['Yp            Hsig'].astype(str).str.split("    ", expand=True)[0].astype(float)
vegptsdf['Hsig'] = vegptsdf['Yp            Hsig'].astype(str).str.split("    ", expand=True)[1].astype(float)
vegptsdf.drop(columns=['Yp            Hsig'], inplace=True)

vegptsdf['Time'] = pd.to_datetime(vegptsdf['Time'], format='%Y%m%d.%H%M%S', utc=True)
vegptsdf['Hsig'] = vegptsdf['Hsig'].astype(float)
vegptsdf['X-Windv'] = vegptsdf['X-Windv'].astype(float)
vegptsdf['Y-Windv'] = vegptsdf['Y-Windv'].astype(float)

## Read COAWST data
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
veginputfile = vegdir+'/upper_ches_his.nc'
vegf = netCDF4.Dataset(veginputfile, 'r')
veglon = vegf.variables['lon_rho'][:][:]
veglat = vegf.variables['lat_rho'][:][:]
vegocean_time = vegf.variables['ocean_time'][:]
vegubar = vegf.variables['ubar_eastward'][:]
vegvbar = vegf.variables['vbar_northward'][:]
vegmud_01 = vegf.variables['mud_01'][:, :, :, :]
vegmud_01 = np.mean(vegmud_01,axis=1)

vegplant_height = vegf.variables['plant_height'][0, 0, :, :]
vegplant_height_orig = np.ma.masked_outside(vegplant_height,0.3,1)
## Do some date conversions ##
vegdatetime_list=[]
for sec in vegocean_time:
    vegdatetime_list.append(
        netCDF4.num2date(sec, units=vegf.variables['ocean_time'].units, calendar=vegf.variables['ocean_time'].calendar))


## pull in noveg data:
novegptsfile = novegdir+"/tripod_wave.pts"
novegptsdf = pd.read_fwf(novegptsfile, header=4)
novegptsdf.drop([0, 1], axis=0, inplace=True)
novegptsdf.rename(columns={'%       Time':'Time'}, inplace=True)
novegptsdf['Yp'] = novegptsdf['Yp            Hsig'].astype(str).str.split("    ", expand=True)[0].astype(float)
novegptsdf['Hsig'] = novegptsdf['Yp            Hsig'].astype(str).str.split("    ", expand=True)[1].astype(float)
novegptsdf.drop(columns=['Yp            Hsig'], inplace=True)
novegptsdf['Time'] = pd.to_datetime(novegptsdf['Time'], format='%Y%m%d.%H%M%S', utc=True)
novegptsdf['Hsig'] = novegptsdf['Hsig'].astype(float)
novegptsdf['X-Windv'] = novegptsdf['X-Windv'].astype(float)
novegptsdf['Y-Windv'] = novegptsdf['Y-Windv'].astype(float)

## Read COAWST data
noveginputfile = novegdir+'/upper_ches_his.nc'
novegf = netCDF4.Dataset(noveginputfile, 'r')
noveglon = novegf.variables['lon_rho'][:][:]
noveglat = novegf.variables['lat_rho'][:][:]
novegocean_time = novegf.variables['ocean_time'][:]
novegubar = novegf.variables['ubar_eastward'][:]
novegvbar = novegf.variables['vbar_northward'][:]
novegmud_01 = novegf.variables['mud_01'][:, :, :, :]
novegmud_01 = np.mean(novegmud_01,axis=1)

## Do some date conversions ##
novegdatetime_list=[]
for sec in novegocean_time:
    novegdatetime_list.append(
        netCDF4.num2date(sec, units=novegf.variables['ocean_time'].units, calendar=novegf.variables['ocean_time'].calendar))





# to subset time series
# start_date = datetime.datetime(2011,8,15)
# end_date = datetime.datetime(2011,9,15)
#
# start_idx = coawstpy.nearest_ind(datetime_list, start_date)
# end_idx = coawstpy.nearest_ind(datetime_list, end_date)

# to do entire time series
start_idx = 1 # zero screws up cause velocities are zero
end_idx = -1

novegdatetime_list = novegdatetime_list[start_idx:end_idx]
novegmud_01 = novegmud_01[start_idx:end_idx, :, :]
novegubar = novegubar[start_idx:end_idx, :, :]
novegvbar = novegvbar[start_idx:end_idx, :, :]
novegptsdf = novegptsdf.iloc[start_idx:end_idx]

vegdatetime_list = vegdatetime_list[start_idx:end_idx]
vegmud_01 = vegmud_01[start_idx:end_idx, :, :]
vegubar = vegubar[start_idx:end_idx, :, :]
vegvbar = vegvbar[start_idx:end_idx, :, :]
vegptsdf = vegptsdf.iloc[start_idx:end_idx]


# set up figure
fig, ax = plt.subplots(ncols=2,figsize=(8, 6))


#set up VEGETATION map
m = Basemap(llcrnrlon=veglon.min(), llcrnrlat=veglat.min(), urcrnrlon=veglon.max(), urcrnrlat=veglat.max(),
    resolution='i', projection='merc', ax=ax[0])

# pcolor variable of interest
vegcax = m.pcolormesh(veglon, veglat, vegmud_01[0, :-1, :-1], latlon=True,
                    vmin=0, vmax=0.6, cmap='terrain', ax=ax[0])
#cbar = fig.colorbar(vegcax)
#cbar.set_label('$\\widebar{SSC}_{f}$ (kg $m^{-3}$)')

#init wind vector
vegx, vegy = m(-75.979,39.38)
vegdax = m.quiver(vegx,vegy,vegptsdf['X-Windv'].iloc[0],vegptsdf['Y-Windv'].iloc[0],ax=ax[0])
#barb_increments = dict({'half':1,'full':5,'flag':15})
#dax = m.barbs(x,y,ptsdf['X-Windv'].iloc[0],ptsdf['Y-Windv'].iloc[0], barb_increments=barb_increments)

#eax = m.quiver(lon[::1], lat[::1], ubar[0,::1], vbar[0,::1], latlon=True,
#             pivot='middle', headaxislength=3, width=0.003)

vegxr, vegyr = m(veglon[23,1], veglat[23,1])
vegeax = m.quiver(vegxr, vegyr, vegubar[0,23,1], vegvbar[0,23,1],ax=ax[0])#, pivot='middle', headaxislength=3, width=0.003)#, latlon=True)
#             pivot='middle', headaxislength=3, width=0.003)

vegfax = m.contour(veglon, veglat, vegplant_height, levels=[0], colors='k', linewidths=.5, latlon=True, corner_mask=True)

vegxt, vegyt = m(veglon[5,90], veglat[5,90])
veggax = m.quiver(vegxt, vegyt, vegubar[0,5,90], vegvbar[0,5,90],ax=ax[0])


#set up NO VEGETATION map
m = Basemap(llcrnrlon=noveglon.min(), llcrnrlat=noveglat.min(), urcrnrlon=noveglon.max(), urcrnrlat=noveglat.max(),
    resolution='i', projection='merc', ax=ax[1])

# pcolor variable of interest
novegcax = m.pcolormesh(noveglon, noveglat, novegmud_01[0, :-1, :-1], latlon=True,
                    vmin=0, vmax=0.6, cmap='terrain', ax=ax[1])
#cbar = fig.colorbar(novegcax)
#cbar.set_label('$\\widebar{SSC}_{f}$ (kg $m^{-3}$)')

#init wind vector
novegx, novegy = m(-75.979,39.38)
novegdax = m.quiver(novegx,novegy,novegptsdf['X-Windv'].iloc[0],novegptsdf['Y-Windv'].iloc[0],ax=ax[1])
#barb_increments = dict({'half':1,'full':5,'flag':15})
#dax = m.barbs(x,y,ptsdf['X-Windv'].iloc[0],ptsdf['Y-Windv'].iloc[0], barb_increments=barb_increments)

#eax = m.quiver(lon[::1], lat[::1], ubar[0,::1], vbar[0,::1], latlon=True,
#             pivot='middle', headaxislength=3, width=0.003)

novegxr, novegyr = m(noveglon[23,1], noveglat[23,1])
novegeax = m.quiver(novegxr, novegyr, novegubar[0,23,1], novegvbar[0,23,1],ax=ax[1])#, pivot='middle', headaxislength=3, width=0.003)#, latlon=True)
#             pivot='middle', headaxislength=3, width=0.003)

#fax = m.contour(lon, lat, plant_height, levels=[0], colors='k', linewidths=.5, latlon=True, corner_mask=True)

novegxt, novegyt = m(noveglon[5,90], noveglat[5,90])
noveggax = m.quiver(novegxt, novegyt, novegubar[0,5,90], novegvbar[0,5,90],ax=ax[1])

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.815, 0.2, 0.03, 0.6])
cbar = fig.colorbar(novegcax, cax=cbar_ax, orientation='vertical', extend='max')
cbar.set_label('$\\widebar{SSC}_{f}$ (kg $m^{-3}$)')

#sys.exit()
def animate(i):
    # pcolor variable
    vegcax.set_array(vegmud_01[i, :-1, :-1].flatten())
    #eax.set_UVC(ubar[i,::1], vbar[i,::1])
    vegeax.set_UVC(vegubar[i,23,1], vegvbar[i,23,1])
    veggax.set_UVC(vegubar[i,5,90], vegvbar[i,5,90])
    # add wind vector
    vegdax.set_UVC(vegptsdf['X-Windv'].iloc[i], vegptsdf['Y-Windv'].iloc[i])
    ax[0].set_title('%s'% (vegdatetime_list[i]))

    # pcolor variable
    novegcax.set_array(novegmud_01[i, :-1, :-1].flatten())
    # eax.set_UVC(ubar[i,::1], vbar[i,::1])
    novegeax.set_UVC(novegubar[i, 23, 1], novegvbar[i, 23, 1])
    noveggax.set_UVC(novegubar[i, 5, 90], novegvbar[i, 5, 90])
    # add wind vector
    novegdax.set_UVC(novegptsdf['X-Windv'].iloc[i], novegptsdf['Y-Windv'].iloc[i])
    ax[1].set_title('%s' % (novegdatetime_list[i]))
    #plt.cla()




anim = animation.FuncAnimation(fig, animate, frames=len(vegdatetime_list))
#plt.draw()
#plt.show()
print("Saving animation...")
outdir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/defense/'
anim.save(outdir+'mud_full_combined.mp4', fps=15)
#anim.save('out.gif', writer='imagemagick', fps=10)

print("Done.")
