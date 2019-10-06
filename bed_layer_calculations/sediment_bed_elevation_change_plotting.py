import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import netCDF4
import coawstpy

runs = ['veg','noveg']
event = 'typical'
#point_data = coawstpy.get_point_data(run)
times = coawstpy.get_time_periods()
locs = coawstpy.get_point_locations()

i=0
fig, ax = plt.subplots(ncols=2,figsize=(8, 6),sharey=True,sharex=True)
for run in runs:
    runs_dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/'
    if run == 'noveg':
        direct = runs_dir + 'Full_20110719T23_20111101_final_noveg'
    elif run == 'veg':
        direct = runs_dir + 'Full_20110719T23_20111101_final'

    inputfile = direct + '/upper_ches_his.nc'
    f = netCDF4.Dataset(inputfile, 'r')
    lon = f.variables['lon_rho'][:]
    lat = f.variables['lat_rho'][:]
    ocean_time = f.variables['ocean_time'][:]
    plant_height = f.variables['mask_rho'][:]

    # build Flats cell locations:
    # This does not include the transect cells themselves. Just all cells between T1 and T2, excluding South River.
    # 5 = good points
    plant_height = np.ma.masked_less(plant_height, 1)
    plant_height.harden_mask() # makes the mask immutable
    plant_height[:, 14:58] = 5 ## Southern Boundary
    plant_height[78:100, 41:59] = 0 # Remove Elk River
    plant_height[72:78, 55:59] = 0 # Remove Elk River
    ## Northern Boundary
    plant_height[48:64,0:14] = 5 # add mill creek/Furnace Bay
    # dealing with angled transect boundary here
    plant_height[30:50, 13] = 5
    plant_height[31:37, 12] = 5
    plant_height[32:35, 11] = 5


    ## Do some date conversions ##
    datetime_list=[]
    for sec in ocean_time:
        datetime_list.append(
            netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))

    start = coawstpy.nearest_ind(datetime_list, times[event][0])
    end = coawstpy.nearest_ind(datetime_list, times[event][1]) + 1

    # pick data to plot
    datetime_list = datetime_list[start:end]
    dvar = 'bed_thickness'
    data = f.variables[dvar][start:end, :, :, :]

    data_init=np.sum(data[0,:,:,:],axis=0) # total from all layers
    data_init = np.ma.masked_where(plant_height != 5,data_init)

    data_final=np.sum(data[-1,:,:,:],axis=0) # total from all layers
    data_final = np.ma.masked_where(plant_height != 5,data_final)

    data_diff = (data_final-data_init)*100 # cm
    # set up figure


    #set up map
    m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
        resolution='i', projection='merc', ax=ax[i])
    #m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax[i])
    #m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax[i])
    #m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

    # pcolor variable of interest
    cax = m.pcolormesh(lon, lat, data_diff, latlon=True,
                        vmin=-0.05,vmax=0.05,cmap='jet', ax=ax[i])
    contour = m.contour(lon, lat, data_diff, 0,
                        colors='k', linestyles='dashed', linewidths=0.5, latlon=True, ax=ax[i])
    #cbar = fig.colorbar(cax)
    #cbar.set_label('Bed evolution [cm]')
    #m.scatter(-76.079,39.414,marker='o',color='k',alpha=0.4,edgecolors='k',linewidths=0.3,latlon=True)
    #m.scatter(-76.088,39.380,marker='o',color='k',alpha=0.4,edgecolors='k',linewidths=0.3,latlon=True)
    ax[i].set_title("%s" % run)
    i+=1
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(cax, cax=cbar_ax)
cbar.set_label('Bed evolution [cm]')
cbar.add_lines(contour)
#m.colorbar(cax)
plt.suptitle("%s %s through %s" % (event, datetime_list[0],datetime_list[-1]))