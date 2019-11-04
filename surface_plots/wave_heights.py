import os
os.environ["PROJ_LIB"] = "/Users/mbiddle/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import netCDF4
import coawstpy
import datetime

runs = ['veg','noveg']
event = 'post-Lee'
point_data = coawstpy.get_point_data('veg')

#date = datetime.datetime(2011, 8, 28, 12, 59, 57) # Irene
date = datetime.datetime(2011, 10, 20, 18, 59, 57) # Post-Lee
#locs = coawstpy.get_point_locations()
u = point_data['CBIBS'].loc[date,'X-Windv']
v = point_data['CBIBS'].loc[date,'Y-Windv']
#coawstpy.stick_plot(ptsdf['Time'],ptsdf['X-Windv'],ptsdf['Y-Windv'], ax=ax[0],scale=200)

i=0
fig, ax = plt.subplots(ncols=2,figsize=(20, 10),sharey=True,sharex=True)
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
    plant_height.harden_mask()  # makes the mask immutable
    plant_height[:, 10:58] = 5  ## Southern Boundary
    plant_height[78:100, 41:59] = 0  # Remove Elk River
    plant_height[72:78, 55:59] = 0  # Remove Elk River
    ## Northern Boundary
    plant_height[48:64, 0:14] = 5  # add mill creek/Furnace Bay
    # dealing with angled transect boundary here
    # plant_height[30:50, 13] = 5
    # plant_height[31:37, 12] = 5
    # plant_height[32:35, 11] = 5


    ## Do some date conversions ##
    datetime_list=[]
    for sec in ocean_time:
        datetime_list.append(
            netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))

    time = coawstpy.nearest_ind(datetime_list, date)
#    end = coawstpy.nearest_ind(datetime_list, times[event][1]) + 1

    # pick data to plot
    datetime_list = datetime_list[time]
    #dvar = 'bed_thickness'
    Hwave = f.variables['Hwave'][time, :, :]

    Hwavem = np.ma.masked_where(plant_height != 5,Hwave)

    #curr_mag = np.sqrt((ubarm**2)+(vbrm**2))
    print(run)
    print('Maximum wave: %f' % Hwavem.max())
#    data_diff = (data_final-data_init)*100 # cm
    # set up figure

    lonm = np.ma.masked_where(plant_height != 5, lon)
    latm = np.ma.masked_where(plant_height != 5, lat)

    #set up map
    m = Basemap(llcrnrlon=lonm.min()-0.01, llcrnrlat=latm.min()-0.01, urcrnrlon=lonm.max()+0.01, urcrnrlat=latm.max()+0.01,
        resolution='i', projection='merc', ax=ax[i])

    # pcolor variable of interest
    caxm = m.pcolormesh(lon, lat, Hwavem, latlon=True, ax=ax[i],vmin=0,vmax=0.5)
    #m.quiver(lon, lat, ubarm, vbarm, latlon=True, ax=ax[i])
#                        vmin=-0.02,vmax=0.02,cmap='jet', ax=ax[i])

    #contour = m.contour(lon, lat, data_diff, 0,
    #                    colors='k', linestyles='dashed', linewidths=0.5, latlon=True, ax=ax[i])

    ax[i].set_title("%s" % run)
    i+=1

x,y=m(lonm.min()+.01,latm.min()+.01)
m.quiver(x, y, u, v, scale=200, ax=ax[0])
m.quiver(x, y, u, v, scale=200, ax=ax[1])
#coawstpy.stick_plot(date,u,v, ax=caxm)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.125, 0.17, 0.675, 0.03])
cbar = fig.colorbar(caxm, cax=cbar_ax, orientation='horizontal')
cbar.set_label('significant wave height %s [m]' % date)
#cbar.add_lines(contour)
#m.colorbar(cax)
#plt.suptitle("%s" % date)


#writedir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/figures/significant_wave_height_maps/'
#image_name = '%s_map.png' % datetime_list
#outfile = writedir+image_name
#print("Saving image to %s" % outfile)
#plt.savefig(outfile, bbox_inches='tight', dpi=500)