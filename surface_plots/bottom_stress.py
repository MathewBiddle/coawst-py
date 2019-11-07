import os
os.environ["PROJ_LIB"] = "/Users/mbiddle/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import netCDF4
import coawstpy
import datetime

runs = ['veg']#,'noveg']
event = 'typical'
#point_data = coawstpy.get_point_data(run)
#date = datetime.datetime(2011, 8, 28, 13, 00) # Irene
date = datetime.datetime(2011, 10, 20, 19, 00) # Post-Lee
#locs = coawstpy.get_point_locations()
dates = [datetime.datetime(2011, 8, 28, 13, 00),datetime.datetime(2011, 10, 20, 19, 00)]
i=0
fig, ax = plt.subplots(ncols=3,nrows=2,figsize=(10, 10),sharey=True,sharex=True)
for date in dates:
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
        bstressmax_mag = f.variables['bstrcwmax'][time, :, :] # max wave and current btm stress

        bustrc = f.variables['bustrc'][time, :, :] # current u stress
        bvstrc = f.variables['bvstrc'][time, :, :] # current v stress
        bstrc_mag = np.sqrt(bustrc**2 + bvstrc**2) # magnitude

        bustrw = f.variables['bustrw'][time, :, :] # wind induced u stress
        bvstrw = f.variables['bvstrw'][time, :, :]  # wind induced v stress
        bstrw_mag = np.sqrt(bustrw ** 2 + bvstrw ** 2) # magnitude

        # apply mask
        bstressmax_magm = np.ma.masked_where(plant_height != 5, bstressmax_mag)
        bstrc_magm = np.ma.masked_where(plant_height != 5, bstrc_mag)
        bstrw_magm = np.ma.masked_where(plant_height != 5, bstrw_mag)

        #curr_mag = np.sqrt((ubarm**2)+(vbarm**2))
        #print('Maximum current: %f' % curr_mag.max())
        #    data_diff = (data_final-data_init)*100 # cm
        # set up figure

        lonm = np.ma.masked_where(plant_height != 5, lon)
        latm = np.ma.masked_where(plant_height != 5, lat)

        #set up map
        m = Basemap(llcrnrlon=lon.min()-0.01, llcrnrlat=lat.min()-0.01, urcrnrlon=lon.max()+0.01,
                    urcrnrlat=lat.max()+0.01, resolution='i', projection='merc', ax=ax[i,0])
        caxm = m.pcolormesh(lon, lat, bstressmax_mag, latlon=True, ax=ax[i,0], vmin=0, vmax=0.2, cmap='jet')
        #ax[i,0].set_title("%s bstress tot" % run)

        m = Basemap(llcrnrlon=lon.min() - 0.01, llcrnrlat=lat.min() - 0.01, urcrnrlon=lon.max() + 0.01,
                    urcrnrlat=lat.max() + 0.01, resolution='i', projection='merc', ax=ax[i,1])
        caxm = m.pcolormesh(lon, lat, bstrc_mag, latlon=True, ax=ax[i, 1], vmin=0, vmax=0.2, cmap='jet')
        #ax[i,1].set_title("%s current bstress" % run)

        m = Basemap(llcrnrlon=lon.min() - 0.01, llcrnrlat=lat.min() - 0.01, urcrnrlon=lon.max() + 0.01,
                    urcrnrlat=lat.max() + 0.01, resolution='i', projection='merc', ax=ax[i,2])
        caxm = m.pcolormesh(lon, lat, bstrw_mag, latlon=True, ax=ax[i, 2], vmin=0, vmax=0.2, cmap='jet')
        #ax[i,2].set_title("%s wave bstress" % run)
        # pcolor variable of interest

        #m.quiver(lon, lat, ubarm, vbarm, latlon=True, ax=ax[i])
        #                        vmin=-0.02,vmax=0.02,cmap='jet', ax=ax[i])

        #contour = m.contour(lon, lat, data_diff, 0,
        #                    colors='k', linestyles='dashed', linewidths=0.5, latlon=True, ax=ax[i])


        i+=1
#fig.subplots_adjust(right=0.8)
fig.subplots_adjust(hspace=0.01,wspace=0.01)
cbar_ax = fig.add_axes([0.125, 0.09, 0.775, 0.02])
cbar = fig.colorbar(caxm, cax=cbar_ax, orientation='horizontal', extend='max')
cbar.set_label('bottom stress magnitude (N/m$^2$)')
ax[0,0].set_title("total")
ax[0,1].set_title("current")
ax[0,2].set_title("wave")
ax[0,0].set_ylabel("Irene %s" % dates[0])
ax[1,0].set_ylabel("Post-Lee %s" % dates[1])
#cbar.add_lines(contour)
#m.colorbar(cax)
#plt.suptitle("%s" % date)


#writedir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/figures/bottom_stress_maps/'
#image_name = 'bottom_stress_compare_map.png'
#outfile = writedir+image_name
#print("Saving image to %s" % outfile)
#plt.savefig(outfile, bbox_inches='tight', dpi=500)