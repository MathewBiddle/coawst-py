'''
Biddle et al 2020 manuscript Figure 8:
Caption:
Top two panels: spatial distribution of depth averaged current speeds (m/s), from low (green) to high (pink/white),
with velocity vectors, to indicate direction, during the peak discharge event from Lee, on 9 September 2011 at
04:12, under vegetative (left panel) and non-vegetative (right panel) conditions. Bottom two panels: Spatial
distribution of the elevation difference (cm) between the final and initial bed thickness over the Lee time period,
under vegetative (left) and non-vegetative (right) conditions. Coloring indicates removal of elevation (blues) and
addition of elevation (reds) over the time period. A dashed line at zero is included in both panels to indicate the
delineation between removal and addition of elevation. A black dot is included in the map to indicate the location of
station FLT.

@author: Mathew Biddle
'''

import os
#os.environ["PROJ_LIB"] = "/Users/mbiddle/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import netCDF4
import coawstpy
import datetime

files = coawstpy.get_file_paths()
runs = ['veg','noveg']
#event = 'Lee'
#point_data = coawstpy.get_point_data(run)
#date = datetime.datetime(2011, 8, 2, 21, 0) # typical
#date = datetime.datetime(2011, 9, 9, 4, 12)  # Lee
#date = datetime.datetime(2011, 10, 21, 0, 0) # post-Lee
locs = coawstpy.get_point_locations()
flt_lon = locs.loc[locs['Site'] == 'FLT', 'lon'].values
flt_lat = locs.loc[locs['Site'] == 'FLT', 'lat'].values

times = coawstpy.get_time_periods()

for event in times:
    #init figure
    fig, ax = plt.subplots(ncols=1,nrows=2, figsize=(12, 10),sharey=True,sharex=True)

    i=0
    for run in runs:
        print(run)
        # runs_dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/'
        # if run == 'noveg':
        #     direct = runs_dir + 'Full_20110719T23_20111101_final_noveg'
        # elif run == 'veg':
        #     direct = runs_dir + 'Full_20110719T23_20111101_final'
        #
        # inputfile = direct + '/upper_ches_his.nc'
        # f = netCDF4.Dataset(inputfile, 'r')
        f = netCDF4.Dataset(files[run], 'r')
        lon = f.variables['lon_rho'][:]
        lat = f.variables['lat_rho'][:]
        ocean_time = f.variables['ocean_time'][:]
        plant_height = f.variables['mask_rho'][:]
        h = f.variables['h'][:]

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

        ## Bed elevation difference
        ## Do some date conversions ##
        datetime_list = []
        for sec in ocean_time:
            datetime_list.append(
                netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))

        start = coawstpy.nearest_ind(datetime_list, times[event][0])
        end = coawstpy.nearest_ind(datetime_list, times[event][1]) + 1

        # pick data to plot
        datetime_list = datetime_list[start:end]
        dvar = 'bed_thickness'
        data = f.variables[dvar][start:end, :, :, :]

        data_init = np.sum(data[0, :, :, :], axis=0)  # total from all layers
        data_init = np.ma.masked_where(plant_height != 5, data_init)

        data_final = np.sum(data[-1, :, :, :], axis=0)  # total from all layers
        data_final = np.ma.masked_where(plant_height != 5, data_final)

        data_diff = (data_final - data_init) * 100  # cm
        # set up figure

        lonm = np.ma.masked_where(plant_height != 5, lon)
        latm = np.ma.masked_where(plant_height != 5, lat)

        # set up map
        m = Basemap(llcrnrlon=lonm.min() - 0.01, llcrnrlat=latm.min() - 0.01, urcrnrlon=lonm.max() + 0.01,
                    urcrnrlat=latm.max() + 0.01,
                    resolution='i', projection='merc', ax=ax[i])
        # set up map
        # m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
        #    resolution='i', projection='merc', ax=ax[i])
        # m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax[i])
        # m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax[i])
        # m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

        # pcolor variable of interest
        cax = m.pcolormesh(lon, lat, data_diff, latlon=True,
                           vmin=-10, vmax=10, cmap='jet', ax=ax[i])
        #contour = m.contour(lon, lat, data_diff, 0,
        #                    colors='k', linestyles='dashed', linewidths=0.5, latlon=True, ax=ax[1,i])

        #add FLT point
        m.scatter(flt_lon, flt_lat, latlon=True, s=40, marker='.', color='k', edgecolors='k', linewidths=0.3, ax=ax[i])
        # add channel contour
        m.contour(lon, lat, h, [3], linewidths=1, linestyles=':', colors='k', latlon=True)

        i+=1
    #cbar.add_lines(contour)
    #m.colorbar(cax)
    #plt.suptitle("%s" % date)

    fig.subplots_adjust(right=0.8)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar_axb = fig.add_axes([0.125, 0.07, 0.675, 0.02])
    cbar = fig.colorbar(cax, cax=cbar_axb, orientation='horizontal')
    #cbar = fig.colorbar(cax, cax=cbar_ax)
    cbar.set_label('$\\Delta h_b$ ($cm$)')
    #cbar.add_lines(contour)

    plt.suptitle('%s' % event)

    # writedir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/Manuscript/figures/'
    image_name = '%s.png' % event
    # outfile = writedir+image_name
    # print("Saving image to %s" % outfile)
    plt.savefig(image_name, bbox_inches='tight', dpi=500)