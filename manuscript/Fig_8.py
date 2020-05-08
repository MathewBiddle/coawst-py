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
os.environ["PROJ_LIB"] = "/Users/mbiddle/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import netCDF4
import coawstpy
import datetime

runs = ['veg','noveg']
event = 'Lee'
#point_data = coawstpy.get_point_data(run)
#date = datetime.datetime(2011, 8, 2, 21, 0) # typical
date = datetime.datetime(2011, 9, 9, 4, 12)  # Lee
#date = datetime.datetime(2011, 10, 21, 0, 0) # post-Lee
locs = coawstpy.get_point_locations()
flt_lon = locs.loc[locs['Site'] == 'FLT', 'lon'].values
flt_lat = locs.loc[locs['Site'] == 'FLT', 'lat'].values

times = coawstpy.get_time_periods()

#init figure
fig, ax = plt.subplots(ncols=2,nrows=2, figsize=(12, 10),sharey=True,sharex=True)

i=0
for run in runs:
    print(run)
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
    ubar = f.variables['ubar_eastward'][time, :, :]
    vbar = f.variables['vbar_northward'][time, :, :]

    ubarm = np.ma.masked_where(plant_height != 5,ubar)
    vbarm = np.ma.masked_where(plant_height != 5,vbar)

    curr_mag = np.sqrt((ubarm**2)+(vbarm**2))
    print('Maximum current: %f' % curr_mag.max())
#    data_diff = (data_final-data_init)*100 # cm
    # set up figure

    lonm = np.ma.masked_where(plant_height != 5, lon)
    latm = np.ma.masked_where(plant_height != 5, lat)

    #set up map
    m = Basemap(llcrnrlon=lonm.min()-0.01, llcrnrlat=latm.min()-0.01, urcrnrlon=lonm.max()+0.01, urcrnrlat=latm.max()+0.01,
        resolution='i', projection='merc', ax=ax[0,i])

#     # pcolor variable of interest
    caxm = m.pcolormesh(lon, lat, curr_mag, latlon=True, ax=ax[0,i],vmin=0,vmax=2.5, cmap='gist_ncar')
#     #m.quiver(lon, lat, ubarm, vbarm, latlon=True, ax=ax[i])
# #                        vmin=-0.02,vmax=0.02,cmap='jet', ax=ax[i])
    m.quiver(lon[::3], lat[::3], ubarm[::3], vbarm[::3], latlon=True, ax=ax[0,i],
             pivot='tail', headaxislength=4, headwidth=4, width=0.002)

    ## testing unit vector
    # caxm = m.quiver(lon[::3], lat[::3], ubarm[::3], vbarm[::3], np.sqrt(ubarm[::3]**2+vbarm[::3]**2),
    #                 latlon=True, ax=ax[0, i], clim=(0, 2.5), pivot='tail', cmap='gist_ncar',
    #                 angles='xy', scale_units='xy', scale=1,
    #                 width=0.01, headwidth=1, headlength=10, minshaft=5)

    # add FLT point
    m.scatter(flt_lon, flt_lat, latlon=True, s=40, marker='.', color='k', edgecolors='k', linewidths=0.3, ax=ax[0, i])

    #contour = m.contour(lon, lat, data_diff, 0,
    #                    colors='k', linestyles='dashed', linewidths=0.5, latlon=True, ax=ax[i])

    ax[0,i].set_title("%s" % run)

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
                resolution='i', projection='merc', ax=ax[1,i])
    # set up map
    # m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
    #    resolution='i', projection='merc', ax=ax[i])
    # m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax[i])
    # m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax[i])
    # m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

    # pcolor variable of interest
    cax = m.pcolormesh(lon, lat, data_diff, latlon=True,
                       vmin=-10, vmax=10, cmap='jet', ax=ax[1,i])
    contour = m.contour(lon, lat, data_diff, 0,
                        colors='k', linestyles='dashed', linewidths=0.5, latlon=True, ax=ax[1,i])

    #add FLT point
    m.scatter(flt_lon, flt_lat, latlon=True, s=40, marker='.', color='k', edgecolors='k', linewidths=0.3, ax=ax[1, i])

    i+=1
fig.subplots_adjust(right=0.8)
cbar_axt = fig.add_axes([0.125, 0.51, 0.675, 0.02])
cbar = fig.colorbar(caxm, cax=cbar_axt, orientation='horizontal', extend='max')
cbar.set_label('$\\widebar{\\upsilon}$ ($m$ $s^{-1}$) on %s' % date)
#cbar.add_lines(contour)
#m.colorbar(cax)
#plt.suptitle("%s" % date)

fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar_axb = fig.add_axes([0.125, 0.07, 0.675, 0.02])
cbar = fig.colorbar(cax, cax=cbar_axb, orientation='horizontal')
#cbar = fig.colorbar(cax, cax=cbar_ax)
cbar.set_label('$\\Delta h_b$ ($cm$)')
cbar.add_lines(contour)

writedir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/Manuscript/figures/'
image_name = 'Fig_8.png'
outfile = writedir+image_name
print("Saving image to %s" % outfile)
plt.savefig(outfile, bbox_inches='tight', dpi=500)