'''
Biddle et al 2020 manuscript Figure 6:
Caption:
Top two panels are the spatial distribution of significant wave heights (m), smaller (blue) to larger (yellow),
during peak winds in Irene, on 29 August 2011 at 13:00, under vegetative (left) and non-vegetative (right) conditions.
On each panel an arrow is included to indicate wind direction and speed at the specified time. The scale bar for
these panels ranges from 0 (purple) to 0.5 (yellow) m.

Bottom two panels are the spatial distribution of the mud mass difference (delta m_f kg/m2) between the final and
initial mass in the bed layer sum over the time period associated with the Irene event, under vegetative
(left) and non-vegetative (right) conditions. Coloring indicates removal of mass (blues and greens) and
addition of mass (reds) over the time period. A dashed line is included in both panels to indicate the
delineation between removal and addition of mass. The scale bar for these panels ranges from -1 (blue) to 0.6 (red) kg m^-2.

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
#f = netCDF4.Dataset(inputfile, 'r')

runs = ['veg','noveg']
event = 'Irene'
point_data = coawstpy.get_point_data('veg')
times = coawstpy.get_time_periods()
transects = coawstpy.get_transect_indexes()
date = datetime.datetime(2011, 8, 28, 12, 59, 57) # Irene
#date = datetime.datetime(2011, 10, 20, 18, 59, 57) # Post-Lee
#locs = coawstpy.get_point_locations()
u = point_data['CBIBS'].loc[date,'X-Windv']
v = point_data['CBIBS'].loc[date,'Y-Windv']
#coawstpy.stick_plot(ptsdf['Time'],ptsdf['X-Windv'],ptsdf['Y-Windv'], ax=ax[0],scale=200)

i=0
fig, ax = plt.subplots(ncols=2,nrows=2, figsize=(12, 10),sharey=True,sharex=True)

for run in runs:
    #runs_dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/'
    #if run == 'noveg':
    #    direct = runs_dir + 'Full_20110719T23_20111101_final_noveg'
    #elif run == 'veg':
    #    direct = runs_dir + 'Full_20110719T23_20111101_final'

    #inputfile = direct + '/upper_ches_his.nc'
    f = netCDF4.Dataset(files[run], 'r')
    lon = f.variables['lon_rho'][:]
    lat = f.variables['lat_rho'][:]
    h = f.variables['h'][:]
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
        resolution='i', projection='merc', ax=ax[0,i])

    # pcolor variable of interest
    caxm = m.pcolormesh(lon, lat, Hwavem, latlon=True, ax=ax[0,i],vmin=0,vmax=0.5,cmap='jet')

    # add 3m isobath
    m.contour(lon, lat, h, [3], linewidths=1, linestyles=':', colors='k', latlon=True, ax=ax[0, i])
    #m.quiver(lon, lat, ubarm, vbarm, latlon=True, ax=ax[i])
#                        vmin=-0.02,vmax=0.02,cmap='jet', ax=ax[i])

    #contour = m.contour(lon, lat, data_diff, 0,
    #                    colors='k', linestyles='dashed', linewidths=0.5, latlon=True, ax=ax[i])

    ax[0,i].set_title("%s" % run)


    ## mud mass difference
    mud_mass = f.variables['mudmass_01'][:]
    datetime_list = []
    for sec in ocean_time:
        datetime_list.append(
            netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))

    start = coawstpy.nearest_ind(datetime_list, times[event][0])
    end = coawstpy.nearest_ind(datetime_list, times[event][1]) + 1
    # Calculate the delta bed mass
    # time indexes: 1202 - 1466 for event
    # see https://www.myroms.org/forum/viewtopic.php?f=20&t=4447

    mud_mass_init = np.sum(mud_mass[start, :, :, :], axis=0)  # total from all 3 bed layers
    mud_mass_final = np.sum(mud_mass[end, :, :, :], axis=0)  # total from all 3 bed layers
    mud_mass_diff = mud_mass_final - mud_mass_init  # kg/m2

    # Susquehanna River mouth
    plant_height[transects['T1']['x'], transects['T1']['y']] = 10
    # Turkey Point to Sandy Point
    plant_height[transects['T1']['x'], transects['T1']['y']] = 10

    ## Plotting
    # apply the mask
    # plant_height = 5 is where the region of interest is.
    # so apply a mask to everything not 5 to the bed thick matrix
    mud_mass_diff_ma = np.ma.masked_where(plant_height != 5, mud_mass_diff)

    ## Plotting
    # set up figure

    # set up map
    m = Basemap(llcrnrlon=lonm.min() - 0.01, llcrnrlat=latm.min() - 0.01, urcrnrlon=lonm.max() + 0.01,
                urcrnrlat=latm.max() + 0.01,
                resolution='i', projection='merc', ax=ax[1, i])

    # pcolor variable of interest
    cax0 = m.pcolormesh(lon, lat, mud_mass_diff_ma, latlon=True,
                        cmap='jet', ax=ax[1, i],
                        vmin=-1, vmax=1)
    # post-Lee and Irene vmin=-1,vmax=0.6)
    # Lee vmin=-4,vmax=4)
    # typical vmin=-0.35,vmax=0.15)
    # contour0 = m.contour(lon, lat, mud_mass_diff_ma, 0,
    #                      colors='k', linestyles='dashed', linewidths=0.5, latlon=True, ax=ax[1, i])
    # cbar0 = m.colorbar(cax0, ax=ax[0, i], location='right', shrink=0.6)
    # cbar0.add_lines(contour0)
    # cbar0.set_label('mud mass diff [kg/m2]')
    # add 3m isobath
    m.contour(lon, lat, h, [3], linewidths=1, linestyles=':', colors='k', latlon=True, ax=ax[1, i])

    ax[1, i].set_title('%s' % (run))

    i+=1

x,y=m(lonm.min()+.01,latm.min()+.01)
m.quiver(x, y, u, v, scale=200, ax=ax[0,0])
m.quiver(x, y, u, v, scale=200, ax=ax[0,1])
#coawstpy.stick_plot(date,u,v, ax=caxm)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.125, 0.51, 0.675, 0.02])
cbar = fig.colorbar(caxm, cax=cbar_ax, orientation='horizontal')
cbar.set_label('$H_s$ ($m$) on %s' % date) # change to H_s at [date]
#cbar.add_lines(contour)
#m.colorbar(cax)
#plt.suptitle("%s" % date)

cbar_axb = fig.add_axes([0.125, 0.09, 0.675, 0.02])
cbarb = fig.colorbar(cax0, cax=cbar_axb, orientation='horizontal')
#cbar = fig.colorbar(cax, cax=cbar_ax)
cbarb.set_label('$\\Delta m_{f}$ ($kg$ $m^{-2}$)') # change in mud mass
#cbarb.add_lines(contour0)
#m.colorbar(cax)
#plt.suptitle("%s %s through %s" % (event, datetime_list[0],datetime_list[-1]))

## maybe look into the horizonatal spacing between boxes, something like hspace is the option



#writedir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/Manuscript/figures/'
image_name = 'Fig_6_revision.png'
#image_name = '%s_map.png' % datetime_list
#outfile = writedir+image_name
#print("Saving image to %s" % outfile)
plt.savefig(image_name, bbox_inches='tight', dpi=500)