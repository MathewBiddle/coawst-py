'''
Biddle et al 2020 manuscript Figure 8:
Caption:
Spatial distribution of the mud mass difference (top two panels) and sand mass difference (kg/m2) (bottom two
panels) between the final and initial mass in the bed layer sum over the time period associated with the Lee
event, under vegetative (left) and non-vegetative (right) conditions. Coloring indicates removal of mass (blues and
greens) and addition of mass (reds) over the time period. A dashed line is included in all panels to indicate the
delineation between removal and addition of mass. A black dot is included in the map to indicate the location of
station FLT.

@author: Mathew Biddle
'''

import os
#os.environ["PROJ_LIB"] = "/Users/mbiddle/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import coawstpy
import numpy as np
import netCDF4
import scipy.integrate as integrate

## Read COAWST data
runs = ['veg','noveg']
event = 'Lee'
transects = coawstpy.get_transect_indexes()
times = coawstpy.get_time_periods()
files = coawstpy.get_file_paths()
#locs = coawstpy.get_point_locations()
locs = coawstpy.get_point_locations()
flt_lon = locs.loc[locs['Site'] == 'FLT', 'lon'].values
flt_lat = locs.loc[locs['Site'] == 'FLT', 'lat'].values

fig, (ax) = plt.subplots(nrows=2, ncols=2,sharex=True, sharey=True, figsize=(12, 10))
i=0
for run in runs:
    # runs_dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/'
    # if run == 'noveg':
    #     direct = runs_dir + 'Full_20110719T23_20111101_final_noveg'
    # elif run == 'veg':
    #     direct = runs_dir + 'Full_20110719T23_20111101_final'
    # inputfile = direct+'/upper_ches_his.nc'
    f = netCDF4.Dataset(files[run], 'r')
    lon = f.variables['lon_rho'][:][:]
    lat = f.variables['lat_rho'][:][:]
    h = f.variables['h'][:] # at rho points
    ocean_time = f.variables['ocean_time'][:]
    bed_thick = f.variables['bed_thickness'][:] # m
    sand_mass = f.variables['sandmass_01'][:] # kg/m2
    mud_mass = f.variables['mudmass_01'][:] # kg/m2 already accounts for cell height.
    mask_rho = f.variables['mask_rho'][:]
    Srho = f.variables['Srho'][0] # kg/m3 2650 kg/m^3 sediment grain density, same for both
    pm = f.variables['pm'][:] #XI --> cell width in x dir. east-west 1/m
    pn = f.variables['pn'][:] #ETA --> cell width in y dir. north-south Want to use this for Surface Area Calcs 1/m
    plant_height = f.variables['mask_rho'][:]

    # build Flats cell locations:
    # This does not include the transect cells themselves. Just all cells between T1 and T2, excluding South River.
    # 5 = good points
    plant_height = np.ma.masked_less(plant_height, 1)
    plant_height.harden_mask() # makes the mask immutable
    plant_height[:, 10:58] = 5 ## Southern Boundary
    plant_height[78:100, 41:59] = 0 # Remove Elk River
    plant_height[72:78, 55:59] = 0 # Remove Elk River
    ## Northern Boundary
    plant_height[48:64,0:14] = 5 # add mill creek/Furnace Bay
    # dealing with angled transect boundary here
    #plant_height[30:50, 13] = 5
    #plant_height[31:37, 12] = 5
    #plant_height[32:35, 11] = 5

    #mud_mass=np.ma.masked_where(plant_height != 5, mud_mass)
    #sand_mass=np.ma.masked_where(plant_height != 5, sand_mass)

    hm = np.ma.masked_where(mask_rho == 0, h)  # initial masked water depth at rho points

    lonm = np.ma.masked_where(plant_height != 5, lon)
    latm = np.ma.masked_where(plant_height != 5, lat)

    SA = (1/pm) * (1/pn) # m2
    SAm = np.ma.masked_where(mask_rho == 0, SA) # mask land
    cell_vol_init = hm * SAm # total initial volume of water in valid cells. assuming all cells are square!

    # checking the volume against what ROMS gives
    print('Initial domain volumes:  TotVolume = 6.1046316405e+08 m3\n   Calculated Initial Total Volume = %e m3' % np.sum(cell_vol_init))

    ## Do some date conversions ##
    datetime_list=[]
    for sec in ocean_time:
        datetime_list.append(
            netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))

    start = coawstpy.nearest_ind(datetime_list, times[event][0])
    end = coawstpy.nearest_ind(datetime_list, times[event][1]) + 1
    # Calculate the delta bed mass
    # time indexes: 1202 - 1466 for event
    # see https://www.myroms.org/forum/viewtopic.php?f=20&t=4447

    #mud_mass_init = np.sum(mud_mass[start,0,:,:],axis=0) # total from all 3 bed layers [t,z,x,y] z = 0 surface
    #mud_mass_final = np.sum(mud_mass[end,0,:,:],axis=0) # total from all 3 bed layers
    mud_mass_init = mud_mass[start, 0, :, :]  # [t,z,x,y] z = 0 surface
    mud_mass_final = mud_mass[end,0,:,:]
    mud_mass_diff = mud_mass_final - mud_mass_init # kg/m2

    #sand_mass_init = np.sum(sand_mass[start,0,:,:],axis=0) # total from all 3 bed layers
    #sand_mass_final = np.sum(sand_mass[end,0,:,:],axis=0) # total from all 3 bed layers
    sand_mass_init = sand_mass[start, 0, :, :]  # surface
    sand_mass_final = sand_mass[end,0,:,:] # surface
    sand_mass_diff = sand_mass_final - sand_mass_init # kg/m2

    # Susquehanna River mouth
    plant_height[transects['T1']['x'], transects['T1']['y']] = 10
    # Turkey Point to Sandy Point
    plant_height[transects['T1']['x'], transects['T1']['y']] = 10

    ## Plotting
    # apply the mask
    # plant_height = 5 is where the region of interest is.
    # so apply a mask to everything not 5 to the bed thick matrix

    sand_mass_diff_ma = np.ma.masked_where(plant_height != 5, sand_mass_diff)
    sand_mass_deposition = np.ma.masked_less(sand_mass_diff_ma, 0) # deposited sand kg/m2
    sand_mass_erosion = np.ma.masked_greater(sand_mass_diff_ma, 0) # eroded sand kg/m2
    sand_mass_deposited = sand_mass_deposition * SAm # kg/m^2 * m^2 = kg
    sand_mass_eroded = sand_mass_erosion * SAm # kg/m^2 * m^2 = kg

    mud_mass_diff_ma = np.ma.masked_where(plant_height != 5, mud_mass_diff)
    mud_mass_deposition = np.ma.masked_less(mud_mass_diff_ma, 0) # deposited sand kg/m2
    mud_mass_erosion = np.ma.masked_greater(mud_mass_diff_ma, 0) # eroded sand kg/m2
    mud_mass_deposited = mud_mass_deposition * SAm # kg/m^2 * m^2 = kg
    mud_mass_eroded = mud_mass_erosion * SAm # kg/m^3 * m^2 = kg


    ## Plotting
    # set up figure


    # set up map
    m = Basemap(llcrnrlon=lonm.min()-0.01, llcrnrlat=latm.min()-0.01, urcrnrlon=lonm.max()+0.01, urcrnrlat=latm.max()+0.01,
        resolution='i', projection='merc', ax=ax[0,i])

    # pcolor variable of interest
    cax0 = m.pcolormesh(lon, lat, mud_mass_diff_ma, latlon=True,
                        cmap='jet', ax=ax[0,i],
                        vmin=-1, vmax=1)
                        # post-Lee and Irene vmin=-1,vmax=0.6)
                        # Lee vmin=-4,vmax=4)
                        # typical vmin=-0.35,vmax=0.15)
    # contour0 = m.contour(lon, lat, mud_mass_diff_ma, 0,
    #                     colors='k', linestyles='dashed', linewidths=0.5, latlon=True, ax=ax[0,i])

    # add channel contour
    m.contour(lon, lat, h, [3], linewidths=1, linestyles=':', colors='k', latlon=True, ax=ax[0,i])
    # cbar0 = m.colorbar(cax0, ax=ax[0, i], location='right', shrink=0.6)
    # cbar0.add_lines(contour0)
    # cbar0.set_label('mud mass diff [kg/m2]')
    ax[0,i].set_title('%s' % (run))
    m.scatter(flt_lon, flt_lat, latlon=True, s=40, marker='.', color='k', edgecolors='k', linewidths=0.3, ax=ax[0, i])

    #i+=1
    # set up map
    m = Basemap(llcrnrlon=lonm.min()-0.01, llcrnrlat=latm.min()-0.01, urcrnrlon=lonm.max()+0.01, urcrnrlat=latm.max()+0.01,
        resolution='i', projection='merc', ax=ax[1,i])

    # pcolor variable of interest
    cax1 = m.pcolormesh(lon, lat, sand_mass_diff_ma, latlon=True,
                        cmap='jet', ax=ax[1,i],
                        vmin=-100, vmax=100)
                        # post-Lee and Irene vmin=-6, vmax=6)
                        # Lee vmin=-100,vmax=100)
                        # typical vmin=-1.5,vmax=0.5)
    # contour1 = m.contour(lon, lat, sand_mass_diff_ma, 0,
    #                     colors='k', linestyles='dashed', linewidths=0.5, latlon=True, ax=ax[1,i])

    # add channel contour
    m.contour(lon, lat, h, [3], linewidths=1, linestyles=':', colors='k', latlon=True, ax=ax[1,i])

    # cbar1 = m.colorbar(cax1, ax=ax[1, i], location='right', shrink=0.6)
    # cbar1.add_lines(contour1)
    # cbar1.set_label('sand mass diff [kg/m$^{2}$]')
    ax[1,i].set_title('%s' % (run))
    m.scatter(flt_lon, flt_lat, latlon=True, s=40, marker='.', color='k', edgecolors='k', linewidths=0.3, ax=ax[1, i])

    #plt.suptitle('%s %s mass evolution' % (event, run))
    i+=1
    ## print out some values
    mass_eroded = sand_mass_eroded + mud_mass_eroded
    mass_deposited = sand_mass_deposited + mud_mass_deposited
    print('\n\nEvent: %s %s' % (event, run))
    print('total mass eroded    = %e tons' % (np.sum(mass_eroded) / 1000))
    print('total mass deposited = %e tons' % (np.sum(mass_deposited) / 1000))
    print('\ntotal sand eroded    = %e tons' % (np.sum(sand_mass_eroded) / 1000))
    print('total sand deposited = %e tons' % (np.sum(sand_mass_deposited) / 1000))
    print('Sand change = %f tons' % ((np.sum(sand_mass_eroded)+np.sum(sand_mass_deposited))/ 1000))
    print('\ntotal mud eroded    = %e tons' % (np.sum(mud_mass_eroded) / 1000))
    print('total mud deposited = %e tons' % (np.sum(mud_mass_deposited) / 1000))
    print('Mud change = %f tons' % ((np.sum(mud_mass_eroded)+np.sum(mud_mass_deposited))/ 1000))

    #cbar = m.colorbar(cax0, ax=ax[i,0], location='bottom')
fig.subplots_adjust(right=0.8)
cbar_ax0 = fig.add_axes([0.125, 0.51, 0.675, 0.02])
cbar0 = fig.colorbar(cax0, cax=cbar_ax0, orientation='horizontal', extend='both')
cbar0.set_label('$\\Delta$ $m_f$ ($kg$ $m^{-2}$)')
# cbar0.add_lines(contour0)

cbar_ax1 = fig.add_axes([0.125, 0.07, 0.675, 0.02])
cbar1 = fig.colorbar(cax1, cax=cbar_ax1, orientation='horizontal', extend='both')
cbar1.set_label('$\\Delta$ $m_c$ ($kg$ $m^{-2}$)')
# cbar1.add_lines(contour1)

# writedir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/Manuscript/figures/'
image_name = 'Fig_9.png'
# outfile = writedir+image_name
# print("Saving image to %s" % outfile)
#plt.savefig(image_name, bbox_inches='tight', dpi=500)
