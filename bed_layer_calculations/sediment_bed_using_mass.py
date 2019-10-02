'''
Calculate change in bed using mass observations.
'''

import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import coawstpy
import numpy as np
import netCDF4
import scipy.integrate as integrate

## Read COAWST data
runs = ['veg','noveg']
event = 'post-Lee'
#point_data = coawstpy.get_point_data(run)
times = coawstpy.get_time_periods()
#locs = coawstpy.get_point_locations()

for run in runs:
    runs_dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/'
    if run == 'noveg':
        direct = runs_dir + 'Full_20110719T23_20111101_final_noveg'
    elif run == 'veg':
        direct = runs_dir + 'Full_20110719T23_20111101_final'
    inputfile = direct+'/upper_ches_his.nc'
    f = netCDF4.Dataset(inputfile, 'r')
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
    plant_height[:, 14:58] = 5 ## Southern Boundary
    plant_height[78:100, 41:59] = 0 # Remove Elk River
    plant_height[72:78, 55:59] = 0 # Remove Elk River
    ## Northern Boundary
    plant_height[48:64,0:14] = 5 # add mill creek/Furnace Bay
    # dealing with angled transect boundary here
    plant_height[30:50, 13] = 5
    plant_height[31:37, 12] = 5
    plant_height[32:35, 11] = 5

    #mud_mass=np.ma.masked_where(plant_height != 5, mud_mass)
    #sand_mass=np.ma.masked_where(plant_height != 5, sand_mass)

    hm = np.ma.masked_where(mask_rho == 0, h)  # initial masked water depth at rho points

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

    mud_mass_init = np.sum(mud_mass[start,:,:,:],axis=0) # total from all 3 bed layers
    mud_mass_final = np.sum(mud_mass[end,:,:,:],axis=0) # total from all 3 bed layers
    mud_mass_diff = mud_mass_final - mud_mass_init # kg/m2

    sand_mass_init = np.sum(sand_mass[start,:,:,:],axis=0) # total from all 3 bed layers
    sand_mass_final = np.sum(sand_mass[end,:,:,:],axis=0) # total from all 3 bed layers
    sand_mass_diff = sand_mass_final - sand_mass_init # kg/m2

    # Susquehanna River mouth
    trans_name = 'T1'
    t1_x = np.array([29,30,31,32])
    t1_y = np.array([13,12,11,10])

    for i in range(len(t1_x)):
        plant_height[t1_x[i], t1_y[i]] = 10

    # Turkey Point to Sandy Point
    trans_name = 'T2'
    t2_x = np.array(list(range(42,67)))  #
    t2_y = np.array([58]*len(t2_x))

    # Verify point location
    for i in range(len(t2_x)):
        #l = i/5
        plant_height[t2_x[i], t2_y[i]] = 10



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

    ## using sed mass
    mass_eroded = sand_mass_eroded + mud_mass_eroded
    mass_deposited = sand_mass_deposited + mud_mass_deposited
    print('Event: %s %s' % (event, run))
    print('mass eroded    = %e tons' % (np.sum(mass_eroded) / 1000))
    print('mass deposited = %e tons' % (np.sum(mass_deposited) / 1000))
    #sys.exit()

    ## Plotting
    #plt.figure()
    #plt.pcolor(f.variables['lon_rho'][:], f.variables['lat_rho'][:], plant_height)
    #plt.title('Transects')

    # set up figure
    fig, (ax) = plt.subplots(nrows=4, ncols=2, figsize=(8, 12))

    # set up map
    m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
        resolution='i', projection='merc', ax=ax[0,0])
    #m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax[0,0])
    #m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax[0,0])
    #m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

    # pcolor variable of interest
    cax0 = m.pcolormesh(lon, lat, mud_mass_diff_ma, latlon=True,
                        cmap='jet', ax=ax[0,0])
    cbar = fig.colorbar(cax0,ax=ax[0,0])
    cbar.set_label('mud mass diff [kg/m2]')
    #plt.title('%s mud mass evolution' % run)
    #plt.title("%s through %s"%(datetime_list[0],datetime_list[-1]))

    # set up figure
    #fig, ax = plt.subplots(figsize=(8, 6))

    # set up map
    m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
        resolution='i', projection='merc', ax=ax[0,1])
    #m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax[0,1])
    #m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax[0,1])
    #m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

    # pcolor variable of interest
    cax1 = m.pcolormesh(lon, lat, sand_mass_diff_ma, latlon=True,
                        cmap='jet', ax=ax[0,1])
    cbar = fig.colorbar(cax1,ax=ax[0,1])
    cbar.set_label('sand mass diff [kg/m2]')
    #plt.title('%s sand mass evolution' % run)
    #plt.title("%s through %s"%(datetime_list[0],datetime_list[-1]))

    # set up map
    m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
        resolution='i', projection='merc', ax=ax[1,0])
    #m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax[1,0])
    #m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax[1,0])
    #m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

    # pcolor variable of interest
    cax1 = m.pcolormesh(lon, lat, mud_mass_eroded, latlon=True,
                        cmap='jet', ax=ax[1,0])
    cbar = fig.colorbar(cax1,ax=ax[1,0])
    cbar.set_label('mud mass eroded [kg]')


    # set up map
    m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
        resolution='i', projection='merc', ax=ax[1,1])
    #m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax[1,1])
    #m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax[1,1])
    #m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

    # pcolor variable of interest
    cax1 = m.pcolormesh(lon, lat, sand_mass_eroded, latlon=True,
                        cmap='jet', ax=ax[1,1])
    cbar = fig.colorbar(cax1,ax=ax[1,1])
    cbar.set_label('sand mass eroded [kg]')

    # set up map
    m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
        resolution='i', projection='merc', ax=ax[2,0])
    #m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax[2,0])
    #m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax[2,0])
    #m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

    # pcolor variable of interest
    cax1 = m.pcolormesh(lon, lat, mud_mass_deposited, latlon=True,
                        cmap='jet', ax=ax[2,0])
    cbar = fig.colorbar(cax1,ax=ax[2,0])
    cbar.set_label('mud mass deposited [kg]')


    # set up map
    m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
        resolution='i', projection='merc', ax=ax[2,1])
    #m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax[2,1])
    #m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax[2,1])
    #m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

    # pcolor variable of interest
    cax1 = m.pcolormesh(lon, lat, sand_mass_deposited, latlon=True,
                        cmap='jet', ax=ax[2,1])
    cbar = fig.colorbar(cax1,ax=ax[2,1])
    cbar.set_label('sand mass deposited [kg]')
    # set up figure
    #fig, ax = plt.subplots(figsize=(8, 6))

    # set up map
    m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
        resolution='i', projection='merc', ax=ax[3,0])
    #m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax[3,0])
    #m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax[3,0])
    #m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

    # pcolor variable of interest
    cax2 = m.pcolormesh(lon, lat, mass_eroded, latlon=True,
                        cmap='jet', ax=ax[3,0])
    cbar = fig.colorbar(cax2,ax=ax[3,0])
    cbar.set_label('Total bed mass eroded [kg]')
    #plt.title('%s total mass eroded = %e tons' % (run, np.sum(mass_eroded) / 1000))
    #plt.title("%s through %s"%(datetime_list[0],datetime_list[-1]))

    # set up figure
    #fig, ax = plt.subplots(figsize=(8, 6))

    # set up map
    m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
        resolution='i', projection='merc', ax=ax[3,1])
    #m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax[3,1])
    #m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax[3,1])
    #m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

    # pcolor variable of interest
    cax3 = m.pcolormesh(lon, lat, mass_deposited, latlon=True,
                        cmap='jet', ax=ax[3,1])
    cbar = fig.colorbar(cax3,ax=ax[3,1])
    cbar.set_label('Total bed mass deposited [kg]')
    #plt.title('%s total mass deposited = %e tons' % (run, np.sum(mass_deposited) / 1000))
    #plt.title("%s through %s"%(datetime_list[0],datetime_list[-1]))
    plt.subplots_adjust(wspace=0.005)
    plt.suptitle('%s %s mass evolution' % (event, run))