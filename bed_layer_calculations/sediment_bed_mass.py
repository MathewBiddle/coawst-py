import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import netCDF4
import scipy.integrate as integrate

## Read COAWST data
# from coupling.out
# Initial domain volumes:  TotVolume =  6.1046316405E+08 m3
#                         MinCellVol =  1.0721942507E+03 m3
#                         MaxCellVol =  1.0719327656E+05 m3
#                            Max/Min =  9.9975612154E+01

dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final_noveg'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'
if dir.split("_")[-1] == 'noveg':
    run = "noveg"
else:
    run = "veg"
inputfile = dir+'/upper_ches_his.nc'
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
#sys.exit()
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

# Calculate the delta bed thickness
# time indexes: 1202 - 1466 for event
bed_thick_init = np.sum(bed_thick[0,:,:,:],axis=0) # total from all 3 bed layers
bed_thick_final = np.sum(bed_thick[-1,:,:,:],axis=0) # total from all 3 bed layers
bed_thick_diff = bed_thick_final - bed_thick_init # meters

# TODO Ask Cindy/Larry if sum converts from kg/m2 -> kg??
# see https://www.myroms.org/forum/viewtopic.php?f=20&t=4447
mud_mass_init = np.sum(mud_mass[0,:,:,:],axis=0) # total from all 3 bed layers
mud_mass_final = np.sum(mud_mass[-1,:,:,:],axis=0) # total from all 3 bed layers
mud_mass_diff = mud_mass_final - mud_mass_init # kg/m2

sand_mass_init = np.sum(sand_mass[0,:,:,:],axis=0) # total from all 3 bed layers
sand_mass_final = np.sum(sand_mass[-1,:,:,:],axis=0) # total from all 3 bed layers
sand_mass_diff = sand_mass_final - sand_mass_init # kg/m2

# Susquehanna River mouth
trans_name = 'T1'
t1_x = np.array([29,30,31,32])
t1_y = np.array([13,12,11,10])

plant_height = f.variables['mask_rho'][:]
plant_height = np.ma.masked_less(plant_height, 1)
plant_height.harden_mask() # makes the mask immutable

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

# build Flats cell locations:
# This does not include the transect cells themselves. Just all cells between T1 and T2, excluding South River.
# 5 = good points
## Southern Boundary
plant_height[:, 14:58] = 5

# Remove Elk River
plant_height[78:100, 41:59] = 0
plant_height[72:78, 55:59] = 0

## Northern Boundary
plant_height[48:64,0:14] = 5 # add mill creek/Furnace Bay
# dealing with angled transect boundary here
plant_height[30:50, 13] = 5
plant_height[31:37, 12] = 5
plant_height[32:35, 11] = 5

#plt.figure()
#plt.plot_date(datetime_list, np.sum(bed_thick[:, :, 32, 15], axis=1), markersize=1, linewidth=1)
#plt.ylabel('total bed thickness [m]')
#plt.title('Bed thickness time series at point index (32,15)')
#plt.plot_date(datetime_list,integrate.trapz(bed_thick[:,:,32,15],axis=1),label='trapz',markersize=1,linewidth=1)
#plt.legend()
#sys.exit()
## Plotting
# apply the mask
# plant_height = 5 is where the region of interest is.
# so apply a mask to everything not 5 to the bed thick matrix

bed_thick_diff_ma = np.ma.masked_where(plant_height != 5, bed_thick_diff)
bed_deposition = np.ma.masked_less(bed_thick_diff_ma, 0) # height in meters
bed_erosion = np.ma.masked_greater(bed_thick_diff_ma, 0) # height in meters

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

## using volume
#bed_dep_vol = bed_deposition * SAm#.mean() # vol deposited m^3
#bed_ero_vol = bed_erosion * SAm#.mean() # vol eroded m^3
#mass_deposited = Srho * bed_dep_vol  # kg/m^3 * m^3 = kg # TODO multiply by dry bulk density using Srho and porosity
#mass_eroded = Srho * bed_ero_vol  # kg/m^3 * m^3 = kg # TODO multiply by dry bulk density using Srho and porosity

## using sed mass
mass_eroded = sand_mass_eroded + mud_mass_eroded
mass_deposited = sand_mass_deposited + mud_mass_deposited

print('mass eroded    = %e tons' % (np.sum(mass_eroded) / 1000))
print('mass deposited = %e tons' % (np.sum(mass_deposited) / 1000))


## Plotting
plt.figure()
plt.pcolor(f.variables['lon_rho'][:], f.variables['lat_rho'][:], plant_height)
plt.title('Transects')

# set up figure
fig, ax = plt.subplots(figsize=(8, 6))

# set up map
m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
    resolution='i', projection='merc', ax=ax)
m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax)
m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax)
#m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

# pcolor variable of interest
cax = m.pcolormesh(lon, lat, mud_mass_diff_ma, latlon=True,
                    cmap='jet', ax=ax)
cbar = fig.colorbar(cax)
cbar.set_label('mud mass [kg/m2]')
plt.title('mud mass evolution')
#plt.title("%s through %s"%(datetime_list[0],datetime_list[-1]))

# set up figure
fig, ax = plt.subplots(figsize=(8, 6))

# set up map
m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
    resolution='i', projection='merc', ax=ax)
m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax)
m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax)
#m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

# pcolor variable of interest
cax = m.pcolormesh(lon, lat, sand_mass_diff_ma, latlon=True,
                    cmap='jet', ax=ax)
cbar = fig.colorbar(cax)
cbar.set_label('sand mass [kg/m2]')
plt.title('sand mass evolution')
#plt.title("%s through %s"%(datetime_list[0],datetime_list[-1]))

# set up figure
fig, ax = plt.subplots(figsize=(8, 6))

# set up map
m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
    resolution='i', projection='merc', ax=ax)
m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax)
m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax)
#m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

# pcolor variable of interest
cax = m.pcolormesh(lon, lat, mass_eroded, latlon=True,
                    cmap='jet', ax=ax)
cbar = fig.colorbar(cax)
cbar.set_label('Bed mass eroded [kg]')
plt.title('mass eroded = %e tons' % (np.sum(mass_eroded) / 1000))
#plt.title("%s through %s"%(datetime_list[0],datetime_list[-1]))

# set up figure
fig, ax = plt.subplots(figsize=(8, 6))

# set up map
m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
    resolution='i', projection='merc', ax=ax)
m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax)
m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax)
#m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

# pcolor variable of interest
cax = m.pcolormesh(lon, lat, mass_deposited, latlon=True,
                    cmap='jet', ax=ax)
cbar = fig.colorbar(cax)
cbar.set_label('Bed mass deposited [kg]')
plt.title('mass deposited = %e tons' % (np.sum(mass_deposited) / 1000))
#plt.title("%s through %s"%(datetime_list[0],datetime_list[-1]))