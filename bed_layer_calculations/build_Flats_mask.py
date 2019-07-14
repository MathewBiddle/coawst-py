import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import netCDF4

## Read COAWST data
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')
lon = f.variables['lon_rho'][:][:]
lat = f.variables['lat_rho'][:][:]
ocean_time = f.variables['ocean_time'][:]
bed_thick = f.variables['bed_thickness'][:]
Srho = f.variables['Srho'][0] # kg/m3 2650 kg/m^3 sediment grain density
pm = f.variables['pm'][:] #XI --> cell width in x dir.
pn = f.variables['pn'][:] #ETA --> cell width in y dir. Want to use this for Surface Area Calcs

SA = (1/pm) * (1/pn)
## Do some date conversions ##
datetime_list=[]
for sec in ocean_time:
    datetime_list.append(
        netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))


bed_thick_init = np.sum(bed_thick[0,:,:,:],axis=0) # total from all 3 bed layers
bed_thick_final = np.sum(bed_thick[-1,:,:,:],axis=0) # total from all 3 bed layers
bed_thick_diff = bed_thick_final - bed_thick_init # meters

# Susquehanna River mouth
trans_name = 'T1'
t1_x = np.array([29,30,31,32])
t1_y = np.array([13,12,11,10])

plant_height = f.variables['plant_height'][0, 0, :, :]
plant_height = np.ma.masked_greater(plant_height, 1)
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


# apply the mask
# plant_height = 5 is where the region of interest is.
# so apply a mask to everything not 5 to the bed thick matrix
bed_thick_diff_ma = np.ma.masked_where(plant_height != 5, bed_thick_diff)

bed_deposition = np.ma.masked_less(bed_thick_diff_ma, 0) # height in meters
bed_erosion = np.ma.masked_greater(bed_thick_diff_ma, 0) # height in meters

bed_dep_vol = bed_deposition * SA.mean() # vol deposited m^3
bed_ero_vol = bed_erosion * SA.mean() # vol eroded m^3

mass_deposited = Srho * bed_dep_vol  # kg/m^3 * m^3 = kg
mass_eroded = Srho * bed_ero_vol  # kg/m^3 * m^3 = kg

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
cax = m.pcolormesh(lon, lat, mass_deposited, latlon=True,
                    cmap='jet', ax=ax)
cbar = fig.colorbar(cax)
cbar.set_label('Bed mass deposited [kg]')
plt.title("%s through %s"%(datetime_list[0],datetime_list[-1]))


