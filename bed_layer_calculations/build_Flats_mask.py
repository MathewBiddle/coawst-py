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
plant_height.harden_mask()
#mk = plant_height.mask
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

# build flats inner cell indexes:
## Southern Boundary
plant_height[:, 13:59] = 5

## Northern Boundary
plant_height[30:37, 12] = 5
plant_height[31:35, 11] = 5
plant_height[32, 10] = 5

# South river
plant_height[78:100, 41:59] = 0
plant_height[72:78, 55:59] = 0

#sf_mask = np.ma.masked_equal(plant_height,5)

bed_thick_diff_ma = np.ma.masked_where(plant_height != 5, bed_thick_diff)

plt.figure()
plt.pcolor(f.variables['lon_rho'][:], f.variables['lat_rho'][:], plant_height)
plt.title('Transects')

# set up figure
fig, ax = plt.subplots(figsize=(8, 6))

#set up map
m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
    resolution='i', projection='merc', ax=ax)
m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax)
m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax)
#m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

# pcolor variable of interest
cax = m.pcolormesh(lon, lat, bed_thick_diff_ma*100, latlon=True,
                    vmin=-10, vmax=10, cmap='jet', ax=ax)
cbar = fig.colorbar(cax)
cbar.set_label('Bed evolution [cm]')
plt.title("%s through %s"%(datetime_list[0],datetime_list[-1]))
