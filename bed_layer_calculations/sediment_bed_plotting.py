import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import netCDF4

# Cindy's locations
locs = pd.DataFrame(columns=['Site', 'lat', 'lon'])
locs['Site'] = ['1','2','3','4','5',
                'Lee7','Lee6','Lee5','Lee2.5','Lee2','Lee0','LeeS2',
                'CBIBS','Tripod']
locs['lat'] = [39.527,39.533,39.515,39.505,39.497,
               39.414,39.380,39.346,39.197,39.135,39.061,38.757,
               39.5396,39.4931]
locs['lon'] = [-76.061,-76.061,-76.051,-76.039,-76.036,
               -76.079,-76.088,-76.197,-76.311,-76.328,-76.328,-76.473,
               -76.0741,-76.0341]
locs['comment'] = ['Russ and Palinkas 2018','','Middle of bed','','',
                   '','','','','','','',
                   'CBIBS Susquehanna Flats','Larry tripod site']

# get coords for Site of choice
site = 'Lee6'
lat_pt, lon_pt = locs.loc[locs['Site'] == site, ['lat', 'lon']].values[0]

## Read COAWST data
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')
lon = f.variables['lon_rho'][:][:]
lat = f.variables['lat_rho'][:][:]
x = np.abs(f.variables['lon_rho'][:, 1] - lon_pt).argmin()
y = np.abs(f.variables['lat_rho'][1, :] - lat_pt).argmin()
ocean_time = f.variables['ocean_time'][:]

## Do some date conversions ##
datetime_list=[]
for sec in ocean_time:
    datetime_list.append(
        netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))

# pick data to plot
dvar = 'bed_thickness'
data = f.variables[dvar][:, :, :, :]

data_init=np.sum(data[1202,:,:,:],axis=0)
data_final=np.sum(data[1466,:,:,:],axis=0)
# set up figure
fig, ax = plt.subplots(figsize=(8, 6))

#set up map
m = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(), urcrnrlon=lon.max(), urcrnrlat=lat.max(),
    resolution='i', projection='merc', ax=ax)
m.drawparallels(np.arange(39.3,39.6,0.05),labels=[1,0,0,0],ax=ax)
m.drawmeridians(np.arange(-76.15,-75.90,0.05),labels=[0,0,0,1],ax=ax)
#m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels = 3000)

# pcolor variable of interest
cax = m.pcolormesh(lon, lat, (data_final-data_init)*100, latlon=True,
                    vmin=-10, vmax=10, cmap='jet', ax=ax)
cbar = fig.colorbar(cax)
cbar.set_label('Bed evolution [cm]')
m.scatter(-76.079,39.414,marker='o',color='k',alpha=0.4,edgecolors='k',linewidths=0.3,latlon=True)
m.scatter(-76.088,39.380,marker='o',color='k',alpha=0.4,edgecolors='k',linewidths=0.3,latlon=True)
plt.title("%s through %s"%(datetime_list[1202],datetime_list[1466]))
print('Sediment Deposition @ %s = %f cm' % (site, (data_final[x, y]-data_init[x,y])*100))
