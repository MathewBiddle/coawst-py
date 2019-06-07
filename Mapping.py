# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 20:25:09 2017

@author: matt
"""

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

# setup Lambert Conformal basemap.
lat_max = locs.loc[locs['Site'] == '3','lat']+.08
lat_min = locs.loc[locs['Site'] == '3','lat']-.16
lon_min = locs.loc[locs['Site'] == '3','lon']-.1
lon_max = locs.loc[locs['Site'] == '3','lon']+.1
#m = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max,
#             resolution='i', projection='merc', lat_ts=locs['lat'].mean())
m = Basemap(llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max,
    resolution='i', projection='merc', lat_ts=locs.loc[locs['Site'] == '3', 'lat'], epsg=3395)

#lon=list(locs['lon'])
#lat=list(locs['lat'])

#x,y = m(lon,lat)
# draw coastlines.
m.arcgisimage(service="Canvas/World_Light_Gray_Base",xpixels = 3000)

#,5,color='r',edgecolors='k',latlon=True,linewidths=0.3)

# for label, xpt, ypt in zip(locs['Site'], locs['lon'], locs['lat']):
#     if label != 'Lee5':
#         plt.annotate(xpt, ypt, label, size='x-small')

#dir='/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110702_20111101'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110906_20110926'
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
#dir = '/Volumes/Documents/COAWST_34_UPPER_CHES_FULL'
#inputfile = dir+'/upper_ches_his.nc'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110714201800_20111031231800'
inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')
lon = f.variables['lon_rho'][:][:]
lat = f.variables['lat_rho'][:][:]
h = f.variables['h'][:][:]
mask = f.variables['mask_rho'][:][:]
mask = np.ma.masked_equal(mask,0)
h.mask = mask.mask
plant_height = f.variables['plant_height'][0, 0, :, :]
plant_height = np.ma.masked_outside(plant_height,0.3,1)
#plant_height = np.ma.masked_greater(plant_height, 1)
#plant_height = np.ma.masked_less(plant_height, 0.3)
m.pcolormesh(lon, lat, h, latlon=True, cmap='viridis_r')
cbar = m.colorbar()
cbar.ax.set_ylabel('NOS MLLW bathymetry [meters]', rotation=270, fontdict=dict(size=5))
cbar.ax.tick_params(labelsize=5)
m.pcolormesh(lon, lat, plant_height, latlon=True, cmap='binary',vmin=0,vmax=0.3, alpha=0.3,linewidth=0)

lon=list(locs['lon'])
lat=list(locs['lat'])

x,y = m(lon,lat)

plt.scatter(x,y,s=50,marker='.',color='r',edgecolors='k',linewidths=0.3)
for label, xpt, ypt in zip(locs['Site'], x, y):
    if label not in ['Lee5','Lee2.5','Lee2','Lee0','LeeS2']:
        if label == 'Tripod':
            plt.text(xpt+500,ypt-200, label, fontdict=dict(size=5))
        else:
            plt.text(xpt+500, ypt, label, fontdict=dict(size=5))


plt.title('Site locations and SAV distribution',fontdict=dict(size=5))

outfile = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/figures/Site_locations.png'

plt.savefig(outfile, bbox_inches='tight', dpi = 1000)