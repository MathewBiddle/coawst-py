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
import coawstpy

# Cindy's locations
# locs = pd.DataFrame(columns=['Site', 'lat', 'lon'])
# locs['Site'] = ['1','2','3','4','5',
#                 'Lee7','Lee6','Lee5','Lee2.5','Lee2','Lee0','LeeS2',
#                 'CBIBS','Tripod']
# locs['lat'] = [39.527,39.533,39.515,39.505,39.497,
#                39.414,39.380,39.346,39.197,39.135,39.061,38.757,
#                39.5396,39.4931]
# locs['lon'] = [-76.061,-76.061,-76.051,-76.039,-76.036,
#                -76.079,-76.088,-76.197,-76.311,-76.328,-76.328,-76.473,
#                -76.0741,-76.0341]
# locs['comment'] = ['Russ and Palinkas 2018','','Middle of bed','','',
#                    '','','','','','','',
#                    'CBIBS Susquehanna Flats','Larry tripod site']

locs = coawstpy.get_point_locations()

# setup Lambert Conformal basemap.
lat_max = locs.loc[locs['Site'] == '3','lat']+.08
lat_min = locs.loc[locs['Site'] == '3','lat']-.16
lon_min = locs.loc[locs['Site'] == '3','lon']-.1
lon_max = locs.loc[locs['Site'] == '3','lon']+.1

# build map
m = Basemap(llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max,
    resolution='i', projection='merc', lat_ts=locs.loc[locs['Site'] == '3', 'lat'], epsg=3395)

# add background
m.arcgisimage(service="Canvas/World_Light_Gray_Base",xpixels = 3000)


dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'

inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')
lon = f.variables['lon_rho'][:][:]
lat = f.variables['lat_rho'][:][:]
h = f.variables['h'][:][:]
mask = f.variables['mask_rho'][:][:]
mask = np.ma.masked_equal(mask, 0)
h.mask = mask.mask # apply mask
plant_height_orig = f.variables['plant_height'][0, 0, :, :]
plant_height_orig = np.ma.masked_outside(plant_height_orig,0.3,1)
#plant_height = np.ma.masked_greater(plant_height, 1)
#plant_height = np.ma.masked_less(plant_height, 0.3)

# shift plant distribution to right by one cell
zeroes = np.zeros((1, 100))
plant_height = np.concatenate((zeroes,plant_height_orig),0)
plant_height = np.delete(plant_height, 100, 0)
plant_height = np.ma.masked_outside(plant_height,0.3,1)

# plot bathymetry
m.pcolormesh(lon, lat, h, latlon=True, cmap='viridis_r')
cbar = m.colorbar()
cbar.ax.set_ylabel('NOS MLLW bathymetry [meters]', rotation=270, fontdict=dict(size=5))
cbar.ax.tick_params(labelsize=5)

# plot pant height
m.pcolormesh(lon, lat, plant_height, latlon=True, cmap='binary',vmin=0,vmax=0.3, alpha=0.3,linewidth=0)
#m.pcolormesh(lon, lat, plant_height_orig, latlon=True, cmap='binary',vmin=0,vmax=0.3, alpha=0.3,linewidth=0)

## Do mapping for points
lon=list(locs['lon'])
lat=list(locs['lat'])

x,y = m(lon,lat)

# plot stations
plt.scatter(x,y,s=50,marker='.',color='r',edgecolors='k',linewidths=0.3)
for label, xpt, ypt in zip(locs['Site'], x, y):
    if label not in ['Lee5','Lee2.5','Lee2','Lee0','LeeS2']:
        if label == 'Tripod':
            plt.text(xpt+500,ypt-200, label, fontdict=dict(size=5))
        else:
            plt.text(xpt+500, ypt, label, fontdict=dict(size=5))


# Transects
transects = coawstpy.get_transect_indexes()
for transect in transects:
    t_name = transect
    x_start = transects[transect]['x'][0]
    y_start = transects[transect]['y'][0]
    x_end = transects[transect]['x'][-1]
    y_end = transects[transect]['y'][-1]
    lon = [f.variables['lon_rho'][x_start,y_start],f.variables['lon_rho'][x_end,y_end]]
    lat = [f.variables['lat_rho'][x_start,y_start],f.variables['lat_rho'][x_end,y_end]]
    xm, ym = m(lon, lat)
    m.plot(xm, ym, '-', color='k', linewidth=2)
    plt.text(xm[0]-1000,ym[0],t_name, fontdict=dict(size=5))

sys.exit()
# identify the points to extract data
# Susquehanna River mouth
t1_name = 'S.R. Mouth'
t1x = np.array([29, 33])
t1y = np.array([13, 9])
t1lon=[f.variables['lon_rho'][t1x[0], t1y[0]], f.variables['lon_rho'][t1x[1], t1y[1]]]
t1lat=[f.variables['lat_rho'][t1x[0], t1y[0]], f.variables['lat_rho'][t1x[1], t1y[1]]]

t1xm, t1ym = m(t1lon, t1lat)
m.plot(t1xm, t1ym, '-', color='k', linewidth=2)
plt.text(np.max(t1xm)+100, np.max(t1ym)+100, 'T1', fontdict=dict(size=5))

# Turkey Point to Sandy Point
t2_name = 'Turkey Pt to Sandy Pt'
t2x = np.array([42,67])  #
t2y = np.array([58, 58])
t2lon=[f.variables['lon_rho'][t2x[0], t2y[0]], f.variables['lon_rho'][t2x[1], t2y[1]]]
t2lat=[f.variables['lat_rho'][t2x[0], t2y[0]], f.variables['lat_rho'][t2x[1], t2y[1]]]

t2xm, t2ym = m(t2lon, t2lat)
m.plot(t2xm, t2ym, '-', color='k', linewidth=2)
plt.text(np.min(t2xm)-1000, np.mean(t2ym), 'T2', fontdict=dict(size=5))


plt.title('Site locations and SAV distribution', fontdict=dict(size=5))

outfile = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/figures/Site_locations.png'

plt.savefig(outfile, bbox_inches='tight', dpi = 1000)