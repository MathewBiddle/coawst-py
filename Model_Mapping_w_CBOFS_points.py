# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 20:25:09 2017

@author: matt
"""

import os
os.environ["PROJ_LIB"] = "/Users/mbiddle/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import netCDF4
import coawstpy

locs = coawstpy.get_point_locations()

# setup Lambert Conformal basemap.
lat_max = locs.loc[locs['Site'] == '3','lat']+.16#08
lat_min = locs.loc[locs['Site'] == 'TOL','lat']-.016
lon_min = locs.loc[locs['Site'] == 'TOL','lon']-.01
lon_max = locs.loc[locs['Site'] == '3','lon']+.1

# build map
#fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))

m = Basemap(llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max,
    resolution='i', projection='merc', lat_ts=locs.loc[locs['Site'] == '3', 'lat'], epsg=3395)

# add background
m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels=3000)

direct = \
'/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'

inputfile = direct+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')
lon = f.variables['lon_rho'][:][:]
lat = f.variables['lat_rho'][:][:]
h = f.variables['h'][:][:]
mask = f.variables['mask_rho'][:][:]
mask = np.ma.masked_equal(mask, 0)
h.mask = mask.mask # apply mask
plant_height_raw = f.variables['plant_height'][0, 0, :, :]
plant_height_orig = np.ma.masked_outside(plant_height_raw,0.3,1)
#plant_height_domain = np.ma.masked_greater(plant_height_raw, 1)
#plant_height = np.ma.masked_less(plant_height, 0.3)

# shift plant distribution to right by one cell
#zeroes = np.zeros((1, 100))
#plant_height = np.concatenate((zeroes,plant_height_orig),0)
#plant_height = np.delete(plant_height, 100, 0)
#plant_height = np.ma.masked_outside(plant_height,0.3,1)
plant_height = plant_height_orig

# plot bathymetry
m.pcolormesh(lon, lat, h, latlon=True, cmap='viridis_r')
cbar = m.colorbar()
cbar.ax.set_ylabel('NOS MLLW bathymetry [meters]', rotation=270, fontdict=dict(size=5))
cbar.ax.tick_params(labelsize=5)

# plot pant height
m.pcolormesh(lon, lat, plant_height, latlon=True, cmap='binary',vmin=0,vmax=0.3, alpha=0.3,linewidth=0)
#plant_height_domain = np.concatenate((zeroes,plant_height_domain),0)
#plant_height_domain = np.delete(plant_height_domain, 100, 0)
#m.contour(lon, lat, plant_height_domain, 0.305, colors='k', linewidths=.5, latlon=True, corner_mask=False)
#m.pcolormesh(lon, lat, plant_height_orig, latlon=True, cmap='binary',vmin=0,vmax=0.3, alpha=0.3,linewidth=0)

## Do mapping for points
#lon=list(locs['lon'])
#lat=list(locs['lat'])

#x,y = m(lon,lat)

# plot stations
#plt.scatter(x,y,s=50,marker='.',color='r',edgecolors='k',linewidths=0.3)
for label in locs['Site']:
    if label in ['CBIBS','Tripod','SUS','FLT','SHAD','TOL']:
        lon = locs.loc[locs['Site'] == label,'lon'].values
        lat = locs.loc[locs['Site'] == label,'lat'].values
        x,y = m(lon,lat)
        plt.scatter(x, y, s=40, marker='.', color='r', edgecolors='k', linewidths=0.3)
        if label == 'Tripod':
            plt.text(x+500,y-200, label, fontdict=dict(size=5))
        elif label == 'SUS':
            plt.text(x-2500,y-500, label, fontdict=dict(size=5))
        elif label == 'SHAD':
            plt.text(x+500, y-500, label, fontdict=dict(size=5))
        else:
            plt.text(x+500, y, label, fontdict=dict(size=5))

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
    plt.text(xm[-1]+500,ym[0],t_name, fontdict=dict(size=5))

# add conowingo sensor location
latc = 39.658
lonc = -76.1744
xc, yc = m(lonc, latc)
plt.scatter(xc, yc, s=40, marker='.', color='r', edgecolors='k', linewidths=0.3)
plt.text(xc-1500,yc-500, 'C', fontdict=dict(size=5))
plt.title('Site locations and SAV distribution', fontdict=dict(size=5))

outfile = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/figures/Site_locations_CBOFS_conowingo.png'

plt.savefig(outfile, bbox_inches='tight', dpi = 500)