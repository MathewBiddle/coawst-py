'''
Biddle et al 2020 manuscript Figure 1:
Caption:
Map of study region with points denoting the locations of USGS Conowingo Dam discharge observations (C),
turbidity observations from Maryland Department of Natural Resources Chesapeake Bay-Segment 1 Susquehanna Flats (FLT)
and Susquehanna River – Havre de Grace (SUS), wind and current observations (CBIBS), wave observations from 2013
(Tripod), surface elevation prediction (SHAD), and surface elevation measurement (TOL). Lines T1 and T2 are the
transects used in the flux calculations at the entrance and exit of the system, respectively. The colors indicate
the NOAA NOS MLLW bathymetry, in meters, scaling from shallow water (yellow) to deep water (green), used for model
initialization. The vegetation patch supplied to the model is the grey shaded region. Model grid included...

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

locs = coawstpy.get_point_locations()

# setup Lambert Conformal basemap.
lat_max = locs.loc[locs['Site'] == '3','lat']+.16#08
lat_min = locs.loc[locs['Site'] == 'TOL','lat']-.016
lon_min = locs.loc[locs['Site'] == 'TOL','lon']-.01
lon_max = locs.loc[locs['Site'] == '3','lon']+.11

# build map
#fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))

m = Basemap(llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max,
    resolution='i', projection='merc', lat_ts=locs.loc[locs['Site'] == '3', 'lat'], epsg=3395)

# add background
# arcgisimage function doesn't work
# image_url = m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels=3000, verbose=True)
# import urllib.request
# urllib.request.urlretrieve(image_url,'map.png')
#
# Download image, load it, and display it on the map
#image_url = 'http://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/export?bbox=-8488667.770441078,4722963.927554283,-8453713.45033199,4791541.555979901&bboxSR=3395&imageSR=3395&size=3000,5885&dpi=96&format=png32&transparent=true&f=image'
import matplotlib.image as mpimg
img = mpimg.imread('C:\\Users\\Mathew.Biddle\\Documents\\GitProjects\\Thesis_data\\map.png')
m.imshow(np.flipud(img))
#mount_point = 'Z:\matt_backups' # windows mapping of external drive
#mount_point = '/Users/mbiddle' # local on mac machine
#root = '\matt_backups\Documents\BCO-DMO\Graduate_School\Thesis\COAWST\COAWST_RUNS\COAWST_OUTPUT\Full_20110719T23_20111102_final'

#'\Documents\BCO-DMO\Graduate_School\Thesis\COAWST\COAWST_RUNS\COAWST_OUTPUT\Full_20110719T23_20111101_final'

#inputfile = direct+'\\upper_ches_his.nc'
inputfile = coawstpy.get_file_paths()['veg']
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
m.contour(lon, lat, h, [3], linewidths=0.5, linestyles=':', colors='k', latlon=True)
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

# Make the plot fancy for manuscript

# labels = [left,right,top,bottom]
#parallels = np.arange(lat_min.get_values(),lat_max.get_values(),0.1)
parallels = np.linspace(lat_min,lat_max-0.001,6).flatten()
m.drawparallels(parallels,labels=[True,False,True,True],fmt='%4.2f', fontdict=dict(size=5),linewidth=0.5,
                dashes=[4, 900])
#meridians = np.arange(lon_min.get_values(),lon_max.get_values(),0.1)
meridians = np.linspace(lon_min,lon_max,4).flatten()
m.drawmeridians(meridians,labels=[True,False,False,True],fmt='%4.2f', fontdict=dict(size=5),linewidth=0.5,
                dashes=[4, 900])

#m.ax.set_xticks([0,1,2])
# scale
m.drawmapscale(lon=lon_max.max()-0.05, lat=lat_min.min()+0.02,
               lon0=lon_max.max(), lat0=lat_min.min(),
               units='km', length=5, yoffset=530, linewidth=0.2, fontsize=5, barstyle='simple')

plt.title('Site locations and SAV distribution', fontdict=dict(size=5))

#outfile = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/figures/Site_locations_CBOFS_conowingo.png'
#outfile = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/Manuscript/figures/Fig_1.png'
outfile = 'Fig_1.png'
plt.savefig(outfile, bbox_inches='tight', dpi = 500)