'''
This script plots the model grids of the various coordinate systems
'''

import os
os.environ["PROJ_LIB"] = "/Users/mbiddle/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import netCDF4

# bring in the data
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final_noveg'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'
if dir.split("_")[-1] == 'noveg':
    run = "noveg"
else:
    run = "veg"
inputfile = dir+'/upper_ches_avg.nc'
print('Reading %s...' % inputfile.split("/")[-1])
f = netCDF4.Dataset(inputfile, 'r')

lon_min = np.min([
    f.variables['lon_u'][:].min(), f.variables['lon_v'][:].min(),
    f.variables['lon_psi'][:].min(), f.variables['lon_rho'][:].min()
    ])
lon_max = np.min([
    f.variables['lon_u'][:].max(), f.variables['lon_v'][:].max(),
    f.variables['lon_psi'][:].max(), f.variables['lon_rho'][:].max()
    ])

lat_min = np.min([
    f.variables['lat_u'][:].min(), f.variables['lat_v'][:].min(),
    f.variables['lat_psi'][:].min(), f.variables['lat_rho'][:].min()
    ])
lat_max = np.min([
    f.variables['lat_u'][:].max(), f.variables['lat_v'][:].max(),
    f.variables['lat_psi'][:].max(), f.variables['lat_rho'][:].max()
    ])
# set up figure

for grid in ['rho']:
    fig, ax = plt.subplots(figsize=(10, 8))
    m = Basemap(llcrnrlon=lon_min-0.01, llcrnrlat=lat_min-0.02, urcrnrlon=lon_max+0.01, urcrnrlat=lat_max+0.02,
                resolution='i', projection='merc', ax=ax, epsg=3395)
    m.arcgisimage(service="Canvas/World_Light_Gray_Base", xpixels=500)#3000)
    #cax = m.pcolor(f.variables['lon_%s' % grid][:], f.variables['lat_%s' % grid][:], f.variables['mask_%s' % grid][:],
    #               latlon=True, facecolor='none', ax=ax, edgecolors='k', cmap='binary_r')
    cax = m.pcolormesh(f.variables['lon_%s' % grid][:], f.variables['lat_%s' % grid][:], np.ones((100,100)),
                   latlon=True, facecolor='none', ax=ax, edgecolor='black', linewidth=0.005, alpha=0.5)
    #plt.title('%s %s %s' % (f.variables['mask_%s' % grid].long_name,
    #                        f.variables['mask_%s' % grid].location,
    #                        f.variables['mask_%s' % grid].shape))
    outfile = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/figures/grid_%s.png' % grid
    print("Saving to %s" % outfile)
    plt.savefig(outfile, bbox_inches='tight')#, dpi=1000)