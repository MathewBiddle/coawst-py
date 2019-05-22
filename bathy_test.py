# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 20:25:09 2017

@author: matt
"""

import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
# setup Lambert Conformal basemap.
lat_max = 40.066667
lat_min = 35.891667
lon_min = -78.266667
lon_max = -74.091667
m = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max,
             resolution='i', projection='merc', lat_0 = 39.5,lon_0 = 1)
# draw coastlines.
m.drawcoastlines()
#m.etopo()
# draw a boundary around the map, fill the background.
# this background will end up being the ocean color, since
# the continents will be drawn on top.
#m.drawmapboundary(fill_color='aqua')
# fill continents, set lake color same as ocean color.
m.fillcontinents(color='coral',lake_color='aqua')
plt.show()
