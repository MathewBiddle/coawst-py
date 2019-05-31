#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 15:41:50 2014

usage: import_netcdf.py [filename,opendap path]

purpose: To generate an initial plot of netCDF data.

example: import_netcdf.py GOSUD_GO_WTED_2013_TRAJ_converted.nc
	 import_netcdf.py http://data.nodc.noaa.gov/thredds/dodsC/nmsp/bml/BOD001_000_20050627_20051109.nc

Updated: 2014-06-10
@author: mbiddle
"""
#%tb
import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
#import warnings
import numpy as np
#import numpy.ma as ma
import datetime
import time
import netCDF4
import sys
#import os

## Define haversine function to claculate distance between two points
def haversine(lon1, lat1, lon2, lat2):
   """
   Calculate the great circle distance between two points
   on the earth (specified in decimal degrees) in meters.
   """
   from math import radians, cos, sin, asin, sqrt
   # convert decimal degrees to radians
   lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
   # haversine formula
   dlon = lon2 - lon1
   dlat = lat2 - lat1
   a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
   c = 2 * asin(sqrt(a))
   # 6378.1 km is the equatorial radius of the Earth from http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
   km = 6378.1 * c
   m = 1000 * km
   return m


#inputfile = sys.argv[1] # get filename
#inputfile = 'https://data.nodc.noaa.gov/thredds/dodsC/testdata/mbiddle/veg_test_his_compressed.nc'
#inputfile = 'https://data.nodc.noaa.gov/thredds/dodsC/testdata/mbiddle/upper_ches_his.nc'
#inputfile = '/Users/mbiddle/Documents/thesis/working_directory/upper_ches_his_translated.nc'
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110702_20111101'
inputfile = dir+'/upper_ches_his.nc'
#inputfile = '\\USERS\mbiddle\umd\thesis\data'
## -------netcdf file processing------#
f = netCDF4.Dataset(inputfile,'r') # open the netcdf file
#print "\n",inputfile
#for v in f.variables: # print each variable of the netCDF file
#   if hasattr(f.variables[v],'units') and hasattr(f.variables[v],'long_name'):
#      print f.variables[v].long_name,';',f.variables[v].units
#   elif hasattr(f.variables[v],'long_name'):
#      print f.variables[v].long_name
obs = 'zeta'#'plant_diameter'#'plant_thickness'#'plant_height'#'Dissip_veg'#'Dissip_fric'#'zeta'#'Hwave'#'Lwave'#'Pwave_bot'#'Pwave_top'#'Dissip_break'#'plant_height'#'h'#'temp'#'zeta'#'temp' # the name of the variable
#ocn_time = 'ocean_time'
plt_type = 'contour'#'cross-shore'
#time_var = f.variables[ocn_time]
time_data = f.variables['ocean_time'][:]
obs_var = f.variables[obs]
x = obs_var.shape[-1] # xi_rho
y = obs_var.shape[-2] # eta_rho

## Do some date conversions ##
epoch_date = '%s %s'%(f.variables['ocean_time'].units.split(' ')[-2],f.variables['ocean_time'].units.split(' ')[-1])
dt_obj = datetime.datetime.strptime(epoch_date, '%Y-%m-%d %H:%M:%S')
time_diff = time.mktime(dt_obj.timetuple())-time.mktime(datetime.datetime(1970,1,1,0,0,0).timetuple())

## set plot ranges according to Beudin 2017.
#if plt_type == 'contour':
#   if obs == 'Hwave':
#      bounds = np.arange(0,0.5,.005) # set min max interval for colorbar
#   elif 'Pwave' in obs:
#      bounds = np.arange(1,2.5,.005)
#   elif obs == 'zeta':
#      bounds = np.arange(-0.5,0.6,0.01)
#elif plt_type == 'cross-shore':
#   if obs == 'Hwave':
#      bounds = [0,0.3]
#   if 'Dissip' in obs:
#      bounds = [0,0.3]
#coord_1_data = range(x)

## Create tide timeseries
#H = np.float(0.5)
#tide=[]
datetime_list=[]
for sec in time_data:
   ts = sec/(12*3600)
   #tide.append(-H*np.sin(2*np.pi*ts))
   datetime_list.append(datetime.datetime.fromtimestamp(sec+time_diff))
lat = f.variables['lat_rho'][:][:]
lon = f.variables['lon_rho'][:][:]
#data = f.variables['tke'][24][20][:][:]
#data = f.variables['h'][:][:]
veg = f.variables['plant_height'][0][0][:][:]
mask_veg = np.zeros(veg.shape)
mask_veg[:] = np.NAN
veg_idx = np.where(veg==0.3048)
mask_veg[veg_idx]=1
veg_lat = lat[veg_idx]
veg_lon = lon[veg_idx]

#upperleft
#lowerleft
#lowerright
#upperleft
veg_verts = [(np.min(veg_lon),np.max(veg_lat)),\
             (np.min(veg_lon),np.min(veg_lat)),\
             (np.max(veg_lon),np.min(veg_lat)),\
             (np.max(veg_lon),np.max(veg_lat)),\
             (np.min(veg_lon),np.max(veg_lat))] 
#print np.transpose(veg_verts)
#print type(np.array(veg_verts))
min_lat = np.min(lat)-.01
max_lat = np.max(lat)+.01
min_lon = np.min(lon)-.01
max_lon = np.max(lon)+.01

# Calculate grid size in meters
grid_distance_meters_y = haversine(np.min(lon),np.max(lat),np.min(lon),np.min(lat))
grid_distance_meters_x = haversine(np.min(lon),np.max(lat),np.max(lon),np.max(lat))
veg_distance_meters_y = haversine(np.min(veg_lon),np.max(veg_lat),np.min(veg_lon),np.min(veg_lat))
veg_distance_meters_x = haversine(np.min(veg_lon),np.max(veg_lat),np.max(veg_lon),np.max(veg_lat))

print("Grid distances (km):")
print('x =',grid_distance_meters_x/1000) #km
print('y =',grid_distance_meters_y/1000) #km
print("Vegetation distances (km):")
print('x =',veg_distance_meters_x/1000) #km
print('y =',veg_distance_meters_y/1000) #km

m = Basemap(projection='merc',llcrnrlat=min_lat,urcrnrlat=max_lat,\
            llcrnrlon=min_lon,urcrnrlon=max_lon,resolution='l')

m.plot(np.min(veg_lon),np.max(veg_lat),'o',latlon=True)
m.plot(np.min(veg_lon),np.min(veg_lat),'o',latlon=True)
m.plot(np.max(veg_lon),np.max(veg_lat),'o',latlon=True)
#m.plot(np.transpose(veg_verts)[0],np.transpose(veg_verts)[1],\
#       '-',latlon=True,color='red',linewidth=2)

m.pcolormesh(lon,lat,np.ones(lon.shape),\
             latlon=True,linewidth=0.005,facecolor='none',edgecolor='black')#,cmap=cm.terrain)#,linewidth=0.0035,color='black')
#print min_veg_lon
#m = Basemap(projection='merc',resolution='f')
#m.pcolormesh(lon,lat,f.variables['mask_rho'][:][:],\
#             latlon=True,cmap=cm.terrain)#,linewidth=0.0035,color='black')

#m.pcolormesh(lon,lat,f.variables['tke'][:][:][25][20]*f.variables['mask_rho'][:][:],\
#             latlon=True,linewidth=0,cmap='terrain')#.007,color='black') #f.variables['h'][:][:]*
#m.pcolormesh(lon,lat,veg,\
#             latlon=True,linewidth=0.007,color='black') #f.variables['h'][:][:]*

#m.drawcoastlines(linewidth=4)
m.drawparallels(np.arange(39.3,39.7,.1),labels=[1,1,0,0])
m.drawmeridians(np.arange(-76.2,-75.8,.1),labels=[0,0,1,1])
m.fillcontinents(color='grey',zorder=0)

#plt.clim(np.min(f.variables['tke'][:][:][25][20]),np.max(f.variables['tke'][:][:][25][20])-0.0235)

#cbar = plt.colorbar(format='%.0e')
#cbar = plt.colorbar()
#cbar.formatter.set_powerlimits((0,0))
#cbar.ax.text(0, 1, r'$\times$10$^{-4}$', va='bottom', ha='left')
#cbar.update_ticks()
#cbar.set_label('%s [%s]' % (f.variables['tke'].long_name,f.variables['tke'].units),fontsize=15) # Label colorbar

#title = '%s @ %3.3f'%(datetime.datetime.fromtimestamp(f.variables['ocean_time'][20]+time_diff).strftime('%x %X'),f.variables['s_rho'][20])
#plt.title(title)

#from matplotlib.patches import Polygon
#patch = Polygon(veg_verts, linewidth=2,color='green')
#plt.gca().add_patch(patch)
#ax.grid(True, which='minor', axis='both', linestyle='-', color='k')
plt.show()
print("Done")
#sys.exit()
## Create surface plot and tide timeseries
#obs='zeta'
for t in range(len(time_data)):
   obs_data = f.variables[obs][t][:][:] # @ time=t,get all (y,x)
   print(obs_data.shape)
## Top subplot ##
   fig = plt.figure(figsize=(10,15))
   ax = fig.add_subplot(2,1,1)#2,1,1)

   if plt_type == 'contour':
#      if 'bounds' in locals():
#         CS = plt.contourf(range(x), range(y), obs_data,bounds,cmap=cm.jet,extend='both') #Pwave_top YlGnBu_r
#      else:
      CS = plt.contourf(range(x), range(y), obs_data,cmap=cm.jet,extend='both')
      for c in CS.collections: # hide contour linestroke
         c.set_edgecolor("face")
      cbar = plt.colorbar(format='%d') # format number display on colorbar
      cbar.set_label('%s' % (obs_var.units),fontsize=10) # Label colorbar
      plt.xlim(0,np.max(x)) # Set x limit to actual max(x)
   elif plt_type == 'cross-shore':
#      print obs_data[:,44].shape
      plt.plot(range(y),obs_data[:,44])
      plt.grid(True)
      plt.ylim(bounds)
      plt.ylabel('%s (%s)'%(obs_var.long_name,obs_var.units))
      plt.xlabel('cross-shore distance (m)')
   title = '%s\n%s'%(obs_var.long_name,datetime.datetime.fromtimestamp(f.variables['ocean_time'][t]+time_diff).strftime('%x %X'),)
   plt.title(title)

## Bottom subplot ##
   ax = fig.add_subplot(2,1,2)
   plt.plot_date(datetime_list,tide,'o') #plot tide function
   ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
   plt.axvline(x=datetime_list[np.int(t)]) # add vertical line to tide plot

   plt.show()
   sys.exit()
#   plt.savefig('%s_%03d.png'%(obs,t))
#   if t == 2:
#      sys.exit()

#   plt.show()
   fig.clf()
   print(t)
#   if t == 20:
   sys.exit()
sys.exit()




print('Coordinates:',obs_var.coordinates)
i=0
for coord in obs_var.coordinates.split(' '):
   print(coord)
   string_meta = 'coord_%s = f.variables[\'%s\']' %(i,coord)
   exec(string_meta)
   string_data = 'coord_%s_data = f.variables[\'%s\'][:]' %(i,coord)
   print(string_data)
   shape = 'print coord_%s.shape' %i
   exec(shape)
   exec(string_data)
   data_print = 'print coord_%s_data'%i
   exec(data_print)
   i+=1
#sys.exit()
#print obs_var.shape # temp dimension = time (121), depth(16), lon(62,202), lat(62,202) 
#print len(obs_data[1,1,1]) # 202 = ?
#print len(obs_data[:,1,1]) # 121 = all time at one depth
#print len(obs_data[1,:,1]) # 16 = all depth at one time
#print len(obs_data[1,1,:]) # 62 = ?
#sys.exit()
print(coord_0_data[1,:].shape)
print(coord_1_data.shape)
print(obs_data[:,:].shape)
#print np.transpose(obs_data[1,:,1,:]).shape

for i in range(0,obs_data.shape[0]): # plot going by ocean_time
   ## need to pick a latitude line to follow for contour plot
   fig = plt.figure()
   ax = fig.add_subplot(111)
   surf = ax.contour(coord_0_data, coord_1_data, obs_data)
   ax.set_xlabel('%s (%s)' %(coord_0.long_name,''),fontsize=10)
   ax.set_ylabel('%s (%s)' % (coord_1.long_name,''),fontsize=10)
   fig.suptitle('%s: %s'%(coord_1.long_name,coord_1.units))#,coord_3.long_name,coord_3[i])) #change time in title
   surf.set_clim(vmin=np.min(obs_data), vmax=np.max(obs_data))
   cbar = fig.colorbar(surf)
#   cbar.set_label('%s (%s)' % (obs_var.long_name,obs_var.units),fontsize=10)
   cbar.set_label('%s (%s)' % (obs_var.long_name,''),fontsize=10)
#   plt.show()
   plt.savefig('test_%03d.png'%i)
   sys.exit()
sys.exit()

## Trying to use code from http://estuary.link/vertical-view-for-roms-output/
#
#stations = np.asarray([13+i for i in range(20)]) # 14 - 33
stations = np.asarray([14,15])
#stations = np.asarray(stations)
time = 100
t=time
vname='salt'
#if type(time) == int:
#   t = time
#   time = netCDF4.num2date(f.variables['ocean_time'][t], romspy.JST)
#elif type(time) == datetime.datetime:
#   time2 = netCDF4.date2num(time, romspy.JST)
#   time3 = f.variables['ocean_time'][:]
#   t = np.where(time3==time2)[0][0]
#else:
#   print 'ERROR: your time type =',type(time)
#print '\n',time, t


        
cs_r = f.variables['Cs_r'][:]
h = f.variables['h'][stations]
zeta = f.variables['zeta'][t,stations]
lon = f.variables['lon_rho'][stations]
lat = f.variables['lat_rho'][stations]
var = f.variables[vname][t,stations,:]
    
depth = np.zeros([len(stations),len(cs_r)])
dist = np.zeros([len(stations),len(cs_r)])
#from vincenty import vincenty
for s in range(len(stations)):
   depth[s,:] = (h[s] + zeta[s]) * cs_r[:]
   print(h[s].shape)
   print(zeta[s].shape)
   print(cs_r[:].shape)
#   depth[s,:] = h[s]
#   depth = h[s]
   if s == 0:
       dist[s,:] = 0
   else:
       back = [lon[s-1],lat[s-1]]
       fore = [lon[s],lat[s]]
   
       R = 6371  #// radius of the earth in km
       x = (fore[0] - back[0]) * np.cos( 0.5*(fore[1]+back[1]) )
       y = fore[1] - back[1]
       dist = R * np.sqrt( x*x + y*y )
       #dist[s,:] = dist[s-1,:] + vincenty(back, fore).meters

#fig, ax = plt.add_subplots(figsize=[12,4])
fig = plt.figure()
ax = fig.add_subplot(111)
if vname == 'temp': 
   cflevels=np.arange(7,30.1,1.0)
   clevels=cflevels
if vname == 'salt': 
   cflevels=np.arange(23,32.2,0.2)
   clevels=np.arange(23,32.2,1.0)
if vname == 'chlorophyll': 
   cflevels=np.arange(0,20.5,0.5)
   clevels=np.arange(0,21,1.0)
origin = 'upper'
    #origin = 'lower'
CF = plt.contourf(dist/1000, depth, var, levels=cflevels, extend='both', origin=origin)
C = plt.contour(dist/1000, depth, var, colors='k', levels=clevels, origin=origin)
plt.clabel(C, fmt = '%2.1f', colors = 'w')
    #plt.pcolor(dist, depth, var) 　　　　　　　　　　# グリッド状がいいなら
CB = plt.colorbar(CF)
CB.ax.set_ylabel(vname)
    
plt.xlim(0,20)
plt.ylim(-20,0)
    
plt.xlabel('distance(km)')
plt.ylabel('depth(m)')
plt.title(time)
#datetime.datetime.strftime(time, '%Y-%m-%d %H:%M:%S'))
plt.show()  

#stafile = 'Z:/roms/Apps/OB500_fennelP/NL10_copy/ob500_sta.nc'
#stations = [13+i for i in range(20)] # 14 - 33




sys.exit()
for i in range(0,obs_data.shape[0]): # plot going by ocean_time
   ## need to pick a latitude line to follow for contour plot
   fig = plt.figure()
   ax = fig.add_subplot(111)
   surf = ax.contour(coord_0_data[1,:], coord_2_data, obs_data[i,:,1,:])
   ax.set_xlabel('%s (%s)' %(coord_0.long_name,''),fontsize=10)
   ax.set_ylabel('%s (%s)' % (coord_2.long_name,''),fontsize=10)
   fig.suptitle('%s: %s and %s: %s'%(coord_2.long_name,coord_2[1],coord_3.long_name,coord_3[i])) #change time in title
   surf.set_clim(vmin=np.min(obs_data), vmax=np.max(obs_data))
   cbar = fig.colorbar(surf)
#   cbar.set_label('%s (%s)' % (obs_var.long_name,obs_var.units),fontsize=10)
   cbar.set_label('%s (%s)' % (obs_var.long_name,''),fontsize=10)
#   plt.show()
   plt.savefig('test_%03d.png'%i)
sys.exit()



for i in range(0,obs_data.shape[0]): # surface plot going by ocean_time
   fig = plt.figure()
   ax = fig.add_subplot(111)
   surf = ax.pcolor(coord_0_data, coord_1_data,obs_data[i,1,:])
   ax.set_xlabel('%s (%s)' %(coord_0.long_name,coord_0.units),fontsize=10)
   ax.set_ylabel('%s (%s)' % (coord_1.long_name,coord_1.units),fontsize=10)
   fig.suptitle('%s: %s and %s: %s'%(coord_2.long_name,coord_2[1],coord_3.long_name,coord_3[i])) #change time in title
   surf.set_clim(vmin=np.min(obs_data), vmax=np.max(obs_data))
   cbar = fig.colorbar(surf)
#   cbar.set_label('%s (%s)' % (obs_var.long_name,obs_var.units),fontsize=10)
   cbar.set_label('%s (%s)' % (obs_var.long_name,'PSU'),fontsize=10)
#   plt.show()
   plt.savefig('salt_%03d.png'%i)
sys.exit()



for i in range(0,obs_data.shape[1]): # surface plot going by s-coordinate at Rho points (z?)
   fig = plt.figure()
   ax = fig.add_subplot(111)
   surf = ax.pcolor(coord_0_data, coord_1_data,obs_data[:,i,1])
   ax.set_xlabel('%s (%s)' %(coord_0.long_name,coord_0.units),fontsize=10)
   ax.set_ylabel('%s (%s)' % (coord_1.long_name,coord_1.units),fontsize=10)
   fig.suptitle('%s: %s and %s: %s'%(coord_2.long_name,coord_2[i],coord_3.long_name,coord_3[1])) #change s-coordinate in title
   surf.set_clim(vmin=np.min(obs_data), vmax=np.max(obs_data))
   cbar = fig.colorbar(surf)
#   cbar.set_label('%s (%s)' % (obs_var.long_name,obs_var.units),fontsize=10)
   cbar.set_label('%s (%s)' % (obs_var.long_name,'PSU'),fontsize=10)
#   plt.show()
   plt.savefig('salt_%03d.png'%i)
sys.exit()
#print coord_0
sys.exit()


time = f.variables['ocean_time']
#print time
time_data = f.variables['ocean_time'][:]

#bath = f.variables['bath']
bath = f.variables['zeta']
#print bath
bath_data = f.variables['zeta'][:]
#bath_data = -1*f.variables['bath'][:]

#x_rho = f.variables['x_rho']
x_rho = f.variables['lon_rho']
#print x_rho
#x_rho_data = f.variables['x_rho'][:]
x_rho_data = f.variables['lon_rho'][:]

#y_rho = f.variables['y_rho']
y_rho = f.variables['lat_rho']
#print y_rho
#y_rho_data = f.variables['y_rho'][:]
y_rho_data = f.variables['lat_rho'][:]
#print '\ntime:\n', time_data,'\nbath:\n',bath_data,'\nx:\n',x_rho_data,'\ny:\n',y_rho_data

print('\n')
print(bath_data[0,:])
print('\n')
print(bath_data[1,:])
print(np.min(bath_data))
print(time_data[3])
print(time_data.shape)
print(bath_data.shape[0])
#for i in range(0,bath_data.shape[0]):
#   fig = plt.figure()
#   ax = Axes3D(fig)
#   surf = ax.plot_surface(x_rho_data, y_rho_data,bath_data[i,:],cmap=cm.copper)
#   ax.set_xlim3d(np.min(x_rho_data),np.max(x_rho_data))
#   ax.set_xlabel('%s (%s)' %(x_rho.long_name,x_rho.units),fontsize=10)

#   ax.set_ylim3d(np.min(y_rho_data),np.max(y_rho_data))
#   ax.set_ylabel('%s (%s)' % (y_rho.long_name,y_rho.units),fontsize=10)

#   ax.set_zlim3d(np.min(bath_data),np.max(bath_data))
#   ax.set_zlabel('%s (%s)' % (bath.long_name,bath.units),fontsize=10)

#   fig.suptitle('timestep: %s'%time_data[i])
#   surf.set_clim(vmin=np.min(bath_data), vmax=np.max(bath_data))
#   cbar = fig.colorbar(surf)
#   cbar.set_label('%s (%s)' % (bath.long_name,bath.units),fontsize=10)
#   ax.dist = 13

#   plt.savefig('3d_%s.png'%i)

#sys.exit()
for i in range(0,bath_data.shape[0]):
   fig = plt.figure()
   ax = fig.add_subplot(111)
   surf = ax.pcolor(x_rho_data, y_rho_data,bath_data[i,:])
   ax.set_xlabel('%s (%s)' %(x_rho.long_name,x_rho.units),fontsize=10)
   ax.set_ylabel('%s (%s)' % (y_rho.long_name,y_rho.units),fontsize=10)
   fig.suptitle('timestep: %s'%time_data[i])
   surf.set_clim(vmin=np.min(bath_data), vmax=np.max(bath_data))
   cbar = fig.colorbar(surf)
   cbar.set_label('%s (%s)' % (bath.long_name,bath.units),fontsize=10)
#   plt.show()
   plt.savefig('%03d.png'%i)
sys.exit()

resume='n'
while resume.lower()=='n':
   print('What variables would you like to plot?\nLeave z empty if you do not want to plot a third variable.')
   var1 = raw_input('X variable: ') # User selected data to plot
   var2 = raw_input('Y variable: ')
   var3 = raw_input('z variable: ')

   var1_var = f.variables[var1]
   var2_var = f.variables[var2]

   dim1 = [str(x) for x in f.variables[var1].dimensions] # convert unicode tuple to tuple
   print('\n******* %s %s attributes:' %(var1,dim1))
   var1_shape = var1_var.shape
   for att in var1_var.ncattrs():
      print('%s = %s' % (att,getattr(var1_var,att))) # printing the attribute information for each selected variable

   dim2 = [str(x) for x in f.variables[var2].dimensions] # convert unicode tuple to tuple
   print('\n******* %s %s attributes:' %(var2,dim2))
   var2_shape = var2_var.shape
   for att in var2_var.ncattrs():
      print('%s = %s' % (att,getattr(var2_var,att)))


   if var3:
      var3_var = f.variables[var3]
      dim3 = [str(x) for x in f.variables[var3].dimensions] # convert unicode tuple to tuple
      print('\n******* %s %s attributes:' %(var3,dim3))
      var3_shape = var3_var.shape
      for att in var3_var.ncattrs():
         print('%s = %s' % (att,getattr(var3_var,att)))


   print('----------------------\n')
   resume = raw_input('Are these the variables you would like to plot? (y/n) ') # verify variables
var1_data = f.variables[var1][:].flatten()
var2_data = f.variables[var2][:].flatten()
if var3:
   var3_data = f.variables[var3][:]#.flatten()

#print var1_data[0][:],'\n',var2_data
#plt.plot(np.transpose(var1_data,var2_data,'o-')

#plt.title('File: %s' % file.split('/')[-1])
#plt.show()
#sys.exit()

print(var1_data,var2_data)
l = raw_input('What depth level [0-%i]?\nLeave depth level empty if you do not want to plot a surface variable.\n'%(len(var2_data)-1))
if l:
   var2_data_full=var2_data
   var2_data=[]
   print(l,len(var2_data_full)-1)
   print(var2_data_full[int(l):int(len(var2_data_full)-1):4])
   var2_data=var2_data_full[int(l):int(len(var2_data_full)-1):4]


#if l: # pulling out the value for each level
#   var3_data_new = []
#   for v3 in var3_data[0][:][:]:
#      var3_data_new.append(v3[int(l)])
#   var3_data_new = np.array(var3_data_new)

#var2_dataa = [val for sublist in var2_data for val in sublist]
#print var2_data.flatten()
#sys.exit()
#if hasattr(var1_data,'units'):
#   print 'units exists'

#------------Trying to work with _FillValue------------------#
#print 'before'
#print var2_data
#print var2_data.mask
#print type(var2_data)
#if hasattr(var2_var,'_FillValue'):
#   print var2_var._FillValue
#   var2_data = ma.masked_where(f.variables[var2][:] == np.float(var2_var._FillValue),f.variables[var2][:]) 
#   var1_data = ma.masked_where(f.variables[var2][:] ==  np.float(var2_var._FillValue), f.variables[var1][:])
#   var2_data = ma.masked_where(f.variables[var2][:] == np.float(-999.9),f.variables[var2][:])
#   var1_data = ma.masked_where(f.variables[var2][:] ==  np.float(-999.9), f.variables[var1][:])

#for val in var1_data:
#   print val
#   print var1_data.mask
#   print var2_data_mask.mask
#   var1_data = var1_var.set_auto_maskandscale(True)
   #var2_data = var2_var.set_auto_maskandscale(True)
   #mask = [var2_data==var2_var._FillValue]
   #print mask[0]
 #  var2_data = var2_data[~var2_data_mask.mask]
 #  var1_data = var1_data[~var2_data_mask.mask]
#print 'after'
#print mask
#------------------------------------------------------------#
#print var1_data,var2_data
#print var1_shape,var2_shape
#print len(var1_shape),len(var2_shape)
#if not len(var1_shape)==len(var2_shape): # if the dimensions don't match
#   if len(var1_shape)>len(var2_shape): # if var1 has more dimensions than var2
#     var1_data = var1_data[0] 

#if not var1_shape[0] == var2_shape[0]:
#   if var1_shape[0] == var2_shape[1]:
#      var1_data = var1_data
#      var2_data = var2_data[0]
#   elif var1_shape[1] == var2_shape[0]:
#      var1_data = var1_data[0]
#      var2_data = var2_data
#   else:
#      print '%s does not have same dimensions as %s\n'%(var1,var2)
#      sys.exit()


## Plot geographic loaction
if hasattr(var1_var,'standard_name') and "longitude" in f.variables[var1].standard_name.lower():
   print("generating map...")
   #print np.min(var2_data)
   m = Basemap(projection='merc',llcrnrlat=np.min(var2_data)-5,urcrnrlat=np.max(var2_data)+5,\
            llcrnrlon=np.min(var1_data)-10,urcrnrlon=np.max(var1_data)+10,lat_ts=20,resolution='i')
#   m = Basemap(projection='cyl',lat_0=0,lon_0=0,resolution='i')
   m.drawcoastlines()
   m.fillcontinents(color='gray')
   x,y = m(var1_data, var2_data)
#   plt.plot(x, y, 'bo',markersize=5,markeredgecolor='k')
   if l:
      plt.scatter(x, y, c=var3_data_new)
      cb = plt.colorbar()
      cb.set_label('%s [%s]' %(f.variables[var3].long_name, f.variables[var3].units))
   elif var3:
      plt.scatter(x, y, c=var3_data)#, markersize=5)
      cb = plt.colorbar()
      cb.set_label('%s [%s]' %(f.variables[var3].long_name, f.variables[var3].units))
   else:
      plt.plot(x, y, 'bo-',markersize=5,markeredgecolor='k')
   plt.title('File: %s' % file.split('/')[-1])
elif hasattr(var1_var,'long_name') and "longitude" in f.variables[var1].long_name.lower():  
#   print "generating map..."
   m = Basemap(projection='merc',llcrnrlat=np.min(var2_data)-5,urcrnrlat=np.max(var2_data)+5,\
            llcrnrlon=np.min(var1_data)-10,urcrnrlon=np.max(var1_data)+10,lat_ts=20,resolution='i')
   m.drawcoastlines()
   m.fillcontinents(color='gray')
   x,y = m(var1_data, var2_data) 
   plt.plot(x, y, 'bo',markersize=5,markeredgecolor='k')
   plt.title('File: %s' % file)
else: # plot geophysical variables
   fig = plt.figure(figsize=(15,10))
   ax = fig.add_subplot(1,1,1)
   if "seconds since 1970-01-01" in f.variables[var1].units: # if time is x variable, convert seconds since 1970 to python datetime
      map_time = map(datetime.datetime.fromtimestamp,var1_data) # map seconds since 1970 to python date
      time_num = mdates.date2num(map_time) # convert python date to number
      plt.plot_date(x=time_num, y=var2_data,fmt='o-') # plot the timeseries
      ax.xaxis.set_major_formatter(mdates.DateFormatter('%b/%y')) # format the xaxis labels
      plt.title('File: %s' % file.split('/')[-1])
   else: # Just plot the data
      plt.plot(var1_data,var2_data,'o-')
#plt.plot(var1_data,var2_data,'o')
      plt.title('File: %s' % file.split('/')[-1])
# Generate the axes labels
if hasattr(var1_var,'long_name'): # check if variable has long_name
   plt.xlabel('%s [%s]' %(f.variables[var1].long_name, f.variables[var1].units)) #use it in xlabel
elif hasattr(var1_var,'standard_name'): # check if variable has standard_name
   plt.xlabel('%s [%s]' %(f.variables[var1].standard_name, f.variables[var1].units)) # use it in xlabel
else: # otherwise use the variable name
   plt.xlabel('%s [%s]' %(var1, f.variables[var1].units))
if hasattr(var2_var,'long_name'): # check if variable has long_name
   plt.ylabel('%s [%s]' %(f.variables[var2].long_name, f.variables[var2].units)) #use it in xlabel
elif hasattr(var2_var,'standard_name'): # check if variable has standard_name
   plt.ylabel('%s [%s]' %(f.variables[var2].standard_name, f.variables[var2].units)) # use it in xlabel
else: # otherwise use the variable name
   plt.ylabel('%s [%s]' %(var2, f.variables[var2].units))
plt.grid(b=True, which='Major', color='k') # add gridding

plt.show()

#time_min_format = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(min(time1))) # Some fancy formating for the date
#time_max_format = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(max(time1)))
#print time_min_format
#print time_max_format

f.close() # close the file
