import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
from mpl_toolkits.basemap import Basemap
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
#import warnings
import numpy as np
import datetime
import time
import netCDF4
import coawstpy
#import pandas as pd

# bring in the data
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')

# Get the data we want
ocean_time = f.variables['ocean_time'][:]
lat = f.variables['lat_rho'][:]
lon = f.variables['lon_rho'][:]
h = f.variables['h'][:]
zeta = f.variables['zeta'][:]
mud_01 = f.variables['mud_01'][:]
sand_01 = f.variables['sand_01'][:]
ubar = f.variables['ubar_eastward'][:]
vbar = f.variables['vbar_northward'][:]
u = f.variables['u_eastward'][:]
v = f.variables['v_northward'][:]
pm = f.variables['pm'][:] #XI --> cell width in x dir.
pn = f.variables['pn'][:] #ETA --> cell width in y dir. Want to use this for Surface Area Calcs
s_rho = f.variables['s_rho'][:] # depth levels

## Do some date conversions ##
datetime_list=[]
for sec in ocean_time:
    datetime_list.append(
        netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))


# figure out the transect:

# identify the points to extract data
# Susquehanna River mouth
trans_name_sr = 'T1'
x_sr = np.array([29,30,31,32])
y_sr = np.array([13,12,11,10])

# Turkey Point to Sandy Point
trans_name_t2p = 'T2'
x_t2p = np.array(list(range(42,67)))  #
y_t2p = np.array([58]*len(x_sr))

# Verify point location
plant_height = f.variables['plant_height'][0, 0, :, :]
plant_height = np.ma.masked_greater(plant_height, 1)
for i in range(len(x_sr)):
    plant_height[x_sr[i], y_sr[i]] = i

#plt.figure()
#plt.pcolor(f.variables['lon_rho'][:][:], f.variables['lat_rho'][:][:], plant_height)
#plt.title('%s Transect' % trans_name_sr)

# calculate the total concentration at each time step for the transect (x,y) and depth.
#
# For the transect comparisons, you should be using axial velocity instead of volume to compare the fluxes. 
# So the flux through T1 is the axial component of the velocity in each cell times the predicted SSC in that
# cell times the cross-sectional area of the cell, summed across all cells in the transect.  This calculation
# is done at each time step.  The total sediment load delivered during the storm is then the time integral of
# the flux over the duration of the storm.  Sometimes the fluxes will be negative, during flood tide before
# or after the storm (especially at T2). 
#
#  rotate the components into along and cross channel?
#   YES.  It is a pretty simple calculation, especially if you know the angle from geometrical considerations. 
# The best way to do this when you are not sure of the angle is to find the major axis of the current ellipse. 
# Which is essentially the same thing as doing a least squares fit to the scatter plot of the E-N (or x-y) velocities. 
# You want one angle for the entire cross-section, not for each cell.

# Gather subset data
mud_01_trans_sr = mud_01[:, :, x_sr, y_sr]
sand_01_trans_sr = sand_01[:, :, x_sr, y_sr]
ubar_trans_sr = ubar[:, x_sr, y_sr]
vbar_trans_sr = vbar[:, x_sr, y_sr]
v_trans_sr = v[:, :, x_sr, y_sr]
u_trans_sr = u[:, :, x_sr, y_sr]

angle=[]
ux = np.empty(ubar_trans_sr.shape)
vy = np.empty(ubar_trans_sr.shape)
for t in range(0,ubar_trans_sr.shape[0]):  # time
    angle.append(coawstpy.maj_ax(ubar_trans_sr[t, :], vbar_trans_sr[t, :])) # angle for entire cross-section
    #  ux, uy = coawstpy.rot2xy(ubar_trans_sr[t, :], vbar_trans_sr[t, :], angle)
    for xy in range(0, ubar_trans_sr.shape[1]):  # cell in xy
        # ux - along channel velocity
        # vy - cross channel velocity (positive to the left of the along channel
        ux[t, xy], vy[t, xy] = coawstpy.rot2xy(ubar_trans_sr[t, xy], vbar_trans_sr[t, xy], angle[t])


# plotting
myFmt = DateFormatter("%m/%d")
dayint = 30

fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))

ax1.plot_date(datetime_list, angle, marker='.', linestyle='', markersize=1)
ax1.set_title('Calculated angle for cross-section at %s' % trans_name_sr)
ax1.set_ylabel('Angle')
ax1.xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
ax1.xaxis.set_major_formatter(myFmt)

axpos = int(np.sqrt(ubar_trans_sr.shape[1]))

fig, (ax2) = plt.subplots(nrows=axpos, ncols=axpos, sharex=True, sharey=True, figsize=(12, 8),
                         gridspec_kw={'hspace': 0, 'wspace': 0})
ax2 = ax2.ravel()
for i in range(ubar_trans_sr.shape[1]):
    ax2[i].plot(ubar_trans_sr[:, i], vbar_trans_sr[:, i], marker='.', linestyle='', markersize=1)
    ax2[i].set_title('cell %i' % i)
plt.suptitle('Depth average velocity at %s (mean angle = %f)' % (trans_name_sr, np.mean(angle[1:])))
fig.text(0.5, 0.04, 'Eastward (m/s)', ha='center', va='center')
fig.text(0.06, 0.5, 'Northward (m/s)', ha='center', va='center', rotation='vertical')


fig, (ax3) = plt.subplots(nrows=axpos, ncols=axpos, sharex=True, sharey=True, figsize=(12, 8),
                         gridspec_kw={'hspace': 0, 'wspace': 0})
ax3 = ax3.ravel()
for i in range(ux.shape[1]):
    ax3[i].plot_date(datetime_list, ux[:, i], marker='.', linestyle='', label='ux', markersize=1)
    ax3[i].plot_date(datetime_list, vy[:, i], marker='.', linestyle='', label='vy', markersize=1)
    ax3[i].xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
    ax3[i].xaxis.set_major_formatter(myFmt)
    ax3[i].set_title('cell %i' % i)
ax3[i].legend(loc='lower right')
plt.suptitle('Along (ux) and Cross (vy) channel depth avg. velocities at %s' % trans_name_sr)


sys.exit()
u = f.variables['u_eastward'][:]
u_trans_sr = u[:, :, x_sr, y_sr]

v = f.variables['v_northward'][:]
v_trans_sr = v[:, :, x_sr, y_sr]

for t in range(0,u_trans_sr.shape[0]):  # time
    for z in range(0,u_trans_sr.shape[1]):  # depth
        for xy in range(0,u_trans_sr.shape[2]):  # cell in xy
            angle = coawstpy.maj_ax(u_trans_sr[t, z, xy], v_trans_sr[t, z, xy])

sys.exit()
# calculate water column volume for points x and y
h = f.variables['h'][:]
zeta = f.variables['zeta'][:]
height_sr = h[x_sr, y_sr]+zeta[:, x_sr, y_sr]

#sys.exit()
pm = f.variables['pm'][:]
pn = f.variables['pn'][:]
coord_sr = pm[x_sr, y_sr]*pn[x_sr, y_sr]
vol_sr = height_sr/coord_sr
#sys.exit()
mud_trans_conc_sr = vol_sr*mud_sr
mud_tot_conc_sr = np.sum(mud_trans_conc_sr, axis=1)/1000  # metric ton

sand_trans_conc_sr = vol_sr*sand_sr
sand_tot_conc_sr = np.sum(sand_trans_conc_sr, axis=1)/1000  # metric ton

total_sediment_sr = np.float(np.sum(sand_tot_conc_sr)+np.sum(mud_tot_conc_sr))  # metric ton
print("Transect at the mouth of S.R. total suspended sediment mass = %e tons" % total_sediment_sr)
# totals to 4.62 x 10^6 metric tons entering the bay
#
# From Palinkas 2014:
# Various estimates have been generated for the amount of sediment delivered to the Chesapeake Bay associated with
# TS Lee, ranging from 6.7 to 19.0 X 10^6 t(Cheng et al., 2013; Hirsch, 2012),
#
# Based upon the model simulation, we estimated that approximately 18 days after the flood (28 September 2011), nearly
# 6.73 X 10^6 t of fluvial sediment deposited inside the Bay, 0.15?106

#sys.exit()
#vel_mag = np.sqrt(ubar**2 + vbar**2)

# initialize the figure
myFmt = DateFormatter("%m/%d")
dayint = 10
fig, (ax) = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(12, 8))
fig.subplots_adjust(hspace=0.25)

# start plotting
ax[0].plot_date(datetime_list,mud_tot_conc_sr, label='mud_01',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1, color='r')
ax[0].plot_date(datetime_list,sand_tot_conc_sr,label='sand_01',
              xdate=True, linestyle='--', linewidth=0.5,
              marker='', markersize=1, color='r')
ax[0].set_ylabel('%s sed [t]'%trans_name_sr)
ax[0].xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
ax[0].xaxis.set_major_formatter(myFmt)
#ax[0].grid(True)
ax[0].legend(loc='upper left')
ax[2].plot_date(datetime_list,(mud_tot_conc_sr+sand_tot_conc_sr), label=trans_name_sr,
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1, color='r')


# Turkey Point to Sandy Point
trans_name_t2s = 'T2'
x_t2s = np.array(list(range(42,67)))  #
y_t2s = np.array([58]*len(x_t2s))

# verify location
plant_height[x_t2s, y_t2s] = 2


# calculate the total concentration at each time step for the transect (x,y) and depth.
# (h+zeta)/(pm*pn) * mud_01_depth_sum(time,x,y)
# sum over depth

mud_01_trans_t2s = mud_01[:, :, x_t2s, y_t2s]
mud_t2s = np.sum(mud_01_trans_t2s, axis=1)

sand_01_trans_t2s = sand_01[:, :, x_t2s, y_t2s]
sand_t2s = np.sum(sand_01_trans_t2s, axis=1)

# calculate water column volume for points x and y
height_t2s = h[x_t2s, y_t2s]+zeta[:, x_t2s, y_t2s]

#sys.exit()
coord_t2s = pm[x_t2s, y_t2s]*pn[x_t2s, y_t2s]
vol_t2s = height_t2s/coord_t2s
#sys.exit()
mud_trans_conc_t2s = vol_t2s*mud_t2s
mud_tot_conc_t2s = np.sum(mud_trans_conc_t2s, axis=1)/1000  # metric ton

sand_trans_conc_t2s = vol_t2s*sand_t2s
sand_tot_conc_t2s = np.sum(sand_trans_conc_t2s, axis=1)/1000  # metric ton

total_sediment_t2s = (np.sum(sand_tot_conc_t2s)+np.sum(mud_tot_conc_t2s))  # metric ton
print("Transect Turkey Point to Sandy Point total suspended sediment mass = %e tons" % total_sediment_t2s)
# totals to 4.62 x 10^6 metric tons exiting the bay
#
# From Palinkas 2014:
# Various estimates have been generated for the amount of sediment delivered to the Chesapeake Bay associated with
# TS Lee, ranging from 6.7 to 19.0 X 10^6 t(Cheng et al., 2013; Hirsch, 2012),
#
# Based upon the model simulation, we estimated that approximately 18 days after the flood (28 September 2011), nearly
# 6.73 X 10^6 t of fluvial sediment deposited inside the Bay, 0.15?106

#sys.exit()
#vel_mag = np.sqrt(ubar**2 + vbar**2)

myFmt = DateFormatter("%m/%d")
dayint = 10

#fig, (ax) = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(12, 8))
#fig.subplots_adjust(hspace=0.05)
ax[1].plot_date(datetime_list,mud_tot_conc_t2s, label='mud_01',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1, color='b')
ax[1].plot_date(datetime_list,sand_tot_conc_t2s,label='sand_01',
              xdate=True, linestyle='--', linewidth=0.5,
              marker='', markersize=1, color='b')
ax[1].set_ylabel('%s sed [t]'%trans_name_t2s)
ax[1].xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
ax[1].xaxis.set_major_formatter(myFmt)
#ax[1].grid(True)
ax[1].legend(loc='upper left')
#ylim = ax[1].get_ylim()

ax[2].plot_date(datetime_list,(mud_tot_conc_t2s+sand_tot_conc_t2s), label=trans_name_t2s,
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1, color='b')
ax[2].legend(loc='upper left')
ax[2].set_ylabel('sand+mud [t]')
plt.suptitle('Total Suspended Sediment Mass')
plt.figure()
plt.pcolor(f.variables['lon_rho'][:][:], f.variables['lat_rho'][:][:], plant_height)
plt.title('Transect Locations')