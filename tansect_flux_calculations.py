import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import numpy as np
import netCDF4
import coawstpy

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

# identify the points to extract data
# Susquehanna River mouth
trans_name = 'T1'
x = np.array([29,30,31,32])
y = np.array([13,12,11,10])

# Turkey Point to Sandy Point
#trans_name = 'T2'
#x = np.array(list(range(42,67)))  #
#y = np.array([58]*len(x))

# Verify point location
plant_height = f.variables['plant_height'][0, 0, :, :]
plant_height = np.ma.masked_greater(plant_height, 1)
for i in range(len(x)):
    plant_height[x[i], y[i]] = i

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
mud_01_trans = mud_01[:, :, x, y]  # kg/m3
sand_01_trans = sand_01[:, :, x, y]  # kg/m3
ubar_trans = ubar[:, x, y]  # m/s
vbar_trans = vbar[:, x, y]  # m/s
t1_v_trans = v[:, :, x, y]  # m/s
t1_u_trans = u[:, :, x, y]  # m/s

water_col_depth = h[x, y]+zeta[:, x, y]  # m
cell_thick = water_col_depth/s_rho.shape[0]  # m

flux_face_y = cell_thick/pn[x, y]  # m^2
flux_face_x = cell_thick/pm[x, y]  # m^2

t1_SSC_flux_yn = np.empty(mud_01_trans.shape)
t1_SSC_flux_xe = np.empty(mud_01_trans.shape)
for i in range(mud_01_trans.shape[1]):  # for each depth level
    # compute flux before rotating axially.
    SSC = mud_01_trans[:, i, :]+sand_01_trans[:, i, :] # kg m^-3
    t1_SSC_flux_yn[:, i, :] = t1_v_trans[:, i, :]*flux_face_x*SSC # (m s^-1) * (m^2) * (kg m^-3) = kg s^-1
    t1_SSC_flux_xe[:, i, :] = t1_u_trans[:, i, :]*flux_face_y*SSC

# variable initialize
angle = []
t1_SSC_flux_ux_rot = np.empty(t1_SSC_flux_yn.shape)
t1_SSC_flux_vy_rot = np.empty(t1_SSC_flux_xe.shape)

# compute angle and rotate
for t in range(0,t1_SSC_flux_yn.shape[0]):  # time
    angle.append(coawstpy.maj_ax(t1_SSC_flux_xe[t, 2, :], t1_SSC_flux_yn[t, 2, :])) # angle for entire cross-section, pick a depth
    #  ux, uy = coawstpy.rot2xy(ubar_trans[t, :], vbar_trans[t, :], angle)
    for xy in range(0, t1_SSC_flux_yn.shape[2]):  # cell in xy
        for z in range(0,t1_SSC_flux_yn.shape[1]):
            # ux - along channel velocity
            # vy - cross channel velocity (positive to the left of the along channel
            t1_SSC_flux_ux_rot[t, z, xy], t1_SSC_flux_vy_rot[t, z, xy] = coawstpy.rot2xy(t1_SSC_flux_xe[t, z, xy], t1_SSC_flux_yn[t, z, xy], angle[t])





# Turkey Point to Sandy Point
trans_name = 'T2'
x = np.array(list(range(42,67)))  #
y = np.array([58]*len(x))

# Verify point location
plant_height = f.variables['plant_height'][0, 0, :, :]
plant_height = np.ma.masked_greater(plant_height, 1)
for i in range(len(x)):
    plant_height[x[i], y[i]] = i

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
mud_01_trans = mud_01[:, :, x, y]  # kg/m3
sand_01_trans = sand_01[:, :, x, y]  # kg/m3
ubar_trans = ubar[:, x, y]  # m/s
vbar_trans = vbar[:, x, y]  # m/s
t2_v_trans = v[:, :, x, y]  # m/s
t2_u_trans = u[:, :, x, y]  # m/s

water_col_depth = h[x, y]+zeta[:, x, y]  # m
cell_thick = water_col_depth/s_rho.shape[0]  # m

flux_face_y = cell_thick/pn[x, y]  # m^2
flux_face_x = cell_thick/pm[x,y]  # m^2

t2_SSC_flux_yn = np.empty(mud_01_trans.shape)
t2_SSC_flux_xe = np.empty(mud_01_trans.shape)
for i in range(mud_01_trans.shape[1]):  # for each depth level
    # compute flux before rotating axially.
    SSC = mud_01_trans[:, i, :]+sand_01_trans[:, i, :] # kg m^-3
    t2_SSC_flux_yn[:, i, :] = t2_v_trans[:, i, :]*flux_face_x*SSC # (m s^-1) * (m^2) * (kg m^-3) = kg s^-1
    t2_SSC_flux_xe[:, i, :] = t2_u_trans[:, i, :]*flux_face_y*SSC

# variable initialize
angle = []
t2_SSC_flux_ux_rot = np.empty(t2_SSC_flux_yn.shape)
t2_SSC_flux_vy_rot = np.empty(t2_SSC_flux_xe.shape)

# compute angle and rotate
for t in range(0,t2_SSC_flux_yn.shape[0]):  # time
    angle.append(coawstpy.maj_ax(t2_SSC_flux_xe[t, 2, :], t2_SSC_flux_yn[t, 2, :])) # angle for entire cross-section, pick a depth
    #  ux, uy = coawstpy.rot2xy(ubar_trans[t, :], vbar_trans[t, :], angle)
    for xy in range(0, t2_SSC_flux_yn.shape[2]):  # cell in xy
        for z in range(0,t2_SSC_flux_yn.shape[1]):
            # ux - along channel velocity
            # vy - cross channel velocity (positive to the left of the along channel
            t2_SSC_flux_ux_rot[t, z, xy], t2_SSC_flux_vy_rot[t, z, xy] = coawstpy.rot2xy(t2_SSC_flux_xe[t, z, xy], t2_SSC_flux_yn[t, z, xy], angle[t])


## Create some plots
# map of transect
plt.figure()
plt.pcolor(f.variables['lon_rho'][:][:], f.variables['lat_rho'][:][:], plant_height)
plt.title('%s Transect' % trans_name)


# plot raw velocity
srho = 3
c = 1
fig, (axd) = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
fig.subplots_adjust(hspace=0.25)
axd[0].plot_date(datetime_list, t1_v_trans[:, srho, c], label='vn',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axd[0].plot_date(datetime_list, t1_u_trans[:, srho, c], label='ue',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axd[0].set_ylabel('velocity (m/s)')
axd[0].set_title('Raw velocity at s-rho=%i, cell=%i, Transect T1' % (srho, c))

c = 12
axd[1].plot_date(datetime_list, t2_v_trans[:, srho, c], label='vn',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axd[1].plot_date(datetime_list, t2_u_trans[:, srho, c], label='ue',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axd[1].set_ylabel('velocity (m/s)')
axd[1].set_title('Raw velocity at s-rho=%i, cell=%i, Transect T2' % (srho, c))
axd[1].xaxis.set_major_locator(mdates.DayLocator(interval=30))
axd[1].xaxis.set_major_formatter(DateFormatter("%m/%d"))
axd[1].legend()


# plot non-rotated flux
srho = 3
c = 1
fig, (axa) = plt.subplots(nrows=2, ncols=1, figsize=(12, 8),sharex=True)
fig.subplots_adjust(hspace=0.25)
axa[0].plot_date(datetime_list, t1_SSC_flux_yn[:, srho, c], label='yn',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axa[0].plot_date(datetime_list, t1_SSC_flux_xe[:, srho, c], label='xe',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axa[0].set_ylabel('Total SSC flux (kg/s)')
axa[0].set_title('SSC flux N-E, non-rotated at s-rho=%i, cell=%i, Transect T1' % (srho, c))

c = 12
axa[1].plot_date(datetime_list, t2_SSC_flux_yn[:, srho, c], label='yn',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axa[1].plot_date(datetime_list, t2_SSC_flux_xe[:, srho, c], label='xe',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axa[1].set_ylabel('Total SSC flux (kg/s)')
axa[1].set_title('SSC flux N-E, non-rotated at s-rho=%i, cell=%i, Transect T2' % (srho, c))
axa[1].xaxis.set_major_locator(mdates.DayLocator(interval=30))
axa[1].xaxis.set_major_formatter(DateFormatter("%m/%d"))
axa[1].legend()


# plot angle calculation
# fig, (axb) = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))
# axb.plot_date(datetime_list, angle,
#               xdate=True, linestyle='', linewidth=1,
#               marker='.', markersize=1)
# axb.set_title('Calculated angle from SSC flux N-E at %s' % trans_name)
# axb.set_ylabel('Angle')
# axb.xaxis.set_major_locator(mdates.DayLocator(interval=30))
# axb.xaxis.set_major_formatter(DateFormatter("%m/%d"))

# plot rotated flux values
fig, (axc) = plt.subplots(nrows=2, ncols=1, figsize=(12, 8),sharex=True)
fig.subplots_adjust(hspace=0.25)

c = 1
axc[0].plot_date(datetime_list, t1_SSC_flux_ux_rot[:, srho, c], label='ux',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axc[0].plot_date(datetime_list,t1_SSC_flux_vy_rot[:, srho, c],label='vy',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axc[0].set_title('Along (ux) and Cross (vy) channel total SSC flux at T1, s-rho=%i, cell=%i,' % (srho,c))
axc[0].set_ylabel('Total SSC flux (kg/s)')


c = 12
axc[1].plot_date(datetime_list, t2_SSC_flux_ux_rot[:, srho, c], label='ux',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axc[1].plot_date(datetime_list,t2_SSC_flux_vy_rot[:, srho, c],label='vy',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axc[1].set_title('Along (ux) and Cross (vy) channel total SSC flux at T2, s-rho=%i, cell=%i,' % (srho,c))
axc[1].set_ylabel('Total SSC flux (kg/s)')

axc[1].xaxis.set_major_locator(mdates.DayLocator(interval=30))
axc[1].xaxis.set_major_formatter(DateFormatter("%m/%d"))
axc[1].legend()

