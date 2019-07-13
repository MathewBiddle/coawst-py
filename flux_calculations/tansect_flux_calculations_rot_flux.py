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
ocean_time = f.variables['ocean_time'][0:131]
lat = f.variables['lat_rho'][:]
lon = f.variables['lon_rho'][:]
h = f.variables['h'][:]
zeta = f.variables['zeta'][0:131,:,:]
mud_01 = f.variables['mud_01'][0:131,:,:,:]
sand_01 = f.variables['sand_01'][0:131,:,:,:]
ubar = f.variables['ubar_eastward'][0:131,:,:]
vbar = f.variables['vbar_northward'][0:131,:,:]
u = f.variables['u_eastward'][0:131,:,:,:]
v = f.variables['v_northward'][0:131,:,:,:]
pm = f.variables['pm'][:] #XI --> cell width in x dir.
pn = f.variables['pn'][:] #ETA --> cell width in y dir. Want to use this for Surface Area Calcs
s_rho = f.variables['s_rho'][:] # depth levels

tx1 = 1202
tx2 = 1466
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
#t1_angle = []
t1_SSC_flux_ux_rot = np.empty(t1_SSC_flux_yn.shape)
t1_SSC_flux_vy_rot = np.empty(t1_SSC_flux_xe.shape)

# compute angle and rotate
srho_angle = 3 # decide on s-rho coordinate for angle computation
#t1_angle = coawstpy.maj_ax(t1_SSC_flux_xe[tx1:tx2, srho_angle, :], t1_SSC_flux_yn[tx1:tx2, srho_angle, :])
t1_angle = 338
for t in range(0,t1_SSC_flux_yn.shape[0]):  # time
    #t1_angle.append(coawstpy.maj_ax(t1_SSC_flux_xe[t, srho_angle, :], t1_SSC_flux_yn[t, srho_angle, :])) # angle for entire cross-section, pick a depth
    #  ux, uy = coawstpy.rot2xy(ubar_trans[t, :], vbar_trans[t, :], angle)
    for xy in range(0, t1_SSC_flux_yn.shape[2]):  # cell in xy
        for z in range(0,t1_SSC_flux_yn.shape[1]):
            # ux - along channel velocity
            # vy - cross channel velocity (positive to the left of the along channel
            t1_SSC_flux_ux_rot[t, z, xy], t1_SSC_flux_vy_rot[t, z, xy] = coawstpy.rot2xy(t1_SSC_flux_xe[t, z, xy], t1_SSC_flux_yn[t, z, xy], t1_angle)


# Turkey Point to Sandy Point
trans_name = 'T2'
x = np.array(list(range(42,67)))  #
y = np.array([58]*len(x))

# Verify point location
for i in range(len(x)):
    l = i/5
    plant_height[x[i], y[i]] = l

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
flux_face_x = cell_thick/pm[x, y]  # m^2

t2_SSC_flux_yn = np.empty(mud_01_trans.shape)
t2_SSC_flux_xe = np.empty(mud_01_trans.shape)
for i in range(mud_01_trans.shape[1]):  # for each depth level
    # compute flux before rotating axially.
    SSC = mud_01_trans[:, i, :]+sand_01_trans[:, i, :] # kg m^-3
    t2_SSC_flux_yn[:, i, :] = t2_v_trans[:, i, :]*flux_face_x*SSC # (m s^-1) * (m^2) * (kg m^-3) = kg s^-1
    t2_SSC_flux_xe[:, i, :] = t2_u_trans[:, i, :]*flux_face_y*SSC

# variable initialize
#t2_angle = []
t2_SSC_flux_ux_rot = np.empty(t2_SSC_flux_yn.shape)
t2_SSC_flux_vy_rot = np.empty(t2_SSC_flux_xe.shape)
#t2_angle = coawstpy.maj_ax(t2_SSC_flux_xe[tx1:tx2, srho_angle, :], t2_SSC_flux_yn[tx1:tx2, srho_angle, :])
t2_angle = 358
# compute angle and rotate
for t in range(0,t2_SSC_flux_yn.shape[0]):  # time
    #t2_angle.append(coawstpy.maj_ax(t2_SSC_flux_xe[t, srho_angle, :], t2_SSC_flux_yn[t, srho_angle, :])) # angle for entire cross-section, pick a depth
    #  ux, uy = coawstpy.rot2xy(ubar_trans[t, :], vbar_trans[t, :], angle)
    for xy in range(0, t2_SSC_flux_yn.shape[2]):  # cell in xy
        for z in range(0,t2_SSC_flux_yn.shape[1]):
            # ux - along channel velocity
            # vy - cross channel velocity (positive to the left of the along channel
            t2_SSC_flux_ux_rot[t, z, xy], t2_SSC_flux_vy_rot[t, z, xy] = coawstpy.rot2xy(t2_SSC_flux_xe[t, z, xy], t2_SSC_flux_yn[t, z, xy], t2_angle)

# magnitude of the rotated flux
t1_mag_ssc = np.sqrt(t1_SSC_flux_ux_rot**2 + t1_SSC_flux_vy_rot**2)
t2_mag_ssc = np.sqrt(t2_SSC_flux_ux_rot**2 + t2_SSC_flux_vy_rot**2)

# sum of magnitude of rotated flux across transect and depth
t1_SSC_ts_sum = np.sum(t1_mag_ssc, axis=(1,2)) ##TODO review this. use t1_mag_ssc?
t2_SSC_ts_sum = np.sum(t2_mag_ssc, axis=(1,2))

# cumulative sum
t1_cumsum_ssc = np.nancumsum(t1_SSC_ts_sum, axis=0)
t1_total_sed = t1_cumsum_ssc[-1]*3600 # sum is integrated over every hour, so we multiply by 3600 seconds
#(datetime_list[-1]-datetime_list[1]).total_seconds() #kg

t2_cumsum_ssc = np.nancumsum(t2_SSC_ts_sum, axis=0)
t2_total_sed = t2_cumsum_ssc[-1]*3600# (datetime_list[-1]-datetime_list[1]).total_seconds()

print('Total Sediment across transect over %s seconds\nT1 = %e kg = %e tons\nT2 = %e kg = %e tons' %
      ((datetime_list[-1]-datetime_list[1]).total_seconds(), t1_total_sed, t1_total_sed/1000, t2_total_sed, t2_total_sed/1000))
sys.exit()
## Create some plots
# pick a depth and cell for investigation
# srho = depth coordinate
# c_* cell number along transect (t1 and t2).
srho = 3
c_t1 = 1
c_t2 = 12


# map of transect
plt.figure()
plt.pcolor(f.variables['lon_rho'][:][:], f.variables['lat_rho'][:][:], plant_height)
plt.title('Transects')

#  plot raw velocity
fig, (axd) = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
fig.subplots_adjust(hspace=0.25)
axd[0].plot_date(datetime_list, t1_v_trans[:, srho, c_t1], label='vn',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axd[0].plot_date(datetime_list, t1_u_trans[:, srho, c_t1], label='ue',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axd[0].set_ylabel('velocity (m/s)')
axd[0].set_title('Raw velocity at s-rho=%i, cell=%i, Transect T1' % (srho, c_t1))


axd[1].plot_date(datetime_list, t2_v_trans[:, srho, c_t2], label='vn',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axd[1].plot_date(datetime_list, t2_u_trans[:, srho, c_t2], label='ue',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axd[1].set_ylabel('velocity (m/s)')
axd[1].set_title('Raw velocity at s-rho=%i, cell=%i, Transect T2' % (srho, c_t2))
axd[1].xaxis.set_major_locator(mdates.DayLocator(interval=30))
axd[1].xaxis.set_major_formatter(DateFormatter("%m/%d"))
axd[1].legend()


# plot non-rotated flux
srho = 3
fig, (axa) = plt.subplots(nrows=2, ncols=1, figsize=(12, 8),sharex=True)
fig.subplots_adjust(hspace=0.25)
axa[0].plot_date(datetime_list, t1_SSC_flux_yn[:, srho, c_t1], label='yn',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axa[0].plot_date(datetime_list, t1_SSC_flux_xe[:, srho, c_t1], label='xe',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axa[0].set_ylabel('Total SSC flux (kg/s)')
axa[0].set_title('SSC flux N-E, non-rotated at s-rho=%i, cell=%i, Transect T1' % (srho, c_t1))

axa[1].plot_date(datetime_list, t2_SSC_flux_yn[:, srho, c_t2], label='yn',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axa[1].plot_date(datetime_list, t2_SSC_flux_xe[:, srho, c_t2], label='xe',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axa[1].set_ylabel('Total SSC flux (kg/s)')
axa[1].set_title('SSC flux N-E, non-rotated at s-rho=%i, cell=%i, Transect T2' % (srho, c_t2))
axa[1].xaxis.set_major_locator(mdates.DayLocator(interval=30))
axa[1].xaxis.set_major_formatter(DateFormatter("%m/%d"))
axa[1].legend()


# plot angle calculation
# fig, (axb) = plt.subplots(nrows=2, ncols=1, figsize=(12, 8),sharex=True)
# axb[0].plot_date(datetime_list, t1_angle, label='T1',
#                xdate=True, linestyle='', linewidth=1,
#                marker='.', markersize=1)
# axb[0].set_title('Calculated angle from SSC flux N-E at Transect T1 from s-rho=%i'% srho_angle)
# axb[0].set_ylabel('Angle')
# axb[1].plot_date(datetime_list, t2_angle, label='T2',
#               xdate=True, linestyle='', linewidth=1,
#               marker='.', markersize=1)
# axb[1].set_title('Calculated angle from SSC flux N-E at Transect T2')
# axb[1].set_ylabel('Angle')
# axb[1].xaxis.set_major_locator(mdates.DayLocator(interval=30))
# axb[1].xaxis.set_major_formatter(DateFormatter("%m/%d"))

# plot rotated flux values at a specific depth and cell
fig, (axc) = plt.subplots(nrows=2, ncols=1, figsize=(12, 8),sharex=True)
fig.subplots_adjust(hspace=0.25)

axc[0].plot_date(datetime_list, t1_SSC_flux_ux_rot[:, srho, c_t1], label='ux',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axc[0].plot_date(datetime_list,t1_SSC_flux_vy_rot[:, srho, c_t1],label='vy',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axc[0].set_title('Along (ux) and Cross (vy) channel total SSC flux at s-rho=%i, cell=%i, Transect T1' % (srho,c_t1))
axc[0].set_ylabel('Total SSC flux (kg/s)')

# plot along and cross channel fluxes at a depth and cell
axc[1].plot_date(datetime_list, t2_SSC_flux_ux_rot[:, srho, c_t2], label='ux',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axc[1].plot_date(datetime_list,t2_SSC_flux_vy_rot[:, srho, c_t2],label='vy',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axc[1].set_title('Along (ux) and Cross (vy) channel total SSC flux at s-rho=%i, cell=%i, Transect T2' % (srho,c_t2))
axc[1].set_ylabel('Total SSC flux (kg/s)')
axc[1].xaxis.set_major_locator(mdates.DayLocator(interval=30))
axc[1].xaxis.set_major_formatter(DateFormatter("%m/%d"))
axc[1].legend()


## Plot magnitude of rotated flux
fig, (axe) = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))
axe.plot_date(datetime_list, t1_mag_ssc[:, srho, c_t1], label='T1 cell %i' % c_t1,
                 xdate=True, linestyle='', linewidth=1,
                 marker='.', markersize=1)
axe.plot_date(datetime_list, t2_mag_ssc[:, srho, c_t2], label='T2 cell %i' % c_t2,
                 xdate=True, linestyle='', linewidth=1,
                 marker='.', markersize=1)
axe.xaxis.set_major_locator(mdates.DayLocator(interval=30))
axe.xaxis.set_major_formatter(DateFormatter("%m/%d"))
axe.legend()
axe.set_title('Magnitude of the rotated flux at s-rho=%i' % srho)
axe.set_ylabel('SSC flux (kg/s)')


# plot the sum of the magnitude of the rotated flux for each transect
fig, (axf) = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))
axf.plot_date(datetime_list, t1_SSC_ts_sum, label='T1',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axf.plot_date(datetime_list, t2_SSC_ts_sum, label='T2',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axf.xaxis.set_major_locator(mdates.DayLocator(interval=30))
axf.xaxis.set_major_formatter(DateFormatter("%m/%d"))
axf.legend()
axf.set_title('Sum (depth and transect) of the magnitude of the rotated flux for each transect')
axf.set_ylabel('SSC flux (kg/s)')


# plot cumulative sum of SSC for transect
fig, (axg) = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))
axg.plot_date(datetime_list, t1_cumsum_ssc, label='T1',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axg.text(datetime_list[-1],t1_cumsum_ssc[-1],'%.2e tons' % (t1_total_sed/1000))
axg.plot_date(datetime_list, t2_cumsum_ssc, label='T2',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axg.text(datetime_list[-1],t2_cumsum_ssc[-1],'%.2e tons' % (t2_total_sed/1000))
axg.xaxis.set_major_locator(mdates.DayLocator(interval=30))
axg.xaxis.set_major_formatter(DateFormatter("%m/%d"))
axg.legend()
axg.set_title('Cumulative sum of the rotated flux for transects')
axg.set_ylabel('Cumulative sum of SSC flux (kg/s)')
