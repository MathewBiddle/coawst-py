import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import numpy as np
import netCDF4
import coawstpy
import scipy.integrate as integrate

# bring in the data
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'
inputfile = dir+'/upper_ches_avg.nc'
f = netCDF4.Dataset(inputfile, 'r')

# Get the data we want
ocean_time = f.variables['ocean_time'][:]
lat = f.variables['lat_rho'][:]
lon = f.variables['lon_rho'][:]
Huon_sand_01 = f.variables['Huon_sand_01'][:] # east-west (from panoply should be S)
Hvom_sand_01 = f.variables['Hvom_sand_01'][:] # north-south (from panoply should be E)
Huon_mud_01 = f.variables['Huon_mud_01'][:] # east-west
Hvom_mud_01 = f.variables['Hvom_mud_01'][:] # north-south
plant_height = f.variables['Huon_sand_01'][0,0,:,:]
s_rho = f.variables['s_rho'][:] # depth levels

tx1 = 1202
tx2 = 1466
srho_angle = 3
datetime_list=[]
for sec in ocean_time:
    datetime_list.append(
        netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))



###############
# SR entrance #
###############
trans_name = 'T0'
x = np.array([27,26,25,24,23])
#x = np.array(list(range(20,28)))
#y = np.array([0]*len(x))
y = np.array([0,1,2,3,4])
# Verify point location
for i in range(len(x)):
    l = i/5
    plant_height[x[i], y[i]] = 10#l*100

# Gather subset data
t0_mud_01_flux_xe = Hvom_mud_01[:, :, x, y]
t0_mud_01_flux_yn = Huon_mud_01[:, :, x, y]
t0_sand_01_flux_xe = Hvom_sand_01[:, :, x, y]
t0_sand_01_flux_yn = Huon_sand_01[:, :, x, y]

# init arrays
t0_mud_01_flux_ux_rot = np.empty(t0_mud_01_flux_xe.shape)
t0_mud_01_flux_vy_rot = np.empty(t0_mud_01_flux_yn.shape)
t0_sand_01_flux_ux_rot = np.empty(t0_sand_01_flux_xe.shape)
t0_sand_01_flux_vy_rot = np.empty(t0_sand_01_flux_yn.shape)

# compute angle and rotate
# angle = maj_ax(ue, vn)
t0_angle = coawstpy.maj_ax(t0_mud_01_flux_xe[tx1:tx2, srho_angle, :], t0_mud_01_flux_yn[tx1:tx2, srho_angle, :])
t0_angle_list = []
for t in range(0,t0_mud_01_flux_yn.shape[0]):  # time
    t0_angle_list.append(coawstpy.maj_ax(t0_mud_01_flux_xe[t, srho_angle, :], t0_mud_01_flux_yn[t, srho_angle, :]))
    for xy in range(0, t0_mud_01_flux_yn.shape[2]):  # cell in xy
        for z in range(0,t0_mud_01_flux_yn.shape[1]):
            # ux,uy = rot2xy(ue, vn, projangle)
            # ux - along channel velocity
            # vy - cross channel velocity (positive to the left of the along channel
            t0_mud_01_flux_ux_rot[t, z, xy], t0_mud_01_flux_vy_rot[t, z, xy] = coawstpy.rot2xy(
                t0_mud_01_flux_xe[t, z, xy], t0_mud_01_flux_yn[t, z, xy], t0_angle)

            t0_sand_01_flux_ux_rot[t, z, xy], t0_sand_01_flux_vy_rot[t, z, xy] = coawstpy.rot2xy(
                t0_sand_01_flux_xe[t, z, xy], t0_sand_01_flux_yn[t, z, xy], t0_angle)

t0_mag_mud_01 = np.sqrt(t0_mud_01_flux_ux_rot**2 + t0_mud_01_flux_vy_rot**2)
t0_mag_sand_01 = np.sqrt(t0_sand_01_flux_ux_rot**2 + t0_sand_01_flux_vy_rot**2)
t0_mag_ssc = t0_mag_mud_01 + t0_mag_sand_01 # total sediment together
t0_SSC_ts_sum = np.sum(t0_mag_ssc, axis=(1,2)) # sum over depth and transect
t0_cumsum_ssc = np.nancumsum(t0_SSC_ts_sum, axis=0)  # cumulative sum over time
t0_cumtrapz_ssc = integrate.cumtrapz(t0_SSC_ts_sum, axis=0)  # cumulative integral over time
t0_total_sed = t0_cumtrapz_ssc[-1]*3600 # sum is integrated for every hour, we multiply by 3600 seconds = total weight


###########################
# Susquehanna River mouth #
###########################
trans_name = 'T1'
#x = np.array([29,30,31,32])
#y = np.array([13,12,11,10])
x = np.array([29,30,31])
y = np.array([13,12,11])

# Verify point location
for i in range(len(x)):
    plant_height[x[i], y[i]] = 10

# Gather subset data
t1_mud_01_flux_xe = Hvom_mud_01[:, :, x, y]
t1_mud_01_flux_yn = Huon_mud_01[:, :, x, y]
t1_sand_01_flux_xe = Hvom_sand_01[:, :, x, y]
t1_sand_01_flux_yn = Huon_sand_01[:, :, x, y]

t1_mud_01_flux_ux_rot = np.empty(t1_mud_01_flux_xe.shape)
t1_mud_01_flux_vy_rot = np.empty(t1_mud_01_flux_yn.shape)
t1_sand_01_flux_ux_rot = np.empty(t1_sand_01_flux_xe.shape)
t1_sand_01_flux_vy_rot = np.empty(t1_sand_01_flux_yn.shape)

# compute angle and rotate
t1_angle = coawstpy.maj_ax(t1_mud_01_flux_xe[tx1:tx2, srho_angle, :], t1_mud_01_flux_yn[tx1:tx2, srho_angle, :])
t1_angle_list = []
for t in range(0,t1_mud_01_flux_yn.shape[0]):  # time
    t1_angle_list.append(coawstpy.maj_ax(t1_mud_01_flux_xe[t, srho_angle, :], t1_mud_01_flux_yn[t, srho_angle, :]))
    for xy in range(0, t1_mud_01_flux_yn.shape[2]):  # cell in xy
        for z in range(0,t1_mud_01_flux_yn.shape[1]):
            # ux - along channel velocity
            # vy - cross channel velocity (positive to the left of the along channel
            t1_mud_01_flux_ux_rot[t, z, xy], t1_mud_01_flux_vy_rot[t, z, xy] = coawstpy.rot2xy(
                t1_mud_01_flux_xe[t, z, xy], t1_mud_01_flux_yn[t, z, xy], t1_angle)

            t1_sand_01_flux_ux_rot[t, z, xy], t1_sand_01_flux_vy_rot[t, z, xy] = coawstpy.rot2xy(
                t1_sand_01_flux_xe[t, z, xy], t1_sand_01_flux_yn[t, z, xy], t1_angle)

# magnitude of the rotated flux
t1_mag_mud_01 = np.sqrt(t1_mud_01_flux_ux_rot**2 + t1_mud_01_flux_vy_rot**2)
t1_mag_sand_01 = np.sqrt(t1_sand_01_flux_ux_rot**2 + t1_sand_01_flux_vy_rot**2)
t1_mag_ssc = t1_mag_mud_01 + t1_mag_sand_01 # total sediment together
t1_SSC_ts_sum = np.sum(t1_mag_ssc, axis=(1,2)) # sum over depth and transect
t1_cumsum_ssc = np.nancumsum(t1_SSC_ts_sum, axis=0)  # cumulative sum over time
t1_cumtrapz_ssc = integrate.cumtrapz(t1_SSC_ts_sum, axis=0)  # cumulative integral over time
t1_total_sed = t1_cumtrapz_ssc[-1]*3600 # sum is integrated for every hour, we multiply by 3600 seconds = total weight


###############################
# Turkey Point to Sandy Point #
###############################
trans_name = 'T2'
x = np.array(list(range(42,66)))  #
y = np.array([58]*len(x))

# Verify point location
for i in range(len(x)):
    l = i/5
    plant_height[x[i], y[i]] = 10#l*100

# Gather subset data
t2_mud_01_flux_xe = Hvom_mud_01[:, :, x, y]
t2_mud_01_flux_yn = Huon_mud_01[:, :, x, y]
t2_sand_01_flux_xe = Hvom_sand_01[:, :, x, y]
t2_sand_01_flux_yn = Huon_sand_01[:, :, x, y]

# init arrays
t2_mud_01_flux_ux_rot = np.empty(t2_mud_01_flux_xe.shape)
t2_mud_01_flux_vy_rot = np.empty(t2_mud_01_flux_yn.shape)
t2_sand_01_flux_ux_rot = np.empty(t2_sand_01_flux_xe.shape)
t2_sand_01_flux_vy_rot = np.empty(t2_sand_01_flux_yn.shape)

# compute angle and rotate
t2_angle = coawstpy.maj_ax(t2_mud_01_flux_xe[tx1:tx2, srho_angle, :], t2_mud_01_flux_yn[tx1:tx2, srho_angle, :])
t2_angle_list = []
for t in range(0,t2_mud_01_flux_yn.shape[0]):  # time
    t2_angle_list.append(coawstpy.maj_ax(t2_mud_01_flux_xe[t, srho_angle, :], t2_mud_01_flux_yn[t, srho_angle, :]))
    for xy in range(0, t2_mud_01_flux_yn.shape[2]):  # cell in xy
        for z in range(0,t2_mud_01_flux_yn.shape[1]):
            # ux - along channel velocity
            # vy - cross channel velocity (positive to the left of the along channel
            t2_mud_01_flux_ux_rot[t, z, xy], t2_mud_01_flux_vy_rot[t, z, xy] = coawstpy.rot2xy(
                t2_mud_01_flux_xe[t, z, xy], t2_mud_01_flux_yn[t, z, xy], t2_angle)

            t2_sand_01_flux_ux_rot[t, z, xy], t2_sand_01_flux_vy_rot[t, z, xy] = coawstpy.rot2xy(
                t2_sand_01_flux_xe[t, z, xy], t2_sand_01_flux_yn[t, z, xy], t2_angle)

# magnitude of rotated flux
t2_mag_mud_01 = np.sqrt(t2_mud_01_flux_ux_rot**2 + t2_mud_01_flux_vy_rot**2)
t2_mag_sand_01 = np.sqrt(t2_sand_01_flux_ux_rot**2 + t2_sand_01_flux_vy_rot**2)
t2_mag_ssc = t2_mag_mud_01 + t2_mag_sand_01 # total sediment together
t2_SSC_ts_sum = np.sum(t2_mag_ssc, axis=(1,2)) # sum across transect and depth
t2_cumsum_ssc = np.nancumsum(t2_SSC_ts_sum, axis=0)  # cumulative sum over time
t2_cumtrapz_ssc = integrate.cumtrapz(t2_SSC_ts_sum, axis=0)  # cumulative integral over time
t2_total_sed = t2_cumtrapz_ssc[-1]*3600



## print out some information
print('Total Sediment across transect over %s seconds' % (datetime_list[-1]-datetime_list[1]).total_seconds())
print('T0 = %e kg = %e tons' % (t0_total_sed, t0_total_sed/1000))
print('T1 = %e kg = %e tons' % (t1_total_sed, t1_total_sed/1000))
print('T2 = %e kg = %e tons' % (t2_total_sed, t2_total_sed/1000))

## Create some plots
# pick a depth and cell for investigation
# srho = depth coordinate
# c_* cell number along transect (t1 and t2).
srho = 3
c_t0 = 3
c_t1 = 1
c_t2 = 12

# map of transect
plt.figure()
plt.pcolor(lon, lat, plant_height)
plt.title('Transects')

#  plot raw velocity
# fig, (axd) = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
# fig.subplots_adjust(hspace=0.25)
# axd[0].plot_date(datetime_list, t1_v_trans[:, srho, c_t1], label='vn',
#               xdate=True, linestyle='', linewidth=1,
#               marker='.', markersize=1)
# axd[0].plot_date(datetime_list, t1_u_trans[:, srho, c_t1], label='ue',
#               xdate=True, linestyle='', linewidth=1,
#               marker='.', markersize=1)
# axd[0].set_ylabel('velocity (m/s)')
# axd[0].set_title('Raw velocity at s-rho=%i, cell=%i, Transect T1' % (srho, c_t1))
#
#
# axd[1].plot_date(datetime_list, t2_v_trans[:, srho, c_t2], label='vn',
#               xdate=True, linestyle='', linewidth=1,
#               marker='.', markersize=1)
# axd[1].plot_date(datetime_list, t2_u_trans[:, srho, c_t2], label='ue',
#               xdate=True, linestyle='', linewidth=1,
#               marker='.', markersize=1)
# axd[1].set_ylabel('velocity (m/s)')
# axd[1].set_title('Raw velocity at s-rho=%i, cell=%i, Transect T2' % (srho, c_t2))
# axd[1].xaxis.set_major_locator(mdates.DayLocator(interval=30))
# axd[1].xaxis.set_major_formatter(DateFormatter("%m/%d"))
# axd[1].legend()


# plot non-rotated flux
srho = 3
fig, (axa) = plt.subplots(nrows=3, ncols=1, figsize=(12, 8),sharex=True)
fig.subplots_adjust(hspace=0.25)
axa[0].plot_date(datetime_list, t0_mud_01_flux_yn[:, srho, c_t0], label='yn',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axa[0].plot_date(datetime_list, t0_mud_01_flux_xe[:, srho, c_t0], label='xe',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axa[0].set_ylabel('mud_01 flux (kg/s)')
axa[0].set_title('mud_01 flux N-E, non-rotated at s-rho=%i, cell=%i, Transect T0' % (srho, c_t0))

axa[1].plot_date(datetime_list, t1_mud_01_flux_yn[:, srho, c_t1], label='yn',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axa[1].plot_date(datetime_list, t1_mud_01_flux_xe[:, srho, c_t1], label='xe',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axa[1].set_ylabel('mud_01 flux (kg/s)')
axa[1].set_title('mud_01 flux N-E, non-rotated at s-rho=%i, cell=%i, Transect T1' % (srho, c_t1))

axa[2].plot_date(datetime_list, t2_mud_01_flux_yn[:, srho, c_t2], label='yn',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axa[2].plot_date(datetime_list, t2_mud_01_flux_xe[:, srho, c_t2], label='xe',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axa[2].set_ylabel('mud_01 flux (kg/s)')
axa[2].set_title('mud_01 flux N-E, non-rotated at s-rho=%i, cell=%i, Transect T2' % (srho, c_t2))
axa[2].xaxis.set_major_locator(mdates.DayLocator(interval=30))
axa[2].xaxis.set_major_formatter(DateFormatter("%m/%d"))
axa[0].legend()

# plot angle calculation
fig, (axb) = plt.subplots(nrows=1, ncols=1, figsize=(12, 8),sharex=True)
axb.plot_date(datetime_list, t0_angle_list, label='T0',
               xdate=True, linestyle='', linewidth=1,
               marker='.', markersize=1)
axb.plot_date(datetime_list, t1_angle_list, label='T1',
               xdate=True, linestyle='', linewidth=1,
               marker='.', markersize=1)
axb.plot_date(datetime_list, t2_angle_list, label='T2',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1)
axb.set_title('Calculated angle from SSC flux N-E at each transect')
axb.set_ylabel('Angle')
axb.xaxis.set_major_locator(mdates.DayLocator(interval=30))
axb.xaxis.set_major_formatter(DateFormatter("%m/%d"))
axb.legend()

# plot rotated flux values at a specific depth and cell
fig, (axc) = plt.subplots(nrows=3, ncols=1, figsize=(12, 8),sharex=True)
fig.subplots_adjust(hspace=0.25)

axc[0].plot_date(datetime_list, t0_mud_01_flux_ux_rot[:, srho, c_t0], label='ux',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axc[0].plot_date(datetime_list,t0_mud_01_flux_vy_rot[:, srho, c_t0],label='vy',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axc[0].set_title('Along (ux) and Cross (vy) channel mud_01 flux at s-rho=%i, cell=%i, Transect T0, angle = %i' %
                 (srho,c_t0,t0_angle))
axc[0].set_ylabel('mud_01 flux (kg/s)')

axc[1].plot_date(datetime_list, t1_mud_01_flux_ux_rot[:, srho, c_t1], label='ux',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axc[1].plot_date(datetime_list,t1_mud_01_flux_vy_rot[:, srho, c_t1],label='vy',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axc[1].set_title('Along (ux) and Cross (vy) channel mud_01 flux at s-rho=%i, cell=%i, Transect T1, angle = %i' %
                 (srho,c_t1,t1_angle))
axc[1].set_ylabel('mud_01 flux (kg/s)')

axc[2].plot_date(datetime_list, t2_mud_01_flux_ux_rot[:, srho, c_t2], label='ux',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axc[2].plot_date(datetime_list,t2_mud_01_flux_vy_rot[:, srho, c_t2],label='vy',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axc[2].set_title('Along (ux) and Cross (vy) channel mud_01 flux at s-rho=%i, cell=%i, Transect T2, angle = %i' %
                 (srho,c_t2, t2_angle))
axc[2].set_ylabel('mud_01 flux (kg/s)')
axc[2].xaxis.set_major_locator(mdates.DayLocator(interval=30))
axc[2].xaxis.set_major_formatter(DateFormatter("%m/%d"))
axc[0].legend()


## Plot magnitude of rotated flux
fig, (axe) = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))
axe.plot_date(datetime_list, t0_mag_ssc[:, srho, c_t0], label='T0 cell %i' % c_t0,
                 xdate=True, linestyle='-', linewidth=0.5,
                 marker='', markersize=1)
axe.plot_date(datetime_list, t1_mag_ssc[:, srho, c_t1], label='T1 cell %i' % c_t1,
                 xdate=True, linestyle='-', linewidth=0.5,
                 marker='', markersize=1)
axe.plot_date(datetime_list, t2_mag_ssc[:, srho, c_t2], label='T2 cell %i' % c_t2,
                 xdate=True, linestyle='-', linewidth=0.5,
                 marker='', markersize=1)
axe.xaxis.set_major_locator(mdates.DayLocator(interval=30))
axe.xaxis.set_major_formatter(DateFormatter("%m/%d"))
axe.legend()
axe.set_title('Magnitude of the rotated total SSC flux at s-rho=%i' % srho)
axe.set_ylabel('SSC flux (kg/s)')


# plot the sum of the magnitude of the rotated flux for each transect
fig, (axf) = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))
axf.plot_date(datetime_list, t0_SSC_ts_sum, label='T0',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axf.plot_date(datetime_list, t1_SSC_ts_sum, label='T1',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axf.plot_date(datetime_list, t2_SSC_ts_sum, label='T2',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axf.xaxis.set_major_locator(mdates.DayLocator(interval=30))
axf.xaxis.set_major_formatter(DateFormatter("%m/%d"))
axf.legend()
axf.set_title('Sum (depth and transect) of the magnitude of the rotated flux for each transect')
axf.set_ylabel('SSC flux (kg/s)')


# plot cumulative sum of SSC for transect
fig, (axg) = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))

# since integrating need to start at second timestep for time.
axg.plot_date(datetime_list[1:], t0_cumtrapz_ssc, label='T0',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axg.text(datetime_list[-1],t0_cumtrapz_ssc[-1],'%.2e tons' % (t0_total_sed/1000))

axg.plot_date(datetime_list[1:], t1_cumtrapz_ssc, label='T1',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axg.text(datetime_list[-1],t1_cumtrapz_ssc[-1],'%.2e tons' % (t1_total_sed/1000))

axg.plot_date(datetime_list[1:], t2_cumtrapz_ssc, label='T2',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axg.text(datetime_list[-1],t2_cumsum_ssc[-1],'%.2e tons' % (t2_total_sed/1000))
axg.xaxis.set_major_locator(mdates.DayLocator(interval=30))
axg.xaxis.set_major_formatter(DateFormatter("%m/%d"))
axg.legend()
axg.set_title('Cumulative integral of the rotated flux for transects')
axg.set_ylabel('Cumulative integral of SSC flux (kg/s)')
