import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import numpy as np
import netCDF4
import coawstpy
import scipy.integrate as integrate

'''
Computes the water flux across transects by computing the flux sum in North and East direction, then translating to
along and cross channel fluxes.
'''

# bring in the data
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final_noveg'
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'
if dir.split("_")[-1] == 'noveg':
    run = "noveg"
else:
    run = "veg"
inputfile = dir+'/upper_ches_avg.nc'
print('Reading %s...' % inputfile.split("/")[-1])
f = netCDF4.Dataset(inputfile, 'r')

# Get the data we want
ocean_time = f.variables['ocean_time'][:]

lat_v = f.variables['lat_v'][:]
lon_v = f.variables['lon_v'][:]
mask_v = f.variables['mask_v'][:]

lat_u = f.variables['lat_u'][:]
lon_u = f.variables['lon_u'][:]
mask_u = f.variables['mask_u'][:]

lat_rho = f.variables['lat_rho'][:]
lon_rho = f.variables['lon_rho'][:]

lat = lat_rho
lon = lon_rho

Huon = f.variables['Huon'][0:-1,:,:,:]  # lon_u east-west (from panoply should be S)
Hvom = f.variables['Hvom'][0:-1,:,:,:]  # lon_v north-south (from panoply should be E)

plant_height = f.variables['mask_rho'][:]

s_rho = f.variables['s_rho'][:]  # depth levels

tx1 =0
tx2 = -1
#tx1 = 1202
#tx2 = 1466
srho_angle = 3
datetime_list=[]
for sec in ocean_time:
    datetime_list.append(
        netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))

##############################
# SR entrance Model boundary #
##############################
trans_name = 'Tb'
print('Extracting data for transect %s...' % trans_name)
xb = np.array(list(range(20,29)))
yb = np.array([0]*len(xb))
# Verify point location
for i in range(len(xb)):
    l = i/5
    plant_height[xb[i], yb[i]] = 5#l*100

# Gather subset data
#tb_flux_xe = np.sum(Hvom[:, :, xb, yb], axis=(1,2))
#tb_flux_yn = np.sum(Huon[:, :, xb, yb], axis=(1,2))

tb_flux_xe = np.sum(Hvom[:, :, xb, yb], axis=(1,2))
tb_flux_yn = np.sum(Huon[:, :, xb, yb], axis=(1,2))



Tb = coawstpy.sediment_flux3(tb_flux_xe, tb_flux_yn, tx1, tx2)
Tb['name'] = trans_name

###############
# SR entrance #
###############
# trans_name = 'T0'
# print('Extracting data for transect %s...' % trans_name)
# x0 = np.array([27,26,25,24,23])
# y0 = np.array([0,0,1,2,3])
# # Verify point location
# for i in range(len(x0)):
#     l = i/5
#     plant_height[x0[i], y0[i]] = 5
#
# # Gather subset data
# t0_flux_xe = np.sum(Hvom[:, :, x0, y0], axis=(1,2))
# t0_flux_yn = np.sum(Huon[:, :, x0, y0], axis=(1,2))
#
# T0 = coawstpy.sediment_flux3(t0_flux_xe, t0_flux_yn, tx1, tx2)
# T0['name'] = trans_name

###########################
# Susquehanna River mouth #
###########################
trans_name = 'T1'
print('Extracting data for transect %s...' % trans_name)
#x1 = np.array(list(range(27,33)))
#y1 = np.array([10]*len(x1))
x1 = np.array([28,29,30,31,32,33])
y1 = np.array([14,13,12,11,10,9])

# Verify point location
for i in range(len(x1)):
    plant_height[x1[i], y1[i]] = 5

# Gather subset data
print('hvom:', Hvom[0, 0, x1, y1])
print('huon:', Huon[0, 0, x1, y1])
#sys.exit()
t1_flux_xe = np.sum(Hvom[:, :, x1, y1], axis=(1,2))
t1_flux_yn = np.sum(Huon[:, :, x1, y1], axis=(1,2))

T1 = coawstpy.sediment_flux3(t1_flux_xe, t1_flux_yn, tx1, tx2)
T1['name'] = trans_name

###############################
# Turkey Point to Sandy Point #
###############################
trans_name = 'T2'
print('Extracting data for transect %s...' % trans_name)
x2 = np.array(list(range(42,67)))  #
y2 = np.array([58]*len(x2))

# Verify point location
for i in range(len(x2)):
    l = i/5
    plant_height[x2[i], y2[i]] = 5#l*100

# Gather subset data
t2_flux_xe = np.sum(Hvom[:, :, x2, y2], axis=(1,2))
t2_flux_yn = np.sum(Huon[:, :, x2, y2], axis=(1,2))

T2 = coawstpy.sediment_flux3(t2_flux_xe, t2_flux_yn, tx1, tx2)
T2['name'] = trans_name


## print out some information
print('Total Water Flux across transect over %s seconds' % (datetime_list[-1]-datetime_list[1]).total_seconds())
trans = [Tb,T1,T2] # group all transects together to loop over
for t in trans:
    mag = np.sqrt(t['xe'] ** 2 + t['yn'] ** 2)
    raw_flux = integrate.cumtrapz(mag, axis=0)
    print('Raw Flux: %s = %e m3/s' % (t['name'], raw_flux[-1]))
    print('Rotated Flux: %s = %e m3/s' % (t['name'], t['cumtrapz'][-1]))
    print('Along channel flux = %e m3/s' % integrate.cumtrapz(t['ux_rot'], axis=0)[-1])
    print('Cross channel flux = %e m3/s\n' % integrate.cumtrapz(t['vy_rot'], axis=0)[-1])

## Create some plots
# pick a depth and cell for investigation
# srho = depth coordinate
# c_* cell number along transect (t1 and t2).
srho = 3
c_t0 = 3
c_t1 = 1
c_t2 = 12

print("Making plots...")
# angle determination
## TODO
plt.figure(); plt.plot(tb_flux_xe,tb_flux_yn,marker='.',linestyle='');plt.title('T0')
xmin = np.min([tb_flux_xe.min(),tb_flux_yn.min()])
xmax = np.min([tb_flux_xe.max(),tb_flux_yn.max()])
plt.xlim(xmin, xmax)
plt.ylim(xmin,xmax)
plt.axis('square')
plt.figure(); plt.plot(t1_flux_xe,t1_flux_yn,marker='.',linestyle='');plt.title('T1')
xmin = np.min([t1_flux_xe.min(),t1_flux_yn.min()])
xmax = np.min([t1_flux_xe.max(),t1_flux_yn.max()])
plt.xlim(xmin, xmax)
plt.ylim(xmin,xmax)
plt.axis('square')
sys.exit()
# map of transect
plt.figure()
plt.pcolor(lon, lat, plant_height, edgecolors='k', cmap='PuBu')
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


# plot fluxes for transects
varlist = ['cumtrapz','mag']
trans = [Tb,T1,T2]
for var in varlist:
    fig, (axg) = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))
    for t in trans:
        if var == 'mag':
            axg.plot_date(datetime_list, t[var], label=t['name'],
                          xdate=True, linestyle='-', linewidth=0.5,
                          marker='', markersize=1)
        else:
            axg.plot_date(datetime_list[1:],t[var],label=t['name'],
                      xdate=True, linestyle='-', linewidth=0.5,
                      marker='', markersize=1)
            axg.text(datetime_list[-1],t[var][-1],'%.2e m3/s' %(t[var][-1]))

    axg.xaxis.set_major_locator(mdates.DayLocator(interval=30))
    axg.xaxis.set_major_formatter(DateFormatter("%m/%d"))
    axg.set_yscale('log')
    axg.legend()
    if var == 'cumtrapz':
        axg.set_title('Cumulative integral of the rotated flux for transects %s' % run)
        axg.set_ylabel('Cumulative integral of %s flux (m3/s)' % var.split("_")[-1])
    else:
        axg.set_title('Magnitude of the rotated flux for transects %s' % run)
        axg.set_ylabel('Magnitude of %s flux (m3/s)' % var.split("_")[-1])

sys.exit()


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
axg.plot_date(datetime_list[1:], tb_cumtrapz_ssc, label='Tb',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
axg.text(datetime_list[-1],tb_cumtrapz_ssc[-1],'%.2e tons' % (tb_total_sed/1000))

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