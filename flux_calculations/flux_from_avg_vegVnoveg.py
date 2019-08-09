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
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final_noveg'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'
if dir.split("_")[-1] == 'noveg':
    run = "noveg"
else:
    run = "veg"
inputfile = dir+'/upper_ches_avg.nc'
print('Reading %s...' % inputfile.split("/")[-1])
f = netCDF4.Dataset(inputfile, 'r')

# Get the data we want
ocean_time = f.variables['ocean_time'][:]
#lat = f.variables['lat_v'][:]
#lon = f.variables['lon_v'][:]

# from https://github.com/trondkr/romstools/blob/master/VolumeFlux/tools.py
## TODO look at some of the scripts at https://github.com/ESMG/pyroms for using lat_v coordinates
## also some scripts at https://github.com/bjornaa/roppy/tree/master/roppy
lat_v = f.variables['lat_v'][:]
lon_v = f.variables['lon_v'][:]
lat_u = f.variables['lat_u'][:]
lon_u = f.variables['lon_u'][:]
mask_v = f.variables['mask_v'][:]
sys.exit()
# adjust to uniform grid
# https://github.com/ESMG/pyroms/blob/df6f90698e6b5903a1d484a738d693595f5cf213/pyroms/pyroms/utility.py#L400
# u ==> psi
# 0.5 * (varin[:,1:,:] + varin[:,:-1,:])
# v ==> psi
# 0.5 * (varin[:,:,1:] + varin[:,:,:-1])
# TODO figure out if that transformation is correct. nervous about non-contiguous grids at SR mouth

lat = 0.5 * (lat_v[:,1:] + lat_v[:,:-1]) # v ==> psi
lon = 0.5 * (lon_v[:,1:] + lon_v[:,:-1]) # v ==> psi

#lat = 0.5 * (lat_u[1:,:] + lat_u[:-1,:])
#lon = 0.5 * (lon_u[1:,:] + lon_u[:-1,:])


Huon_sand_01 = f.variables['Huon_sand_01'][:]  # lon_u east-west (from panoply should be S)
Hvom_sand_01 = f.variables['Hvom_sand_01'][:]  # lon_v north-south (from panoply should be E)
Huon_mud_01 = f.variables['Huon_mud_01'][:]  # lon_u east-west
Hvom_mud_01 = f.variables['Hvom_mud_01'][:]  # lon_v north-south

Huon_sand_01 = 0.5 * (Huon_sand_01[:,:,1:,:] + Huon_sand_01[:,:,:-1,:]) # u ==> psi
Hvom_sand_01 = 0.5 * (Hvom_sand_01[:,:,:,1:] + Hvom_sand_01[:,:,:,:-1]) # v ==> psi

Huon_mud_01 = 0.5 * (Huon_mud_01[:,:,1:,:] + Huon_mud_01[:,:,:-1,:]) # u ==> psi
Hvom_mud_01 = 0.5 * (Hvom_mud_01[:,:,:,1:] + Hvom_mud_01[:,:,:,:-1]) # v ==> psi

plant_height = Hvom_sand_01[0,0,:,:]
#plant_height = f.variables['Hvom_sand_01'][0, 0, :, :]
#plant_height = plant_height[:, 1:]

s_rho = f.variables['s_rho'][:]  # depth levels

tx1 = 0#1202
tx2 = -1#1466
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
x = np.array(list(range(20,28)))
y = np.array([0]*len(x))
# Verify point location
for i in range(len(x)):
    l = i/5
    plant_height[x[i], y[i]] = 5#l*100

# Gather subset data
tb_mud_01_flux_xe = Hvom_mud_01[:, :, x, y]
tb_mud_01_flux_yn = Huon_mud_01[:, :, x, y]
tb_sand_01_flux_xe = Hvom_sand_01[:, :, x, y]
tb_sand_01_flux_yn = Huon_sand_01[:, :, x, y]

Tb = coawstpy.sediment_flux(tb_mud_01_flux_xe,tb_mud_01_flux_yn,tb_sand_01_flux_xe,tb_sand_01_flux_yn,srho_angle,tx1,tx2)
Tb['name'] = trans_name
Tb['c_t'] = 5

###############
# SR entrance #
###############
trans_name = 'T0'
print('Extracting data for transect %s...' % trans_name)
x = np.array([27,26,25,24,23])
#x = np.array(list(range(20,28)))
#y = np.array([0]*len(x))
y = np.array([0,0,1,2,3])
# Verify point location
for i in range(len(x)):
    l = i/5
    plant_height[x[i], y[i]] = 10

# Gather subset data
t0_mud_01_flux_xe = Hvom_mud_01[:, :, x, y]
t0_mud_01_flux_yn = Huon_mud_01[:, :, x, y]
t0_sand_01_flux_xe = Hvom_sand_01[:, :, x, y]
t0_sand_01_flux_yn = Huon_sand_01[:, :, x, y]

T0 = coawstpy.sediment_flux(t0_mud_01_flux_xe,t0_mud_01_flux_yn,t0_sand_01_flux_xe,t0_sand_01_flux_yn,srho_angle,tx1,tx2)
T0['name'] = trans_name
T0['c_t'] = 3

###########################
# Susquehanna River mouth #
###########################
trans_name = 'T1'
print('Extracting data for transect %s...' % trans_name)
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

T1 = coawstpy.sediment_flux(t1_mud_01_flux_xe,t1_mud_01_flux_yn,t1_sand_01_flux_xe,t1_sand_01_flux_yn,srho_angle,tx1,tx2)
T1['name'] = trans_name
T1['c_t'] = 1

###############################
# Turkey Point to Sandy Point #
###############################
trans_name = 'T2'
print('Extracting data for transect %s...' % trans_name)
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

T2 = coawstpy.sediment_flux(t2_mud_01_flux_xe,t2_mud_01_flux_yn,t2_sand_01_flux_xe,t2_sand_01_flux_yn,srho_angle,tx1,tx2)
T2['name'] = trans_name
T2['c_t'] = 12
## print out some information

print('Total Sediment across transect over %s seconds' % (datetime_list[-1]-datetime_list[1]).total_seconds())
trans = [Tb,T0,T1,T2] # group all transects together to loop over
for t in trans:
    print('%s = %e kg = %e tons' % (t['name'],t['total_sed'], t['total_sed']/1000))
#sys.exit()

## Create some plots
# pick a depth and cell for investigation
# srho = depth coordinate
# c_* cell number along transect (t1 and t2).
srho = 3
c_t0 = 3
c_t1 = 1
c_t2 = 12

print("Making plots...")
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

# plot cumulative sum of SSC for transect
varlist = ['cumtrapz_ssc','cumtrapz_mud','cumtrapz_sand','mag_mud', 'mag_sand', 'mag_ssc']
trans = [Tb,T1,T2]
for var in varlist:
    fig, (axg) = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))
    for t in trans:
        if var.split("_")[0] == 'mag':
            axg.plot_date(datetime_list, np.sum(t[var],axis=(1,2)), label=t['name'],
                          xdate=True, linestyle='-', linewidth=0.5,
                          marker='', markersize=1)
            #axg.text(datetime_list[-1], np.sum(t[var],axis=(1,2))[-1], '%.2e tons' % ((np.sum(t[var],axis=(1,2))[-1] * 3600) / 1000))
        else:
            axg.plot_date(datetime_list[1:],t[var],label=t['name'],
                      xdate=True, linestyle='-', linewidth=0.5,
                      marker='', markersize=1)
            axg.text(datetime_list[-1],t[var][-1],'%.2e tons' %((t[var][-1] * 3600)/1000))

    axg.xaxis.set_major_locator(mdates.DayLocator(interval=30))
    axg.xaxis.set_major_formatter(DateFormatter("%m/%d"))
    axg.set_yscale('log')
    axg.legend()
    if var.split("_")[0] == 'cumtrapz':
        axg.set_title('Cumulative integral of the rotated flux for transects %s' % run)
        axg.set_ylabel('Cumulative integral of %s flux (kg/s)' % var.split("_")[-1])
    else:
        axg.set_title('Magnitude of the rotated flux for transects %s' % run)
        axg.set_ylabel('Magnitude of %s flux (kg/s)' % var.split("_")[-1])

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

# plot rotated flux values at a specific depth and cell
i=0
fig, (axc) = plt.subplots(nrows=3, ncols=1, figsize=(12, 8),sharex=True)
fig.subplots_adjust(hspace=0.25)
for t in trans:
    axc[i].plot_date(datetime_list, t['mud_flux_ux_rot'][:, srho, t['c_t']], label='ux',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
    axc[i].plot_date(datetime_list,t['mud_flux_vy_rot'][:, srho, t['c_t']],label='vy',
              xdate=True, linestyle='-', linewidth=0.5,
              marker='', markersize=1)
    axc[i].set_title('Along (ux) and Cross (vy) channel mud_01 flux at s-rho=%i, cell=%i, Transect %s, angle = %i' %
                 (srho, c_t0, t['name'], t['angle']))
    axc[i].set_ylabel('mud_01 flux (kg/s)')
    i+=1
axc[2].xaxis.set_major_locator(mdates.DayLocator(interval=30))
axc[2].xaxis.set_major_formatter(DateFormatter("%m/%d"))
axc[0].legend()

sys.exit()
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