import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import numpy as np
import netCDF4
import coawstpy
import scipy.integrate as integrate
from scipy import stats
import math

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

lat = f.variables['lat_rho'][:]
lon = f.variables['lon_rho'][:]

Huon = f.variables['Huon'][0:-1,:,:,:]  # lon_u east-west (from panoply should be S)
Hvom = f.variables['Hvom'][0:-1,:,:,:]  # lon_v north-south (from panoply should be E)

mask_rho = f.variables['mask_rho'][:]

#s_rho = f.variables['s_rho'][:]  # depth levels

tx1 =0
tx2 = -1
#tx1 = 1202
#tx2 = 1466
#srho_angle = 3
datetime_list=[]
for sec in ocean_time:
    datetime_list.append(
        netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))


transect = dict()
##############################
# SR entrance Model boundary #
##############################
transect['Tb'] = dict()
transect['Tb']['x'] = np.array(list(range(20,29)))
transect['Tb']['y'] = np.array([0]*len(transect['Tb']['x']))

###########################
# Susquehanna River mouth #
###########################
transect['T1'] = dict()
transect['T1']['x'] = np.array([29,30,31,32])
transect['T1']['y'] = np.array([13,12,11,10])
#transect['T1']['x'] = np.array([28,29,30,31,32,33])
#transect['T1']['y'] = np.array([14,13,12,11,10,9])

###############################
# Turkey Point to Sandy Point #
###############################
transect['T2'] = dict()
transect['T2']['x'] = np.array(list(range(42,67)))
transect['T2']['y'] = np.array([58]*len(transect['T2']['x']))
#fig, (ax) = plt.subplots(nrows=3, ncols=1, figsize=(12, 12))
for t in transect: # each transect
    fig, (ax) = plt.subplots(nrows=3, ncols=1, figsize=(12, 12))
    x = transect[t]['x']
    y = transect[t]['y']
    for j in range(len(x)):
        mask_rho[x[j], y[j]] = 5
    transect[t]['Hvom'] = Hvom[:, :, x, y]
    transect[t]['Huon'] = Huon[:, :, x, y]
    transect[t]['fs'] = np.sum(Huon[:, :, x, y], axis=(1,2))
    transect[t]['fe'] = np.sum(Hvom[:, :, x, y], axis=(1,2))
    ax[0].plot(transect[t]['fs'],label='fs')
    ax[0].plot(transect[t]['fe'],label='fe')
    ax[0].set_ylabel('Raw Flux [m3/s]')
    ax[0].legend()
    # TODO plotting
    ax[1].plot(transect[t]['fs'], transect[t]['fe'], marker='.', linestyle='')

    ax[1].set_xlabel('South Flux (fs) [m3/s]', fontsize=8)
    ax[1].set_ylabel('East Flux (fe) [m3/s]', fontsize=8)
    xmin = np.min([transect[t]['fs'].min(), transect[t]['fe'].min()])
    xmax = np.max([transect[t]['fs'].max(), transect[t]['fe'].max()])
    ax[1].set_xlim(xmin, xmax)
    ax[1].set_ylim(xmin, xmax)
    ax[1].grid(axis='both')
    ax[1].set_aspect('equal', 'box')
    if t is 'Tb':
        ax[1].set_yticks([5000, 10000, 15000])
    slope, intercept, r_value, p_value, std_err = stats.linregress(transect[t]['fs'], transect[t]['fe'])
    xs = np.array([xmin, xmax])
    ax[1].plot(xs, slope*xs+intercept, linestyle='-', color='k')
#    theta = math.degrees(math.atan(slope))
    theta = math.atan(slope)
    ax[1].text(8000, 8000, 'theta = %.1f' % theta)
    # TODO work through trigonometry here:
    transect[t]['fd'] = transect[t]['fs']*np.cos(theta) + transect[t]['fe']*np.sin(theta)
    transect[t]['fa'] = -1*transect[t]['fs']*np.sin(theta) + transect[t]['fe']*np.cos(theta)
    ax[2].plot(transect[t]['fa'],label='fa = -fs*sin(theta) + fe*cos(theta)')
    ax[2].plot(transect[t]['fd'], label='fd = fs*cos(theta) + fe*sin(theta)')
    ax[2].set_ylabel('Rotated Flux [m/s]')
    ax[2].legend()
    fig.suptitle(t)

#fig.tight_layout()


plt.figure()
plt.pcolor(lon, lat, mask_rho, edgecolors='k', cmap='PuBu')
plt.title('Transects')