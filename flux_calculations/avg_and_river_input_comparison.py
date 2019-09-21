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
import numpy.ma as ma

'''
Computes the water flux across transects by computing the flux sum in North and East direction, then translating to
along and cross channel fluxes.

Also computes the input flux from the ROMS river file

The figures produced from this script show that the raw South flux values we get from ROMS are the flux values we 
should be using during our sediment flux calculations
'''

# bring in the data
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final_noveg'
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final_post_lee'
if dir.split("_")[-1] == 'noveg':
    run = "noveg"
elif dir.split("_")[-1] == 'lee':
    run = 'post-lee'
else:
    run = "veg"
inputfile = dir+'/upper_ches_avg.nc'
print('Reading %s %s...' % (inputfile.split("/")[-1],run))
f = netCDF4.Dataset(inputfile, 'r')

tstart = 0
tend = -1 # 1000

# Get the data we want
ocean_time = f.variables['ocean_time'][tstart:tend]
lat = f.variables['lat_rho'][:]
lon = f.variables['lon_rho'][:]
mask_u = f.variables['mask_u'][:]
mask_v = f.variables['mask_v'][:]
Huon = f.variables['Huon'][tstart:tend,:,:,:]  # lon_u east-west (from panoply should be S)
Hvom = f.variables['Hvom'][tstart:tend,:,:,:]  # lon_v north-south (from panoply should be E)

mask_rho = f.variables['mask_rho'][:]

datetime_list=[]
for sec in ocean_time:
    datetime_list.append(
        netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))

## river data
river_frc = dir+'/river_frc.nc'
f_river = netCDF4.Dataset(river_frc, 'r')
river_time = f_river.variables['river_time'][:]
river_transport = f_river.variables['river_transport'][:] # middle of water column
river_datetime_list=[]
for sec in river_time:
    river_datetime_list.append(
        netCDF4.num2date(sec, units=f_river.variables['river_time'].units, calendar='standard'))

initial_riv_idx = coawstpy.nearest_ind(river_datetime_list,datetime_list[0])
final_riv_idx = coawstpy.nearest_ind(river_datetime_list,datetime_list[-1])+1 # have to add one for slicing to include last number
river_datetime_list_subset = river_datetime_list[initial_riv_idx : final_riv_idx : 120]
river_transport_subset = river_transport[initial_riv_idx : final_riv_idx : 120,:]

## Now build transects
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
#transect['T1']['x'] = np.array([27,28,29,30,31,32,33,34])
#transect['T1']['y'] = np.array([15,14,13,12,11,10,9,8])
transect['T1']['x'] = np.array(list(range(27,34)))
transect['T1']['y'] = np.array([10]*len(transect['T1']['x']))

###############################
# Turkey Point to Sandy Point #
###############################
transect['T2'] = dict()
transect['T2']['x'] = np.array(list(range(40,69)))
transect['T2']['y'] = np.array([58]*len(transect['T2']['x']))

###############################
#  North-East Boundary Check  #
###############################
#transect['T3'] = dict()
#transect['T3']['x'] = np.array(list(range(87,97)))
#transect['T3']['y'] = np.array([1]*len(transect['T3']['x']))

fig1, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
fig2, (ax2) = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
## Iterate through each transect
for t in transect: # each transect
    fig, (ax) = plt.subplots(nrows=3, ncols=1, figsize=(12, 12))

    # get index of interest
    x = transect[t]['x']
    y = transect[t]['y']
    transect[t]['Hvom'] = np.ma.empty((Hvom.shape[0],Hvom.shape[1],len(x)))
    transect[t]['Huon'] = np.ma.empty((Hvom.shape[0],Hvom.shape[1],len(x)))
    for j in range(len(x)):
        mask_rho[x[j], y[j]] = 5
        transect[t]['Hvom'][:,:,j] = Hvom[:,:,x[j],y[j]]
        transect[t]['Huon'][:,:,j] = Huon[:,:,x[j],y[j]]

    # collect data at index of interest
    #transect[t]['Hvom'] = Hvom[:, :, x, y]
    #transect[t]['Huon'] = Huon[:, :, x, y]

    # sum across transect and depth
    transect[t]['fs'] = np.sum(transect[t]['Huon'], axis=(1, 2))
    transect[t]['fe'] = np.sum(transect[t]['Hvom'], axis=(1, 2))

    # make a plot of raw south and raw east on top subplot
    ax[0].plot_date(datetime_list,transect[t]['fs'],label='fs',
                    linestyle='-', linewidth=1, marker='')
    ax[0].plot_date(datetime_list,transect[t]['fe'],label='fe',
                    linestyle='-', linewidth=1, marker='')
    ax[0].xaxis.set_major_locator(mdates.DayLocator(interval=10))
    ax[0].xaxis.set_major_formatter(DateFormatter("%m/%d"))
    ax[0].set_ylabel('Raw Flux [m3/s]')
    ax[0].legend()

    # make plot of raw south vs raw east
    ax[1].plot(transect[t]['fs'], transect[t]['fe'], marker='.', linestyle='')
    ax[1].set_xlabel('South Flux (fs) [m3/s]', fontsize=8)
    ax[1].set_ylabel('East Flux (fe) [m3/s]', fontsize=8)
    xmin = np.min([transect[t]['fs'].min(), transect[t]['fe'].min()])
    xmax = np.max([transect[t]['fs'].max(), transect[t]['fe'].max()])
    ax[1].set_xlim(xmin, xmax)
    ax[1].set_ylim(xmin, xmax)
    ax[1].grid(axis='both')
    ax[1].set_aspect('equal', 'box')
    #if t is 'Tb':
    #    ax[1].set_yticks([5000, 10000, 15000])

    # compute linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(transect[t]['fs'], transect[t]['fe'])
    xs = np.array([xmin, xmax])
    ax[1].plot(xs, slope*xs+intercept, linestyle='-', color='k')

    # compute angle linear regression is ccw from x axis
    theta = math.atan(slope)
    ax[1].set_title('theta = %.2f' % theta)

    # rotate to theta in radians
    transect[t]['fd'] = transect[t]['fs']*np.cos(theta) + transect[t]['fe']*np.sin(theta)
    transect[t]['fa'] = -1*transect[t]['fs']*np.sin(theta) + transect[t]['fe']*np.cos(theta)
    ax[2].plot_date(datetime_list,transect[t]['fa'], label='fa = -fs*sin(theta) + fe*cos(theta)',
                    linestyle='-', linewidth=1, marker='')
    ax[2].plot_date(datetime_list,transect[t]['fd'], label='fd = fs*cos(theta) + fe*sin(theta)',
                    linestyle='-', linewidth=1, marker='')
    ax[2].set_ylabel('Rotated Flux [m3/s]')
    ax[2].legend()
    ax[2].xaxis.set_major_locator(mdates.DayLocator(interval=10))
    ax[2].xaxis.set_major_formatter(DateFormatter("%m/%d"))
    fig.suptitle(t)

    print('%s max south flux = %g m3/s' % (t, transect[t]['fs'].max()))
    print('%s integrated south flux = %g m3/s\n' % (t, integrate.trapz(transect[t]['fs'])))

    print('%s max downstream flux = %g m3/s' % (t, transect[t]['fd'].max()))
    print('%s integrated downstream flux = %g m3/s\n\n' % (t, integrate.trapz(transect[t]['fd'])))

    ax1.plot_date(datetime_list, transect[t]['fd'], label=t, linestyle='-', linewidth=1, marker='')
    ax1.legend()
    ax1.xaxis.set_major_locator(mdates.DayLocator(interval=10))
    ax1.xaxis.set_major_formatter(DateFormatter("%m/%d"))
    ax1.set_ylabel('Total downstream flux m3/s')

    ax2.plot_date(datetime_list, transect[t]['fs'], label=t, linestyle='-', linewidth=1, marker='')
    ax2.legend()
    ax2.xaxis.set_major_locator(mdates.DayLocator(interval=10))
    ax2.xaxis.set_major_formatter(DateFormatter("%m/%d"))
    ax2.set_ylabel('Total south flux m3/s')
    #ax1.set_yscale('log')
#fig.tight_layout()

ax2.plot_date(river_datetime_list_subset,np.sum(river_transport_subset,1),label='river transport m3/s',linestyle='-', linewidth=1, marker='')
ax2.legend()

ax1.plot_date(river_datetime_list_subset,np.sum(river_transport_subset,1),label='river transport m3/s',linestyle='-', linewidth=1, marker='')
ax1.legend()

print('River max flux (discharge) = %g m3/s' % (np.sum(river_transport_subset,1).max()))
print('River integrated flux (discharge) = %g m3/s\n' % (integrate.trapz(np.sum(river_transport_subset,1))))

plt.figure()
plt.pcolor(lon, lat, mask_rho, edgecolors='k', cmap='PuBu')
plt.title('Transects')
