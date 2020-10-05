'''
Biddle et al 2020 manuscript Figure 3:
Caption:
Timeseries (top panel) of the magnitude of the water velocity (m/s) from CBIBS bin 2 ADCP observations (dots)
and the respective location within the model grid (black line). The bottom two panels are comparisons between
the observations (y-axis) and the model predictions (x-axis) (m/s), for the North (left panel) and East (right
panel) components of the water velocity (m/s). The red dots indicate the velocities between 6 September and 20
September and the black dots are the observations outside that time interval. The solid line is the linear
regression of the observed and predicted velocity components for the entire simulation (black and red dots
combined) and the dashed line is the 1:1 line.

@author: Mathew Biddle
'''
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
import numpy as np
import datetime
import netCDF4
import coawstpy
import pandas as pd
from scipy import stats

# bring in the data
# run = 'veg'
#
# if run == 'noveg':
#     dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final_noveg'
# else:
#     dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'


#inputfile = dir+'/upper_ches_his.nc'
inputfile = coawstpy.get_file_paths()['veg']
print("Reading history file %s..." % inputfile)
f = netCDF4.Dataset(inputfile, 'r')
ocean_time = f.variables['ocean_time'][:]
lat = f.variables['lat_rho'][:]
lon = f.variables['lon_rho'][:]
## Do some date conversions ##
datetime_list=[]
for sec in ocean_time:
    datetime_list.append(netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))

## Raw CBIBS data
CBIBS_file = coawstpy.get_file_paths()['cbibs']
print("Reading CBIBS file %s" % CBIBS_file)
#CBIBS_file = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/Initialization_data/CBIBS_insitu_obs/NCEI_copy/S_2011.nc'
fcbibs = netCDF4.Dataset(CBIBS_file, 'r')
cbibs_lon = -76.0736
cbibs_lat = 39.5404
cbibs_u = fcbibs.variables['eastward_water_velocity'][:, 2] # time,depth
cbibs_v = fcbibs.variables['northward_water_velocity'][:, 2] # time,depth
cbibs_flag = fcbibs.variables['current_flag'][:, 2]
#cbibs_flag[cbibs_flag == 9] = 4
flag_d = {1: 'black', 2: 'blue', 3: 'orange', 4: 'red', 5: 'yellow'}
#flag_values: [1 2 3 4 9]
#flag_meanings: good not_evaluated suspect bad missing
cbibs_mag = np.sqrt(cbibs_u**2 + cbibs_v**2)
cbibs_time = fcbibs.variables['time3'][:]
cbibs_date=[]
for days in cbibs_time:
    cbibs_date.append(datetime.datetime.fromordinal(int(days)) + datetime.timedelta(days=days%1) - datetime.timedelta(days = 366))

lat_pt = cbibs_lat
lon_pt = cbibs_lon
x = np.abs(f.variables['lon_rho'][:, 1]-lon_pt).argmin()
y = np.abs(f.variables['lat_rho'][1, :]-lat_pt).argmin()
ubar = f.variables['ubar_eastward'][:, x, y]
vbar = f.variables['vbar_northward'][:, x, y]
mag_vel = np.sqrt(ubar**2 + vbar**2)

# load subset CBIBS data
mid_water_current = coawstpy.get_file_paths()['mid_wtr_current_veg']
print("Reading mid-water currents %s" % mid_water_current)
#df = pd.read_csv(dir+'/mid_water_currents.csv', index_col=0, parse_dates=True)
df = pd.read_csv(mid_water_current, index_col=0, parse_dates=True)

start_date = '2011-08-01'
end_date = '2011-10-31'

start_date_1 = '2011-08-01'
end_date_1 = '2011-09-06'

start_date_2 = '2011-09-06'
end_date_2 = '2011-09-20'

start_date_3 = '2011-09-20'
end_date_3 = '2011-10-31'

print('Generating figues...')
#fig, ((ax),(ax2,ax3)) = plt.subplots(nrows=2, ncols=2, figsize=(10, 5))
plt.figure(figsize=(12, 8))
ax1 = plt.subplot(211)
ax2 = plt.subplot(223)
ax3 = plt.subplot(224)

ax1.plot_date(datetime_list, mag_vel, label='Predicted',
              xdate=True, linestyle='-', linewidth=2,
              marker='', markersize=1, color='k')

xlim = ax1.get_xlim()
ax1.scatter(cbibs_date, cbibs_mag, label='Observed', c='grey', s=3)
indx1 = coawstpy.nearest_ind(cbibs_date,datetime.datetime.strptime(start_date_2,'%Y-%m-%d'))
indx2 = coawstpy.nearest_ind(cbibs_date,datetime.datetime.strptime(end_date_2,'%Y-%m-%d'))
ax1.scatter(cbibs_date[indx1:indx2],cbibs_mag[indx1:indx2], c='red', s=3)
ax1.xaxis.set_major_locator(mdates.DayLocator(interval=10))
ax1.xaxis.set_major_formatter(DateFormatter("%m/%d"))

ax1.set_ylim([0, 2.5])
ax1.set_xlim(xlim)
ax1.set_ylabel('Water Velocity ($m$ $s^{-1}$)')
ax1.legend(loc='upper left')
#a = ax2.scatter(x=df.loc[start_date:end_date, 'COAWST_Vbar'], y=df.loc[start_date:end_date,'CBIBS_V'],
#                c='grey', s=5)#, colormap='viridis')
a = ax2.scatter(x=df.loc[start_date_1:end_date_1, 'COAWST_Vbar'],y=df.loc[start_date_1:end_date_1,'CBIBS_V'],marker='.',c='grey',s=15,label='other times')
a = ax2.scatter(x=df.loc[start_date_3:end_date_3, 'COAWST_Vbar'],y=df.loc[start_date_3:end_date_3,'CBIBS_V'],marker='.',c='grey',s=15)
a = ax2.scatter(x=df.loc[start_date_2:end_date_2, 'COAWST_Vbar'],y=df.loc[start_date_2:end_date_2,'CBIBS_V'],marker='.',c='red',s=15,label='09/06 - 09/20')

#cm = df.loc['2011-08-01':'2011-09-09', ['COAWST_Ubar','CBIBS_U']].plot.scatter(
#    x='COAWST_Ubar', y='CBIBS_U', c=df.loc['2011-08-01':'2011-09-09'].index, colormap='viridis', ax=ax)
#cbar = fig.colorbar(a, ax=ax[0])
#cbar.ax.set_yticklabels(pd.to_datetime(cbar.get_ticks()).strftime(date_format='%b %d'))
xmin = np.min([df.loc[start_date:end_date, 'COAWST_Vbar'].min(), df.loc[start_date:end_date,'CBIBS_V'].min()])
xmax = np.max([df.loc[start_date:end_date, 'COAWST_Vbar'].max(), df.loc[start_date:end_date,'CBIBS_V'].max()])
ax2.set_xlim(xmin, xmax)
ax2.set_ylim(xmin, xmax)
ax2.grid(axis='both')
ax2.set_aspect('equal', 'box')
ax2.set_ylabel('Observed')
ax2.set_xlabel('Predicted')
ax2.set_title('V - North')
#ax2.legend()
slope, intercept, r_value, p_value, std_err = stats.linregress(
    df.loc[start_date:end_date, 'COAWST_Vbar'].values, df.loc[start_date:end_date,'CBIBS_V'].values)
xs = np.array([xmin, xmax])
ax2.plot(xs, slope*xs+intercept, linestyle='-', color='k')
ax2.plot([-2,2],[-2,2],linestyle=':', color='k')

#a = ax3.scatter(x=df.loc[start_date:end_date, 'COAWST_Ubar'],y=df.loc[start_date:end_date,'CBIBS_U'],
#                c='grey', s=5)#, colormap='viridis')
a = ax3.scatter(x=df.loc[start_date_1:end_date_1, 'COAWST_Ubar'],y=df.loc[start_date_1:end_date_1,'CBIBS_U'],marker='.',c='grey',s=15,label='other times')
a = ax3.scatter(x=df.loc[start_date_3:end_date_3, 'COAWST_Ubar'],y=df.loc[start_date_3:end_date_3,'CBIBS_U'],marker='.',c='grey',s=15)
a = ax3.scatter(x=df.loc[start_date_2:end_date_2, 'COAWST_Ubar'],y=df.loc[start_date_2:end_date_2,'CBIBS_U'],marker='.',c='red',s=15,label='09/06 - 09/20')

xmin = np.min([df.loc[start_date:end_date, 'COAWST_Ubar'].min(), df.loc[start_date:end_date,'CBIBS_U'].min()])
xmax = np.max([df.loc[start_date:end_date, 'COAWST_Ubar'].max(), df.loc[start_date:end_date,'CBIBS_U'].max()])
ax3.set_xlim(xmin, xmax)
ax3.set_ylim(xmin, xmax)
ax3.grid(axis='both')
ax3.set_aspect('equal', 'box')
ax3.set_ylabel('Observed')
ax3.set_xlabel('Predicted')
ax3.set_title('U - East')
ax3.set_yticks(ticks=np.arange(-0.5,1.5,0.5))
#ax3.yticks(ticks=[])

slope, intercept, r_value, p_value, std_err = stats.linregress(
    df.loc[start_date:end_date, 'COAWST_Ubar'].values, df.loc[start_date:end_date,'CBIBS_U'].values)
xs = np.array([xmin, xmax])
ax3.plot(xs, slope*xs+intercept, linestyle='-', color='k')
ax3.plot([-2,2],[-2,2],linestyle=':',color='k')

ax2.legend(bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0.)
#plt.suptitle('%s to %s' % (start_date,end_date))
#ax.grid(True)
#fig.colorbar(format=DateFormatter('%b %d'))
#cb.ax.set_yticklabels(df.loc['2011-08-01':'2011-09-09'].index)
#cb = fig.colorbar(smap,format=DateFormatter('%d %b %y'))

# writedir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/Manuscript/figures/'
# image_name = 'Fig_3.png'
# outfile = writedir+image_name
# print("Saving image to %s" % outfile)
#plt.savefig(outfile, bbox_inches='tight', dpi=500)

print("Done")
