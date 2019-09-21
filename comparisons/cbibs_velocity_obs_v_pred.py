import matplotlib.pyplot as plt
import numpy as np
import datetime
import netCDF4
import coawstpy
import pandas as pd
from scipy import stats

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

# read CBIBS data
print('Reading CBIBS data...')
CBIBS_file = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/Initialization_data/CBIBS_insitu_obs/NCEI_copy/S_2011.nc'
fcbibs = netCDF4.Dataset(CBIBS_file, 'r')
# subset indexes for time. Probably don't need this after we do the searching for matching data.
#itime_start = 1059
#itime_end = 3618

cbibs_u = fcbibs.variables['eastward_water_velocity'][:, 2] # time,depth
cbibs_v = fcbibs.variables['northward_water_velocity'][:, 2] # time,depth
cbibs_mag = np.sqrt(cbibs_u**2 + cbibs_v**2)
cbibs_time = fcbibs.variables['time3'][:]

#Geolocation is 39.5404, -76.0736
cbibs_lon = -76.0736 # using coords from website
cbibs_lat = 39.5404

# To get all coordinates (which varies)
#cbibs_lat = fcbibs.variables['latitude'][:]
#cbibs_lon = fcbibs.variables['longitude'][:]

cbibs_date=[]
for days in cbibs_time:
    #cbibs_date.append(
    #    netCDF4.num2date(days, units=fcbibs.variables['time3'].units))
    cbibs_date.append(datetime.datetime.fromordinal(int(days)) + datetime.timedelta(days=days%1) - datetime.timedelta(days = 366))

lat_pt = cbibs_lat
lon_pt = cbibs_lon

# read ROMS history file
print('Reading ROMS data...')
inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')

ocean_time = f.variables['ocean_time'][:]
lat = f.variables['lat_rho'][:]
lon = f.variables['lon_rho'][:]

## Do some date conversions ##
datetime_list=[]
for sec in ocean_time:
    datetime_list.append(netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))
# find closest point in grid
try:
    lat_pt
    lon_pt
except NameError:
    print("Using indexes x, y = (%i, %i)" % (x, y))
else:
    print("Using geo-coords lat, lon = (%f, %f)" % (lat_pt, lon_pt))
    x = coawstpy.nearest_ind(f.variables['lon_rho'][:, 1],lon_pt)
    y = coawstpy.nearest_ind(f.variables['lat_rho'][1, :],lat_pt)

ubar = f.variables['ubar_eastward'][:, x, y]
vbar = f.variables['vbar_northward'][:, x, y]
mag_vel = np.sqrt(ubar**2 + vbar**2)


print('Matching nearest data...')
# TODO index = coawstpy.nearest_ind(datelist,value)
df = pd.DataFrame(columns=['CBIBS_time','CBIBS_U','CBIBS_V','COAWST_time','COAWST_Ubar','COAWST_Vbar'],
                  index=pd.date_range(start='2011-07-19T23:00:00', end='2011-11-01T00:00:00',freq='1H'))
## TODO build out data frame here!
i=0
idx=[]
for time in df.index:
    print("searching for time: %s" % time)
    cbibs_idx = coawstpy.nearest_ind(cbibs_date,time)
    print("Found cbibs time: %s" % cbibs_date[cbibs_idx])
    coawst_idx = coawstpy.nearest_ind(datetime_list,time)
    if cbibs_idx in idx:
        continue
    idx.append(cbibs_idx)
    print("Found ROMS time: %s\n" % datetime_list[coawst_idx])
    df.loc[time] = [cbibs_date[cbibs_idx],cbibs_u[cbibs_idx],cbibs_v[cbibs_idx],
                  datetime_list[coawst_idx],ubar[coawst_idx],vbar[coawst_idx]]
    i+=1

df.dropna(how='any', inplace=True)
df['CBIBS_U'] = pd.to_numeric(df['CBIBS_U'])
df['CBIBS_V'] = pd.to_numeric(df['CBIBS_V'])
df['COAWST_Vbar'] = pd.to_numeric(df['COAWST_Vbar'])
df['COAWST_Ubar'] = pd.to_numeric(df['COAWST_Ubar'])

start_date = '2011-08-01'
end_date = '2011-11-01'

fig, (ax) = plt.subplots(nrows=2, ncols=1, figsize=(6, 12))

a = ax[0].scatter(x=df.loc[start_date:end_date, 'COAWST_Vbar'],y=df.loc[start_date:end_date,'CBIBS_V'],
                c=df.loc[start_date:end_date].index)#, colormap='viridis')
#cm = df.loc['2011-08-01':'2011-09-09', ['COAWST_Ubar','CBIBS_U']].plot.scatter(
#    x='COAWST_Ubar', y='CBIBS_U', c=df.loc['2011-08-01':'2011-09-09'].index, colormap='viridis', ax=ax)
cbar = fig.colorbar(a, ax=ax[0])
cbar.ax.set_yticklabels(pd.to_datetime(cbar.get_ticks()).strftime(date_format='%b %d'))
xmin = np.min([df.loc[start_date:end_date, 'COAWST_Vbar'].min(), df.loc[start_date:end_date,'CBIBS_V'].min()])
xmax = np.max([df.loc[start_date:end_date, 'COAWST_Vbar'].max(), df.loc[start_date:end_date,'CBIBS_V'].max()])
ax[0].set_xlim(xmin, xmax)
ax[0].set_ylim(xmin, xmax)
ax[0].grid(axis='both')
ax[0].set_aspect('equal', 'box')
ax[0].set_ylabel('CBIBS')
ax[0].set_xlabel('COAWST')
ax[0].set_title('V - North')
slope, intercept, r_value, p_value, std_err = stats.linregress(
    df.loc[start_date:end_date, 'COAWST_Vbar'].values, df.loc[start_date:end_date,'CBIBS_V'].values)
xs = np.array([xmin, xmax])
ax[0].plot(xs, slope*xs+intercept, linestyle='-', color='k')
ax[0].plot([-1,1],[-1,1],linestyle=':')

a = ax[1].scatter(x=df.loc[start_date:end_date, 'COAWST_Ubar'],y=df.loc[start_date:end_date,'CBIBS_U'],
                c=df.loc[start_date:end_date].index)#, colormap='viridis')
#cm = df.loc['2011-08-01':'2011-09-09', ['COAWST_Ubar','CBIBS_U']].plot.scatter(
#    x='COAWST_Ubar', y='CBIBS_U', c=df.loc['2011-08-01':'2011-09-09'].index, colormap='viridis', ax=ax)
cbar = fig.colorbar(a, ax=ax[1])
cbar.ax.set_yticklabels(pd.to_datetime(cbar.get_ticks()).strftime(date_format='%b %d'))
xmin = np.min([df.loc[start_date:end_date, 'COAWST_Ubar'].min(), df.loc[start_date:end_date,'CBIBS_U'].min()])
xmax = np.max([df.loc[start_date:end_date, 'COAWST_Ubar'].max(), df.loc[start_date:end_date,'CBIBS_U'].max()])
ax[1].set_xlim(xmin, xmax)
ax[1].set_ylim(xmin, xmax)
ax[1].grid(axis='both')
ax[1].set_aspect('equal', 'box')
ax[1].set_ylabel('CBIBS')
ax[1].set_xlabel('COAWST')
ax[1].set_title('U - East')
slope, intercept, r_value, p_value, std_err = stats.linregress(
    df.loc[start_date:end_date, 'COAWST_Ubar'].values, df.loc[start_date:end_date,'CBIBS_U'].values)
xs = np.array([xmin, xmax])
ax[1].plot(xs, slope*xs+intercept, linestyle='-', color='k')
ax[1].plot([-1,1],[-1,1],linestyle=':')
plt.suptitle('%s to %s' % (start_date,end_date))
#ax.grid(True)
#fig.colorbar(format=DateFormatter('%b %d'))
#cb.ax.set_yticklabels(df.loc['2011-08-01':'2011-09-09'].index)
#cb = fig.colorbar(smap,format=DateFormatter('%d %b %y'))

sys.exit()
## Make some plots.
fig, (ax) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(12, 8))

ax[0].plot_date(datetime_list,mag_vel,label='COAWST',
              xdate=True, linestyle='', linewidth=1,
              marker='.', markersize=1, color='r')
#ax.spines['left'].set_color('red')
#ax.tick_params(axis='y', colors='red')
#ax.set_ylabel('COAWST vbar [m/s]')
#ax.yaxis.label.set_color('red')

xlim = ax[0].get_xlim()
#ax.grid(True)
#ax2v = ax.twinx()
ax[0].plot_date(cbibs_date,cbibs_mag, label='CBIBS',
              xdate=True, linestyle='', linewidth=0.5,
              marker='.', markersize=1, color='b')
#ax2v.set_ylabel('CBIBS v [m/s]')
#ax2v.tick_params(axis='y', colors='blue')
#ax2v.yaxis.label.set_color('blue')
#ax2v.grid(True)
#ax2v.set_ylim([-2, 0.25])
#ax2v.set_xlim(xlim)
ax[0].xaxis.set_major_locator(mdates.DayLocator(interval=10))
ax[0].xaxis.set_major_formatter(DateFormatter("%m/%d"))

ax[0].set_ylim([0, 2.5])
ax[0].set_xlim(xlim)
ax[0].set_ylabel('Water Velocity [m/s]')
ax[0].legend(loc='upper left')
plt.suptitle('%s'%run)

ax[1].plot(mag_vel,cbibs_mag,marker='.',linestyle='')
ax[1].set_xlabel('Predicted mag [m/s]')
ax[1].set_ylabel('CBIBS mag [m/s]')


#outfile = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/figures/CBIBS_Velocity_comparison.png'

#plt.savefig(outfile, bbox_inches='tight', dpi = 1000)
