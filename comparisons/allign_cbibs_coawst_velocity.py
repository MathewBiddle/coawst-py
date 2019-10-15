import matplotlib.pyplot as plt
import numpy as np
import datetime
import netCDF4
import coawstpy
import pandas as pd
from scipy import stats

# bring in the data
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final_noveg'
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'
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

print("Writing to "+dir+"/mid_water_currents.csv")
df.to_csv(dir+'/mid_water_currents.csv')
