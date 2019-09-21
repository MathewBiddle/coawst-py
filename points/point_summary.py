import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import numpy as np
import netCDF4
import pandas as pd
import coawstpy

# Cindy's locations
locs = pd.DataFrame(columns=['Site', 'lat', 'lon','comment'])

locations = [
    ['1', 39.527, -76.061, 'Russ and Palinkas 2018'],
    ['2', 39.533, -76.061, ''],
    ['3', 39.515, -76.051, 'Middle of bed'],
    ['4', 39.505, -76.039, ''],
    ['5', 39.497, -76.036, ''],
    ['Lee7', 39.414, -76.079, ''],
    ['Lee6', 39.38, -76.088, ''],
    ['Lee5', 39.346, -76.197, ''],
    ['Lee2.5', 39.197, -76.311, ''],
    ['Lee2', 39.135, -76.328, ''],
    ['Lee0', 39.061, -76.328, ''],
    ['LeeS2', 38.757, -76.473, ''],
    ['CBIBS', 39.5396, -76.0741, 'CBIBS Susquehanna Flats'],
    ['Tripod', 39.4931, -76.0341, 'Larry tripod site']
    ]

i = 0
for location in locations:
    locs.loc[i] = location
    i += 1

z_pt = -1 # 0=bottom 4=surface (-1)

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
inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')
ocean_time = f.variables['ocean_time'][:]
lat = f.variables['lat_rho'][:]
lon = f.variables['lon_rho'][:]
## Do some date conversions ##
datetime_list=[]
for sec in ocean_time:
    if sec == 0.0:
        datetime_list.append(
            netCDF4.num2date(sec + 0.0000000000000001, units=f.variables['ocean_time'].units,
                             calendar=f.variables['ocean_time'].calendar))
    else:
        datetime_list.append(
            netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))


#river_frc = dir+'/river_frc.nc'
#f_river = netCDF4.Dataset(river_frc, 'r')
#river_time = f_river.variables['river_time'][:]
#river_transport = (f_river.variables['river_transport'][:, 0] + (0.2 * f_river.variables['river_transport'][:, 0])) * \
#                  f_river.variables['river_transport'].shape[1]
#river_datetime_list=[]
#for sec in river_time:
#    ts = sec/(12*3600)
#    river_datetime_list.append(
#        netCDF4.num2date(sec, units=f_river.variables['river_time'].units, calendar='standard'))

#initial_riv_idx = coawstpy.nearest_ind(river_datetime_list,datetime_list[0])
#final_riv_idx = coawstpy.nearest_ind(river_datetime_list,datetime_list[-1])+1 # have to add one for slicing to include last number
#river_datetime_list_subset = river_datetime_list[initial_riv_idx : final_riv_idx : 120]
#river_transport_subset = river_transport[initial_riv_idx : final_riv_idx : 120]


ptsfile = dir+"/tripod_wave.pts"
ptsdf = pd.read_fwf(ptsfile, header=4)
ptsdf.drop([0,1],axis=0,inplace=True)
ptsdf.rename(columns={'%       Time':'Time'},inplace=True)
ptsdf['Yp'] = ptsdf['Yp            Hsig'].astype(str).str.split("    ",expand=True)[0].astype(float)
ptsdf['Hsig'] = ptsdf['Yp            Hsig'].astype(str).str.split("    ",expand=True)[1].astype(float)
ptsdf.drop(columns=['Yp            Hsig'],inplace=True)
ptsdf['Time']=pd.to_datetime(ptsdf['Time'],format='%Y%m%d.%H%M%S',utc=True)
ptsdf['Hsig']=ptsdf['Hsig'].astype(float)
ptsdf['X-Windv']=ptsdf['X-Windv'].astype(float)
ptsdf['Y-Windv']=ptsdf['Y-Windv'].astype(float)

point_data = dict()
times = coawstpy.time_periods()
for event in times:
    start = coawstpy.nearest_ind(datetime_list,times[event][0])
    end = coawstpy.nearest_ind(datetime_list,times[event][1])
    point_data[event] = dict()
    # collect data for site of choice
    sites = ['CBIBS','3','4','Tripod','Lee7']
    for site in sites:
        point_data[event][site] = pd.DataFrame(columns=['swan_time','X-Windv','Y-Windv',
                                             'ocean_time','Pwave_Top','Hwave','mud','sand','bed_thickness','ubar','vbar'])
        #                                              'river_time','river_transport',

        lat_pt, lon_pt = locs.loc[locs['Site'] == site, ['lat', 'lon']].values[0]
        try:
            lat_pt
            lon_pt
        except NameError:
            print("Using indexes x, y = (%i, %i)" % (x, y))
        else:
            print("Using geo-coords lat, lon = (%f, %f)" % (lat_pt, lon_pt))
            x = np.abs(lon[:, 1]-lon_pt).argmin()
            y = np.abs(lat[1, :]-lat_pt).argmin()

        point_data[event][site]['swan_time'] = ptsdf['Time'][start:end]
        point_data[event][site]['X-Windv'] = ptsdf['X-Windv'][start:end]
        point_data[event][site]['Y-Windv'] = ptsdf['Y-Windv'][start:end]
        #point_data[event][site]['river_time'] = river_datetime_list_subset[start:end]
        #point_data[event][site]['river_transport'] = river_transport_subset[start:end]
        point_data[event][site]['ocean_time'] = datetime_list[start:end]
        point_data[event][site]['Pwave_Top'] = f.variables['Pwave_top'][start:end, x, y]
        point_data[event][site]['Hwave'] = f.variables['Hwave'][start:end,x,y]
        point_data[event][site]['mud'] = f.variables['mud_01'][start:end, z_pt, x, y]
        point_data[event][site]['sand'] = f.variables['sand_01'][start:end, z_pt, x, y]
        point_data[event][site]['bed_thickness'] = np.sum(f.variables['bed_thickness'][start:end,:,x,y], axis=1)
        point_data[event][site]['ubar_eastward'] = f.variables['ubar_eastward'][start:end, x, y]
        point_data[event][site]['vbar_northward'] = f.variables['vbar_northward'][start:end, x, y]
        # Verify point location
        #plant_height = f.variables['plant_height'][0, 0, :, :]
        #plant_height = np.ma.masked_greater(plant_height, 1)
        #plant_height[x, y] = 1
        #plt.figure(1)
        #plt.pcolor(lon, lat, plant_height)
        #plt.title('Site %s' % site)
        # f.variables[var2plot[i]][:, z_pt, x, y]


