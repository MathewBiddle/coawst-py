import netCDF4
import pandas as pd
import matplotlib.pyplot as plt
import coawstpy
import matplotlib.dates as mdates
import numpy as np
import datetime

'''
This script recreates the figure in Gurbisz et al 2015, adding model data
Citation:
Gurbisz, C., Kemp, W.M., Sanford, L.P., Orth, R.J., 2016. Mechanisms of Storm-Related Loss and Resilience in a 
Large Submersed Plant Bed. Estuaries and Coasts 39, 951–966. https://doi.org/10.1007/s12237-016-0074-4
'''
dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'

# Get SSC data
inputfile = dir + '/upper_ches_his.nc'
f = netCDF4.Dataset(inputfile, 'r')
print("Retrieving %s" % (inputfile.split("/")[-1]))
ocean_time = f.variables['ocean_time'][:]
datetime_list = []
for sec in ocean_time:
    if sec == 0.0:
        datetime_list.append(
            netCDF4.num2date(sec + 0.0000000000000001, units=f.variables['ocean_time'].units,
                             calendar=f.variables['ocean_time'].calendar))
    else:
        datetime_list.append(
            netCDF4.num2date(sec, units=f.variables['ocean_time'].units,
                             calendar=f.variables['ocean_time'].calendar))
lat = f.variables['lat_rho'][:]
lon = f.variables['lon_rho'][:]
point_data = dict()
locs = coawstpy.get_point_locations()
#sites = ['CBIBS', '3', 'SUS', 'FLT']
for site in locs['Site']:
    point_data[site] = pd.DataFrame()
    lat_pt, lon_pt = locs.loc[locs['Site'] == site, ['lat', 'lon']].values[0]
    x = np.abs(lon[:, 1] - lon_pt).argmin()
    y = np.abs(lat[1, :] - lat_pt).argmin()
    point_data[site]['mud_bar'] = np.average(f.variables['mud_01'][:, :, x, y], axis=1)
    point_data[site]['sand_bar'] = np.average(f.variables['sand_01'][:, :, x, y], axis=1)


## river data
print("Reading river data...")
river_frc = dir+'/river_frc.nc'
f_river = netCDF4.Dataset(river_frc, 'r')
river_time = f_river.variables['river_time'][:]
river_transport = f_river.variables['river_transport'][:, 0]
river_datetime_list=[]
for sec in river_time:
    river_datetime_list.append(
        netCDF4.num2date(sec, units=f_river.variables['river_time'].units, calendar='standard'))

## wind data
print("Reading wind data...")
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

# ## tide data
# print("Reading tide data...")
# bry_file = dir + '/upper_ches_bry.nc'
# f_bry = netCDF4.Dataset(bry_file, 'r')
# bry_time = f_bry.variables['zeta_time'][:]
# bry_zeta = f_bry.variables['zeta_south'][:, 50]
# bry_datetime_list=[]
# for sec in bry_time:
#     bry_datetime_list.append(netCDF4.num2date(sec, units=f_bry.variables['zeta_time'].units, calendar='standard'))

# EOTB data
print('Reading EOTB data...')
eotb_file = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/Initialization_data/Eyes_on_the_bay/EOTBData_HavredeGrace_Flats_01Jul11_TO_01Nov11.csv'
feotb = pd.read_csv(eotb_file, header=0, parse_dates=['DateTime'], infer_datetime_format=True)
# eotb_turb = feotb['Turb_NTU'].values
tshift = pd.DateOffset(hours=0)
feotb['DateTimeUTC'] = feotb['DateTime']-tshift
# get only the data from out simulation
eotb_FLT_turb = feotb.loc[(feotb['DateTimeUTC'] < datetime_list[-1]) & (feotb['DateTimeUTC'] > datetime_list[0]) &
                      (feotb[' Station'] == ' Chesapeake Bay Segment 1 - Susquehanna Flats'), [' Turb_NTU']].values
eotb_SUS_turb = feotb.loc[(feotb['DateTimeUTC'] < datetime_list[-1]) & (feotb['DateTimeUTC'] > datetime_list[0]) &
                      (feotb[' Station'] == ' Susquehanna River - Havre de Grace'), [' Turb_NTU']].values
eotb_SUS_date = feotb.loc[(feotb['DateTimeUTC'] < datetime_list[-1]) & (feotb['DateTimeUTC'] > datetime_list[0]) &
                      (feotb[' Station'] == ' Susquehanna River - Havre de Grace'), ['DateTimeUTC']].values
eotb_FLT_date = feotb.loc[(feotb['DateTimeUTC'] < datetime_list[-1]) & (feotb['DateTimeUTC'] > datetime_list[0]) &
                      (feotb[' Station'] == ' Chesapeake Bay Segment 1 - Susquehanna Flats'), ['DateTimeUTC']].values
eotb_lon = -76.0848
eotb_lat = 39.5478
eotb_depth = 1.0


## plotting
print("Creating plots...")
fig, (ax) = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(10, 6))
fig.subplots_adjust(hspace=0.1)
myFmt = mdates.DateFormatter("%b")
months = mdates.MonthLocator()  # every month

dayint=10
xlim= (ptsdf['Time'].min(), ptsdf['Time'].max())

##
# time_periods = coawstpy.get_time_periods()
# # add verical bars:
# for i in ax:
#     for time_period in time_periods:
#         i.axvspan(time_periods[time_period][0],time_periods[time_period][1], facecolor='0.5', alpha=0.3)  # Typical low-flow conditions
#     #i.axvspan('2011-08-01','2011-08-06', facecolor='0.5', alpha=0.3) # Typical low-flow conditions
    #i.axvspan('2011-08-27','2011-08-30', facecolor='0.5', alpha=0.3) # Irene wind event
    #i.axvspan('2011-09-07', '2011-09-16', facecolor='0.5', alpha=0.3) # Lee discharge
    #i.axvspan('2011-10-13', '2011-10-24', facecolor='0.5', alpha=0.3)  # End wind event


# ax[0].plot_date(river_datetime_list,
#                 np.sum(f_river.variables['river_sand_01'],axis=(1,2)), label='sand',
#                 xdate=True, linestyle='-', linewidth=0.5,
#                 marker='', markersize=1)
# ax[0].plot_date(river_datetime_list,
#                 np.sum(f_river.variables['river_mud_01'],axis=(1,2)), label='mud',
#                 xdate=True, linestyle='-', linewidth=0.5,
#                 marker='', markersize=1)
# ax[0].xaxis.set_major_locator(months)
# #ax[0].xaxis.set_major_formatter(myFmt)
# ax[0].set_xlim(xlim)
# ax[0].legend()
# ax[0].set_ylabel('River SSC (kg/m3)')


ax[1].plot_date(river_datetime_list,
                (river_transport + (0.2 * river_transport)) * f_river.variables['river_transport'].shape[1],
                xdate=True, linestyle='-', linewidth=2,
                marker='', markersize=1, c='black')
ax[1].xaxis.set_major_locator(months)
#ax[1].xaxis.set_major_formatter(myFmt)
ax[1].set_xlim(xlim)
ax[1].set_ylabel('$Q$ ($m^3$ $s^{-1}$)')
ax[1].set_yticks([0,8000,16000])
ax[1].set_yticklabels([0,8000,16000], va='center', rotation=90)
ax[1].patch.set_facecolor('white')
#ax[0].set_yticks(rotation='horizontal')

q = coawstpy.stick_plot(ptsdf['Time'],ptsdf['X-Windv'],ptsdf['Y-Windv'], ax=ax[2],scale=200)
ax[2].xaxis.set_major_locator(months)
#ax[2].xaxis.set_major_formatter(myFmt)
ax[2].set_xlim(xlim)
ax[2].set_ylabel('$U_{10}$',labelpad=20)
ref = 10
qk = ax[2].quiverkey(q, 0.04, 0.1, ref,
                  "%s m/s" % ref,
                  labelpos='N', coordinates='axes', fontproperties={'size': 'xx-small'})
ax[2].patch.set_facecolor('white')
#ax[2].plot_date(ptsdf['Time'], np.sqrt(ptsdf['X-Windv']**2 + ptsdf['Y-Windv']**2),
#                xdate=True, linestyle='-', linewidth=0.5, marker='', markersize=1)
#ax[2].xaxis.set_major_locator(mdates.DayLocator(interval=dayint))
#ax[2].xaxis.set_major_formatter(myFmt)
#ax[2].set_xlim(xlim)
ax[2].set_xlabel('2011',fontdict=dict(size='large'))


# ax[3].plot_date(bry_datetime_list,bry_zeta, xdate=True, linestyle='-', linewidth=0.5,
#                      marker='', markersize=1)
# ax[3].set_xlim(xlim)
# ax[3].set_ylabel('Water surface [m]')

# ax[0].plot_date(datetime_list,point_data['SUS']['mud_bar']+point_data['SUS']['sand_bar'],label='SUS',
#                 linestyle='-', linewidth=0.5, marker='', markersize=1, c='black')
# ax[0].plot_date(datetime_list,point_data['FLT']['mud_bar']+point_data['FLT']['sand_bar'],label='FLT',
#                 linestyle='-', linewidth=0.5, marker='', markersize=1, c='grey')
# #ax[0].legend(loc='upper left')
# #ax[4].set_yscale('log')
# ax[0].set_ylabel('$\\widebar{SSC}_{tot}$')
# ax[0].xaxis.set_major_locator(months)
# #ax[4].xaxis.set_major_formatter(myFmt)

# eyes on the bay data
ax[0].plot_date(eotb_SUS_date,eotb_SUS_turb,label='EOTB_SUS',
                linestyle='-', linewidth=2, marker='', markersize=1, c='blue')
ax[0].plot_date(eotb_FLT_date,eotb_FLT_turb,label='EOTB_FLT',
                linestyle='-', linewidth=2, marker='', markersize=1, c='green')
#ax[1].legend(loc='upper left')
#ax[4].set_yscale('log')
ax[0].set_ylim([0,650])
ax[0].set_ylabel('$Turbidity$ $(NTU)$')
ax[0].set_yticks([0,200,400,600])
ax[0].set_yticklabels([0,200,400,600], va='center', rotation=90)
ax[0].xaxis.set_major_locator(months)
ax[0].xaxis.set_major_formatter(myFmt)
ax[0].patch.set_facecolor('white')

ax[0].set_title('Effect of Wind on Turbidity')
# ax[0].text('2011-07-30 00:00',2.5,'Before Irene')
# ax[0].text('2011-08-26 12:00',2.50000,'Irene')
# ax[0].text('2011-09-10',2.50000,'Lee')
# ax[0].text('2011-10-15',2.50000,'post-Lee')

# ax[0].text('2011-07-21',2.6,'a',fontsize=18)
# ax[1].text('2011-07-21',525,'b',fontsize=18)
# ax[2].text('2011-07-21',18000,'c',fontsize=18)
# ax[3].text('2011-07-21',0.035,'d',fontsize=18)

#ax[4].set_xlim(time_periods['post-Lee'][0],time_periods['post-Lee'][1])
#outfile = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/defense/Gurbisz_2016_turbidity_fig_recreation_obs.png'
#plt.savefig(outfile, bbox_inches='tight', dpi = 500, facecolor=ax[0].get_facecolor(), edgecolor='none')#transparent=True)