import matplotlib.pyplot as plt
import numpy as np
import datetime
import netCDF4
import coawstpy
import pandas as pd
from scipy import stats

# bring in the data
run = 'veg'

if run == 'noveg':
    dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final_noveg'
else:
    dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101_final'

df = pd.read_csv(dir+'/mid_water_currents.csv', index_col=0, parse_dates=True)

#sys.exit()
start_date = '2011-08-01'
end_date = '2011-10-31'

start_date_1 = '2011-08-01'
end_date_1 = '2011-09-06'

start_date_2 = '2011-09-06'
end_date_2 = '2011-09-20'

start_date_3 = '2011-09-20'
end_date_3 = '2011-10-31'

fig, (ax) = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

#a = ax[0].scatter(x=df.loc[start_date:end_date, 'COAWST_Vbar'],y=df.loc[start_date:end_date,'CBIBS_V'],
#                c=df.loc[start_date:end_date].index)#, colormap='viridis')
a = ax[0].scatter(x=df.loc[start_date_1:end_date_1, 'COAWST_Vbar'],y=df.loc[start_date_1:end_date_1,'CBIBS_V'],c='r',s=0.5,label='other times')
a = ax[0].scatter(x=df.loc[start_date_2:end_date_2, 'COAWST_Vbar'],y=df.loc[start_date_2:end_date_2,'CBIBS_V'],c='b',s=0.5,label='09/06 - 09/20')
a = ax[0].scatter(x=df.loc[start_date_3:end_date_3, 'COAWST_Vbar'],y=df.loc[start_date_3:end_date_3,'CBIBS_V'],c='r',s=0.5)
#cm = df.loc['2011-08-01':'2011-09-09', ['COAWST_Ubar','CBIBS_U']].plot.scatter(
#    x='COAWST_Ubar', y='CBIBS_U', c=df.loc['2011-08-01':'2011-09-09'].index, colormap='viridis', ax=ax)
#cbar = fig.colorbar(a, ax=ax[0])
#cbar.ax.set_yticklabels(pd.to_datetime(cbar.get_ticks()).strftime(date_format='%b %d'))
xmin = np.min([df.loc[start_date:end_date, 'COAWST_Vbar'].min(), df.loc[start_date:end_date,'CBIBS_V'].min()])
xmax = np.max([df.loc[start_date:end_date, 'COAWST_Vbar'].max(), df.loc[start_date:end_date,'CBIBS_V'].max()])
ax[0].set_xlim(xmin, xmax)
ax[0].set_ylim(xmin, xmax)
ax[0].grid(axis='both')
ax[0].set_aspect('equal', 'box')
ax[0].set_ylabel('CBIBS')
ax[0].set_xlabel('COAWST')
ax[0].set_title('V - North')
ax[0].legend()
slope, intercept, r_value, p_value, std_err = stats.linregress(
    df.loc[start_date:end_date, 'COAWST_Vbar'].values, df.loc[start_date:end_date,'CBIBS_V'].values)
xs = np.array([xmin, xmax])
ax[0].plot(xs, slope*xs+intercept, linestyle='-', color='k')
ax[0].plot([-2,2],[-2,2],linestyle=':')

#a = ax[1].scatter(x=df.loc[start_date:end_date, 'COAWST_Ubar'],y=df.loc[start_date:end_date,'CBIBS_U'],
#                c=df.loc[start_date:end_date].index)#, colormap='viridis')
a = ax[1].scatter(x=df.loc[start_date_1:end_date_1, 'COAWST_Ubar'],y=df.loc[start_date_1:end_date_1,'CBIBS_U'],c='r',s=0.5,label='other times')
a = ax[1].scatter(x=df.loc[start_date_2:end_date_2, 'COAWST_Ubar'],y=df.loc[start_date_2:end_date_2,'CBIBS_U'],c='b',s=0.5,label='09/06 - 09/20')
a = ax[1].scatter(x=df.loc[start_date_3:end_date_3, 'COAWST_Ubar'],y=df.loc[start_date_3:end_date_3,'CBIBS_U'],c='r',s=0.5)
#cm = df.loc['2011-08-01':'2011-09-09', ['COAWST_Ubar','CBIBS_U']].plot.scatter(
#    x='COAWST_Ubar', y='CBIBS_U', c=df.loc['2011-08-01':'2011-09-09'].index, colormap='viridis', ax=ax)
#cbar = fig.colorbar(a, ax=ax[1])
#cbar.ax.set_yticklabels(pd.to_datetime(cbar.get_ticks()).strftime(date_format='%b %d'))
xmin = np.min([df.loc[start_date:end_date, 'COAWST_Ubar'].min(), df.loc[start_date:end_date,'CBIBS_U'].min()])
xmax = np.max([df.loc[start_date:end_date, 'COAWST_Ubar'].max(), df.loc[start_date:end_date,'CBIBS_U'].max()])
ax[1].set_xlim(xmin, xmax)
ax[1].set_ylim(xmin, xmax)
ax[1].grid(axis='both')
ax[1].set_aspect('equal', 'box')
ax[1].set_ylabel('CBIBS')
ax[1].set_xlabel('COAWST')
ax[1].set_title('U - East')
ax[1].legend()
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
