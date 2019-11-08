import scipy.io as sio
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import coawstpy
import datetime

# Get pts file
#ptsdir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
ptsdir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/SWAN_20130705_20130715_FRICTION_NOVEG_30SEC_KOMAN_pt4+Bathy'
ptsfile = ptsdir+"/tripod_wave.pts"

SWAN_df = pd.read_fwf(ptsfile, header=4)#,widths=[18,16,16,16,16,16,16,16,16])#,skiprows=range(0,5))

SWAN_df.drop([0,1],axis=0,inplace=True)
SWAN_df.rename(columns={'%       Time':'Time'},inplace=True)
SWAN_df['Yp'] = SWAN_df['Yp            Hsig'].astype(str).str.split("    ",expand=True)[0].astype(float)
SWAN_df['Hsig'] = SWAN_df['Yp            Hsig'].astype(str).str.split("    ",expand=True)[1].astype(float)
SWAN_df.drop(columns=['Yp            Hsig'],inplace=True)

SWAN_df['Time']=pd.to_datetime(SWAN_df['Time'],format='%Y%m%d.%H%M%S',utc=True)
SWAN_df['Time']=SWAN_df['Time']
SWAN_df['Hsig']=SWAN_df['Hsig'].astype(float)
SWAN_df['X-Windv']=SWAN_df['X-Windv'].astype(float)
SWAN_df['Y-Windv']=SWAN_df['Y-Windv'].astype(float)

# Get observational data
obs_file = str('/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/Initialization_data/'
               'Larry_Flats_data_2013/SF2013JulyData4Matt/sftripod1_advo_diwasp_MKS_LWT_lowpass_results.mat')

obs_data = sio.loadmat(obs_file, struct_as_record=False, squeeze_me=True)

obs_df = pd.DataFrame(columns=obs_data['DWV']._fieldnames)

for var in obs_data['DWV']._fieldnames:
    if var == 'metadata':
        exec("metadata=obs_data['DWV'].%s" % var)
        continue
    exec("obs_df['%s']=obs_data['DWV'].%s" % (var, var))

obs_df.drop(columns='metadata', inplace=True)
obs_df['Tp'].loc[obs_df['Tp'] > 20] = np.nan

obs_df['mtime'] = pd.to_datetime(obs_df['mtime']-719529, unit='D')
#obs_df['mtime']

#obs_df = obs_df.set_index('mtime')
obs_df['mtime']=obs_df['mtime'].dt.tz_localize('US/Eastern')# obs_data['DWV'].metadata.tbase


#obs_df.plot(kind='line', x='mtime', y=['Hsig', 'Tp'])

fig, (ax) = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(12, 8))
ax[0].plot_date(obs_df['mtime'], obs_df['Hsig'], label='Observed', xdate=True, linestyle='', linewidth=0.5, c='grey',
                     marker='.', markersize=1)
ax[0].plot_date(SWAN_df['Time'], SWAN_df['Hsig'],label='Predicted', xdate=True, linestyle='-', linewidth=0.5, c='k',
                     marker='', markersize=1)
#ax[0].legend()
ax[0].set_ylabel('$H_s$ (m)')
ax[0].text(datetime.datetime(2013,7,4,20,00),0.32,'a',fontdict=dict(weight="bold"),bbox=dict(fc="none",lw=2))

ax[1].plot_date(obs_df['mtime'], obs_df['Tp'], label='Observed', xdate=True, linestyle='', linewidth=0.5, c='grey',
                     marker='.', markersize=1)
ax[1].plot_date(SWAN_df['Time'], SWAN_df['RTpeak'].astype(float), label='Predicted', xdate=True, linestyle='-',
                     linewidth=0.5, c='k', marker='', markersize=1)
ax[1].legend()
ax[1].set_ylim(0, 6)
ax[1].set_ylabel('$T_p$ (s)')
ax[1].text(datetime.datetime(2013,7,4,20,00),5,'b',fontdict=dict(weight="bold"),bbox=dict(fc="none",lw=2))

q = coawstpy.stick_plot(SWAN_df['Time'],SWAN_df['X-Windv'],SWAN_df['Y-Windv'], ax=ax[2])
ax[2].set_ylabel('$U_{10}$')
ref = 10
ax[2].quiverkey(q, 0.06, 0.15, ref,
                  "%s m/s" % ref,
                  labelpos='N', coordinates='axes', fontproperties={'size': 'medium'})

ax[2].text(datetime.datetime(2013,7,4,20,00),0.04,'c',fontdict=dict(weight="bold"),bbox=dict(fc="none",lw=2))
#plt.xlabel(plt.rcParams['timezone'])
#triplon = metadata.location[1]
#triplon = -1*(float(triplon.split(" ")[0])+(float(triplon.split(" ")[2])/60))
#triplat = metadata.location[0]
#triplat = float(triplat.split(" ")[0])+(float(triplat.split(" ")[1])/60)
#plt.suptitle("Tripod @ %.4f %.3f\nSWAN  @ %s %s" % (triplat, triplon, SWAN_df['Yp'].unique()[0],SWAN_df['Xp'].unique()[0]))
plt.subplots_adjust(hspace=0.05)
outfile = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/figures/2013_SWAN_wave_comparison.png'
plt.savefig(outfile, bbox_inches='tight', dpi = 500)