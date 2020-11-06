import scipy.io as sio
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import coawstpy
import datetime

# Get pts file
#ptsdir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
# ptsdir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/SWAN_20130705_20130715_FRICTION_NOVEG_30SEC_KOMAN_pt4+Bathy'
# ptsfile = ptsdir+"/tripod_wave.pts"
ptsfile = coawstpy.get_file_paths()['ptsfile']

print("Reading SWAN pts file %s" % ptsfile)
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
#obs_file = str('/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/Initialization_data/'
#               'Larry_Flats_data_2013/SF2013JulyData4Matt/sftripod1_advo_diwasp_MKS_LWT_lowpass_results.mat')
obs_file = coawstpy.get_file_paths()['sftripod']
print("Reading tripod data %s" % obs_file)
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
obs_df['mtime'] = obs_df['mtime'].dt.tz_localize('US/Eastern')# obs_data['DWV'].metadata.tbas
obs_df['mtime'] = obs_df['mtime'].dt.tz_convert('UTC')

#swan_hsig = SWAN_df[SWAN_df['Time'] > obs_df['time'].min(), 'Hsig']
swan_hsig = SWAN_df.loc[
    (SWAN_df['Time'] > obs_df['mtime'].min() - datetime.timedelta(minutes=1,seconds=7)) &
    (SWAN_df['Time'] < obs_df['mtime'].max()),
    ['Time','Hsig']]

obs_hsig = obs_df[['mtime','Hsig']]

skill_score = coawstpy.skill_score(swan_hsig['Hsig'], obs_hsig['Hsig'])

