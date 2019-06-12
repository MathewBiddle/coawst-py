import scipy.io as sio
import numpy as np
import pandas as pd

# Get observational data
obs_file = str('/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/Initialization_data/'
               'Larry_Flats_data_2013/SF2013JulyData4Matt/sftripod1_advo_diwasp_MKS_LWT_lowpass_results.mat')

obs_data = sio.loadmat(obs_file, struct_as_record=False, squeeze_me=True)

ldf = pd.DataFrame(columns=obs_data['DWV']._fieldnames)

for var in obs_data['DWV']._fieldnames:
    if var == 'metadata':
        exec("metadata=obs_data['DWV'].%s" % var)
        continue
    exec("ldf['%s']=obs_data['DWV'].%s" % (var, var))

ldf.drop(columns='metadata', inplace=True)
ldf['Tp'].loc[ldf['Tp'] > 20] = np.nan

ldf['mtime'] = pd.to_datetime(ldf['mtime']-719529, unit='D', utc=True)

ldf = ldf.set_index('mtime')

ldf.plot(kind='line', y=['Hsig', 'Tp'])
