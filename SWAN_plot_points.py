import pandas as pd

file='/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/SWAN_20130705_20130715_FRICTION_NOVEG_30SEC_JANSSEN_pt4+Bathy/tripod_wave.pts'

df = pd.read_fwf(file, header=4,widths=[18,16,16,16,16,16,16,16])#,skiprows=range(0,5))

df.drop([0,1],axis=0,inplace=True)
df.rename(columns={'%       Time':'Time'},inplace=True)
df['Time']=pd.to_datetime(df['Time'],format='%Y%m%d.%H%M%S',utc=True)
df['Hsig']=df['Hsig'].astype(float)
df.plot(x='Time',y='Hsig')