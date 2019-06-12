import pandas as pd

dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
file = dir+"/tripod_wave.pts"
#file='/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/SWAN_20130705_20130715_FRICTION_NOVEG_30SEC_JANSSEN_pt4+Bathy/tripod_wave.pts'


df = pd.read_fwf(file, header=4)#,widths=[18,16,16,16,16,16,16,16,16])#,skiprows=range(0,5))

df.drop([0,1],axis=0,inplace=True)
df.rename(columns={'%       Time':'Time'},inplace=True)
df['Yp'] = df['Yp            Hsig'].astype(str).str.split("    ",expand=True)[0].astype(float)
df['Hsig'] = df['Yp            Hsig'].astype(str).str.split("    ",expand=True)[1].astype(float)
df.drop(columns=['Yp            Hsig'],inplace=True)

df['Time']=pd.to_datetime(df['Time'],format='%Y%m%d.%H%M%S',utc=True)
df['Hsig']=df['Hsig'].astype(float)
df['X-Windv']=df['X-Windv'].astype(float)
df['Y-Windv']=df['Y-Windv'].astype(float)
#df.plot(x='Time',y='Hsig')