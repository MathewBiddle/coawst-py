import pandas as pd
import matplotlib.pyplot as plt

url = 'https://tidesandcurrents.noaa.gov/api/datagetter?product=air_temperature&application=NOS.COOPS.TAC.MET&begin_date=20181121&end_date=20191120&station=8638901&time_zone=GMT&units=english&interval=h&format=csv'

df = pd.read_csv(url)#,header=[0,1],sep=' ')

#df['WDIR WSPD GST'].split(' ',expand=True)[1]