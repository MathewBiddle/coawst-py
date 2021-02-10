import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import matplotlib.animation as animation
import coawstpy

# bring in the data
files = coawstpy.get_file_paths()
#dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/Full_20110719T23_20111101'
#inputfile = dir+'/upper_ches_his.nc'
f = netCDF4.Dataset(files['veg'], 'r')

ocean_time = f.variables['ocean_time'][:]
## Do some date conversions ##
datetime_list=[]
for sec in ocean_time:
    datetime_list.append(
        netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))

# figure out the transect:

# identify the points to extract data
# Susquehanna River mouth
#trans_name = 'S.R. Mouth'
#x = np.array([29,30,31,32])
#y = np.array([13,12,11,10])
#z = 2 # 4 surface, 0 bottom

# Turkey Point to Sandy Point
trans_name = 'Turkey Pt to Sandy Pt'
x = np.array(list(range(42,67)))  #
y = np.array([58]*len(x))
t = 100
lat = f.variables['lat_rho'][x,y]
lon = f.variables['lon_rho'][x,y]
depth = f.variables['z_rho'][:]

u = f.variables['u_eastward'][:]
#v = f.variables['v_northward'][:]

v = f.variables['mud_01'][:]

fig, ax = plt.subplots(figsize=(8, 6))
cax=plt.pcolormesh(lon[:,0],depth[0,:-1,x,58].T,v[0,:-1,x,58].T,vmin=0, vmax=0.5, cmap='jet')
cbar = fig.colorbar(cax)
def animate(i):
    # pcolor variable
    cax.set_array(v[i, :-1, x,58].T.flatten())
    plt.title('%s' % (datetime_list[i]))

anim = animation.FuncAnimation(fig, animate, frames=len(datetime_list))
