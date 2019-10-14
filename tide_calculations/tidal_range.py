from scipy.signal import find_peaks
import pandas as pd
from pytides.tide import Tide
import numpy as np
import netCDF4
import coawstpy
import matplotlib.pyplot as plt
import scipy.integrate as integrate

point_data = coawstpy.get_point_data('veg')


min_peaks, _ = find_peaks(point_data['CBIBS']['depth'] * -1, height=-6)
max_peaks, _ = find_peaks(point_data['CBIBS']['depth'], height=0)

#plt.plot(point_data['CBIBS']['depth'])
#plt.plot(point_data['CBIBS']['depth'][max_peaks], "x")
#plt.plot(point_data['CBIBS']['depth'][min_peaks], 'o')
## calculate the tidal amplitude between high and low tides.
amplitude = []
for i in range(0,len(min_peaks)):
    amplitude.append(point_data['CBIBS']['depth'][max_peaks[i]] - point_data['CBIBS']['depth'][min_peaks[i]])

#point_data['CBIBS']['amplitude'] = amplitude