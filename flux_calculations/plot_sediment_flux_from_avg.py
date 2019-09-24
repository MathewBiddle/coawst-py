import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import coawstpy
import scipy.integrate as integrate
'''
Plots the sediment flux across transects.
'''

transect_indexes = coawstpy.get_transect_indexes()
times = coawstpy.get_time_periods()
transect_data = coawstpy.get_sed_flux_data('veg',transect_indexes)

sediment_type = ['mud','sand']

for event in times:
    start = coawstpy.nearest_ind(transect_data['time'],times[event][0])
    end = coawstpy.nearest_ind(transect_data['time'],times[event][1])+1
    #fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
    for t in transect_data:
        if t == 'time':
            continue
        for sed in sediment_type:
            #print('%s %s max south flux = %g %s' % (t, sed, transect_data[t]['Huon_sum_%s' % sed][start:end].max(),
            #                                        transect_data[t]['Huon_sum_%s_units' % sed]))

            print('%s %s %s integrated south flux = %i tons' % (event, t, sed,
                                    (integrate.trapz(transect_data[t]['Huon_sum_%s' % sed][start:end]) * 3600)/1000))

            # ax.plot_date(transect_data['time'][start:end], transect_data[t]['Huon_sum_%s' % sed][start:end],
            #              label="%s %s %s" % (t,sed, event), linestyle='-', linewidth=1, marker='')
            # ax.legend()
            # ax.xaxis.set_major_locator(mdates.DayLocator(interval=10))
            # ax.xaxis.set_major_formatter(DateFormatter("%m/%d"))
            # ax.set_ylabel('Total south flux %s' % transect_data[t]['Huon_sum_%s_units' % sed])

    T1_sand_integral = (integrate.trapz(transect_data['T1']['Huon_sum_sand'][start:end]) * 3600)/1000
    T1_mud_integral = (integrate.trapz(transect_data['T1']['Huon_sum_mud'][start:end]) * 3600)/1000

    T2_sand_integral = (integrate.trapz(transect_data['T2']['Huon_sum_sand'][start:end]) * 3600)/1000
    T2_mud_integral = (integrate.trapz(transect_data['T2']['Huon_sum_mud'][start:end]) * 3600)/1000

    print('%s total sediment through T1: %i tons' % (event,T1_sand_integral+T1_mud_integral))
    print('%s total sediment through T2: %i tons\n' % (event, T2_sand_integral+T2_mud_integral))

print('Total mud through T1: %i tons' % ((integrate.trapz(transect_data['T1']['Huon_sum_mud'][:]) * 3600)/1000))
print('Total sand through T1: %i tons' % ((integrate.trapz(transect_data['T1']['Huon_sum_sand'][:]) * 3600)/1000))
print('Total sum through T1: %i tons\n' % (((integrate.trapz(transect_data['T1']['Huon_sum_sand'][:]) * 3600)/1000)+
      (integrate.trapz(transect_data['T1']['Huon_sum_mud'][:]) * 3600)/1000))

print('Total mud through T2: %i tons' % ((integrate.trapz(transect_data['T2']['Huon_sum_mud'][:]) * 3600)/1000))
print('Total sand through T2: %i tons' % ((integrate.trapz(transect_data['T2']['Huon_sum_sand'][:]) * 3600)/1000))
print('Total sum through T2: %i tons' % (((integrate.trapz(transect_data['T2']['Huon_sum_sand'][:]) * 3600)/1000)+
      (integrate.trapz(transect_data['T2']['Huon_sum_mud'][:]) * 3600)/1000))

