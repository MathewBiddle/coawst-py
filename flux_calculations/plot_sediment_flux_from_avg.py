import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import coawstpy

'''
Plots the sediment flux across transects.
'''

transect_indexes = coawstpy.get_transect_indexes()

transect_data = coawstpy.get_sed_flux_data('veg',transect_indexes)
sediment_type = ['mud','sand']
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
for t in transect_data:
    if t == 'time':
        continue
    for sed in sediment_type:
        print('%s %s max south flux = %g %s' % (t, sed, transect_data[t]['Huon_sum_%s' % sed].max(),
                                                transect_data[t]['Huon_sum_%s_units' % sed]))
        print('%s %s integrated south flux = %g %s\n' % (t, sed, transect_data[t]['Huon_%s_integrated' % sed],
                                                         transect_data[t]['Huon_%s_integrated_units' % sed]))

        ax.plot_date(transect_data['time'], transect_data[t]['Huon_sum_%s' % sed], label="%s %s" % (t,sed),
                      linestyle='-', linewidth=1, marker='')
        ax.legend()
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=10))
        ax.xaxis.set_major_formatter(DateFormatter("%m/%d"))
        ax.set_ylabel('Total south flux %s' % transect_data[t]['Huon_sum_%s_units' % sed])

print('Total sediment through T1: %e tons' %
      (transect_data['T1']['Huon_mud_integrated']+transect_data['T1']['Huon_sand_integrated']))
print('Total sediment through T2: %e tons' %
      (transect_data['T2']['Huon_mud_integrated']+transect_data['T2']['Huon_sand_integrated']))

