import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
import coawstpy
import matplotlib.pyplot as plt

run = 'noveg'
point_data = coawstpy.get_point_data(run)
times = coawstpy.get_time_periods()
vars2plot = ['Pwave_Top', 'Hwave', 'mud_bar', 'sand_bar',
       'bed_thickness', 'river_transport',
       'Windv', 'depth', 'current_bar']

for site in point_data:
       for event in times:
              #if event != 'post-Lee':
              #       continue
              start = coawstpy.nearest_ind(point_data['CBIBS'].index, times[event][0])
              end = coawstpy.nearest_ind(point_data['CBIBS'].index, times[event][1]) + 1
              print('%s Site: %s %s' % (run, site, event))
              print(point_data[site][vars2plot][start:end].describe().T[['min','max','mean']])
              print('')

#for event in times:
#    start = coawstpy.nearest_ind(point_data['CBIBS'].index, times[event][0])
#    end = coawstpy.nearest_ind(point_data['CBIBS'].index, times[event][1]) + 1
#    for site in point_data:
#        point_data[site][vars2plot][start:end].plot()
#        plt.suptitle('%s %s' % (event, site))
