import os
os.environ["PROJ_LIB"] = "/anaconda3/envs/coawst/share/proj/"
import coawstpy
import matplotlib.pyplot as plt

run = 'veg'
point_data_veg = coawstpy.get_point_data(run)
point_data_noveg = coawstpy.get_point_data('noveg')
times = coawstpy.get_time_periods()
vars2plot = ['Pwave_Top', 'Hwave', 'mud_bar', 'sand_bar',
       'bed_thickness', 'river_transport',
       'Windv', 'depth', 'current_bar']
i=0
for site in point_data_veg:
       for event in times:
              if event != 'Lee':
                     continue
              if site != 'S':
                     continue
              start = coawstpy.nearest_ind(point_data_veg['CBIBS'].index, times[event][0])
              end = coawstpy.nearest_ind(point_data_veg['CBIBS'].index, times[event][1]) + 1
              print('Site: %s %s' % (site, event))
              print(point_data_veg[site][vars2plot][start:end].describe().T[['min','max','mean']])
              ax = point_data_veg[site][vars2plot][start:end].plot(subplots=True,sharex=True,linewidth=0.5,legend=False,linestyle='-')
              point_data_noveg[site][vars2plot][start:end].plot(ax=ax, subplots=True, sharex=True,linestyle='--',linewidth=0.5,legend=False)
              ax[0].legend(['veg','noveg'],loc=3,ncol=2, bbox_to_anchor=(.25,1.02,.5,.102), mode='expand',borderaxespad=0)
              plt.suptitle('%s %s' % (event, site))
              for v in vars2plot:
                     ax[i].set_ylabel(v,rotation=0,labelpad=20)
                     ax[i].yaxis.set_label_position("right")
                     i+=1
              print('')

#for event in times:
#    start = coawstpy.nearest_ind(point_data['CBIBS'].index, times[event][0])
#    end = coawstpy.nearest_ind(point_data['CBIBS'].index, times[event][1]) + 1
#    for site in point_data:
#        point_data[site][vars2plot][start:end].plot()
#        plt.suptitle('%s %s' % (event, site))
