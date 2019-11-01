import os
os.environ["PROJ_LIB"] = "/User/mbiddle/anaconda3/envs/coawst/share/proj/"
import coawstpy
import matplotlib.pyplot as plt
#import matplotlib
#matplotlib.use('MacOSX')


run = 'veg'
point_data_veg = coawstpy.get_point_data(run)
point_data_noveg = coawstpy.get_point_data('noveg')
times = coawstpy.get_time_periods()
vars2plot = ['Pwave_Top', 'Hwave', 'mud_bar', 'sand_bar',
        'river_transport','bed_thickness',
        'depth', 'current_bar','Uwave_rms','bstress_mag','Windv'] #'X-Windv','Y-Windv',
#i=0
#fig, ax = plt.subplots(figsize=(20, 18))#,sharey=True,sharex=True)
#plt.figure(figsize=(20,18))
for site in point_data_veg:
       for event in times:
              if event != 'post-Lee':
                     continue
              if site != 'FLT':
                     continue
              start = coawstpy.nearest_ind(point_data_veg['CBIBS'].index, times[event][0])
              end = coawstpy.nearest_ind(point_data_veg['CBIBS'].index, times[event][1]) + 1
              print('Site: %s %s' % (site, event))
              print("\nveg:\n",point_data_veg[site][vars2plot][start:end].describe().T[['min','max','mean','std']])
              print("\nnoveg:\n", point_data_noveg[site][vars2plot][start:end].describe().T[['min', 'max', 'mean','std']])
              ax = point_data_veg[site][vars2plot][start:end].plot(subplots=True, sharex=True, linewidth=0.5,
                                                                   legend=False, linestyle='-')
              point_data_noveg[site][vars2plot][start:end].plot(ax=ax, subplots=True, sharex=True, linestyle='--',
                                                                linewidth=0.5, legend=False)
              ax[0].legend(['veg','noveg'],loc=3,ncol=2, bbox_to_anchor=(.25,1.02,.5,.102), mode='expand',borderaxespad=0)
              plt.suptitle('%s %s' % (event, site))
              evnt=event
              st=site
              i=0
              for v in vars2plot:
                     ax[i].set_ylabel(v,rotation=0,labelpad=40)
                     ax[i].yaxis.set_label_position("right")
                     i+=1
              print('')
              #writedir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/for_larry_20191007/'
              #writevegfile = writedir+"%s_%s_veg.csv" % (site,event)
              # writedir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/figures/timeseries/test/'
              # image_name = '%s_%s_timeseries.png' % (event, site)
              # outfile = writedir+image_name
              # print("Saving image to %s" % outfile)
              # plt.savefig(outfile, bbox_inches='tight', dpi=500)
              #writenovegfile = writedir+"%s_%s_noveg.csv" % (site,event)
              #point_data_veg[site][vars2plot][start:end].to_csv(writevegfile)
              #point_data_noveg[site][vars2plot][start:end].to_csv(writenovegfile)

## Put figure in the Windv axis #6
ticks = ax[10].get_xticklabels()
ax[10].clear()
q = coawstpy.stick_plot(point_data_veg[site].index[start:end],point_data_veg[site]['X-Windv'][start:end],
                    point_data_veg[site]['Y-Windv'][start:end], ax=ax[10], scale=500)
ref = 10
qk = ax[10].quiverkey(q, -0.04, 0.15, ref,
                  "%s m/s" % ref,
                  labelpos='N', coordinates='axes', fontproperties={'size': 'xx-small'})
ax[10].set_ylabel('wind',rotation=0,labelpad=40)
for label in ax[10].get_xmajorticklabels():
    label.set_rotation(30)
    label.set_horizontalalignment("right")
#ax[10].tick_params(axis='x', rotation=45, labelright=True)
#ticks = ax[10].get_xticklabels()
#ax[6].set_xticklabels(ticks)#x[7].get_xticklabels())
#ax[10].set_xticklabels(ticks)
#for event in times:
#    start = coawstpy.nearest_ind(point_data['CBIBS'].index, times[event][0])
#    end = coawstpy.nearest_ind(point_data['CBIBS'].index, times[event][1]) + 1
#    for site in point_data:
#        point_data[site][vars2plot][start:end].plot()
#        plt.suptitle('%s %s' % (event, site))

writedir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/figures/timeseries/'
image_name = '%s_%s_timeseries.png' % (evnt, st)
outfile = writedir+image_name
print("Saving image to %s" % outfile)
plt.savefig(outfile, bbox_inches='tight', dpi=500)