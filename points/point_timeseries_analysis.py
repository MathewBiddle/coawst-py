import os
os.environ["PROJ_LIB"] = "/User/mbiddle/anaconda3/envs/coawst/share/proj/"
import coawstpy
import matplotlib.pyplot as plt
#import matplotlib
#matplotlib.use('MacOSX')
import numpy as np


run = 'veg'
point_data_veg = coawstpy.get_point_data(run)
point_data_noveg = coawstpy.get_point_data('noveg')
times = coawstpy.get_time_periods()
#vars2plot = ['Pwave_Top', 'Hwave', 'mud_bar', 'sand_bar',
#        'river_transport','bed_thickness',
#        'depth', 'current_bar','Uwave_rms','bstress_mag','Windv'] #'X-Windv','Y-Windv',
vars2plot = ['river_transport','mud_bar','sand_bar','current_bar','bstress_mag']
#i=0
#fig, ax = plt.subplots(figsize=(20, 18))#,sharey=True,sharex=True)
#plt.figure(figsize=(20,18))
#labels =['$T_p$','$H_s$','$\\widebar{SSC_{f}}$','$\\widebar{SSC_{c}}$','$Q$','$h_b$','$d$','$|\\widebar{\\upsilon}|$',
#         '$\\upsilon_{bo}$','$|\\tau_b|$','$U_{10}$']
labels = ['$Q$ ($m^3$ $s^{-1}$)','$\\widebar{SSC}_{f}$ (kg $m^{-3}$)','$\\widebar{SSC}_{c}$ (kg $m^{-3}$)',
          '$|\\widebar{\\upsilon}|$ (m $s^{-1}$)','$|\\tau_b|$ (N $m^{-2}$)']
alpha = 'a'
for site in point_data_veg:
       for event in times:
              if event != 'typical':
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
              ax[0].legend(['veg','noveg'],loc=3,ncol=2, bbox_to_anchor=(.25,1.05,.5,.102), mode='expand',borderaxespad=0)
              plt.suptitle('%s %s' % (event, site))
              evnt=event
              st=site
              i=0
              for v in vars2plot:

                  #ax[i].set_ylabel(v,rotation=0,labelpad=40)
                  #ax[i].yaxis.set_label_position("right")
                  ax[i].set_ylabel(labels[i],size=6)
                  #ax[i].set_ylabel(alpha,rotation=0,labelpad=5)
                  ax[i].tick_params(labelsize=6)
                  alpha = chr(ord(alpha) + 1)
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
# ticks = ax[10].get_xticklabels()
# ax[10].clear()
# q = coawstpy.stick_plot(point_data_veg[site].index[start:end],point_data_veg[site]['X-Windv'][start:end],
#                     point_data_veg[site]['Y-Windv'][start:end], ax=ax[10], scale=500)
# ref = 10
# qk = ax[10].quiverkey(q, -0.04, 0.15, ref,
#                   "%s m/s" % ref,
#                   labelpos='N', coordinates='axes', fontproperties={'size': 'xx-small'})
# #ax[10].set_ylabel('k')
# ax[10].set_ylabel(labels[i-1])
# #ax[10].set_ylabel('wind',rotation=0,labelpad=40)
# for label in ax[10].get_xmajorticklabels():
#     label.set_rotation(0)
#     label.set_horizontalalignment("center")

# remove the exponent for bed thickness
#ax[5].get_yaxis().get_major_formatter().set_useOffset(False)

# add critical shear stress line
ax[4].plot_date(point_data_veg[site][vars2plot][start:end].index,np.ones(end-start)*0.049,':',linewidth=0.5,c='black')
# squash all x-axis together
plt.subplots_adjust(hspace=0.2)
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

#writedir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/figures/timeseries/'
writedir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/Paper/defense/'
image_name = '%s_%s_timeseries.png' % (evnt, st)
outfile = writedir+image_name
print("Saving image to %s" % outfile)
plt.savefig(outfile, bbox_inches='tight', dpi=500)