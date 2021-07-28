import coawstpy
import matplotlib.pyplot as plt
#import matplotlib
#matplotlib.use('MacOSX')
import numpy as np
import pandas as pd


#run = 'veg'
print("Collecting data...")
sites = ['1', '2', '3', '4', '5']
point_data_veg = coawstpy.get_point_data('veg', sites)
point_data_noveg = coawstpy.get_point_data('noveg', sites)
#times = coawstpy.get_time_periods()

out_veg = pd.DataFrame()
out_noveg = pd.DataFrame()
for site in sites:

    columns = ['Latitude', 'Longitude', 'bed_thickness', 'sandmass_sum', 'mudmass_sum']
    # in_veg = point_data_veg[site][columns]
    # in_veg['site'] = site
    # out_veg = pd.concat([out_veg, in_veg], axis=0)

    fname = 'data/site_%s.csv' % site
    point_data_veg[site].to_csv(fname,
                                columns=columns,
                                header=['Latitude', 'Longitude', 'bed_thickness (m)', 'sandmass_sum (kg m^-2)',
                                        'mudmass_sum (kgm^-2)'],
                                index_label=['ISO Datetime'],
                                date_format='%Y-%m-%dT%H:%M:%SZ',
                                )

    # fname = 'data/site_%s_noveg.csv' % site
    # point_data_noveg[site].to_csv(fname,
    #                               columns=['Latitude', 'Longitude', 'bed_thickness', 'sandmass_sum', 'mudmass_sum'],
    #                               header=['Latitude', 'Longitude', 'bed_thickness (m)', 'sandmass_sum (kg m^-2)',
    #                                       'mudmass_sum (kgm^-2)'],
    #                               index_label=['ISO Datetime'],
    #                               date_format='%Y-%m-%dT%H:%M:%SZ',
    #                               )

# Make time-series plots
for site in point_data_veg:
    var2plot = ['bed_thickness', 'sandmass_sum', 'mudmass_sum']
    lat = point_data_veg[site]['Latitude'].unique()
    lon = point_data_veg[site]['Longitude'].unique()
    point_data_veg[site][var2plot].plot(subplots=True,
                                        sharex=True,
                                        title='Site %s (%s, %s)' % (site, lat[0], lon[0]),
                                        )
    image_name = 'data/site_%s.png' % site
    # outfile = writedir+image_name
    # print("Saving image to %s" % outfile)
    plt.savefig(image_name, bbox_inches='tight', dpi=500)