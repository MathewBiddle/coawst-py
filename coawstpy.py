import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
#import scipy.integrate as integrate


def stick_plot(time, u, v, **kw):
    '''
    :param time:
    :param u:
    :param v:
    :param kw:
    :return:
    '''

    width = kw.pop('width', 0.001)
    headwidth = kw.pop('headwidth', 0)
    headlength = kw.pop('headlength', 0)
    headaxislength = kw.pop('headaxislength', 0)
    angles = kw.pop('angles', 'uv')
    ax = kw.pop('ax', None)

    if angles != 'uv':
        raise AssertionError("Stickplot angles must be 'uv' so that"
                             "if *U*==*V* the angle of the arrow on"
                             "the plot is 45 degrees CCW from the *x*-axis.")

    time, u, v = map(np.asanyarray, (time, u, v))
    if not ax:
        fig, ax = plt.subplots()

    q = ax.quiver(mdates.date2num(time), [[0] * len(time)], u, v,
                  angles='uv', width=width, headwidth=headwidth,
                  headlength=headlength, headaxislength=headaxislength,
                  **kw)

    ax.axes.get_yaxis().set_ticks([])
    #ax.axes.get_yaxis().set_visible(False)
    ax.xaxis_date()
    return q


def maj_ax(ue, vn):
    '''
    # script to calculate major axis direction in compass direction from
    # east and north components. By LP Sanford, Nov 27, 2006.
    #  translated to python by Mathew Biddle June 26, 2019
    # ue = east velocity
    # vn = north velocity
    # angle = compass direction of major axis of demeaned current vectors
    '''
    #uem = np.mean(ue)
    #vnm = np.mean(vn)
    #uep = ue-uem
    #vnp = vn-vnm
    #n = len(uep)

    cross = np.nansum(np.multiply(ue, vn))
    uesq = np.nansum(np.square(ue))
    vnsq = np.nansum(np.square(vn))

    angle = 0.5*np.arctan(2*cross/(uesq-vnsq))*180/np.pi  # minor axis geometric coords
    angle = 360-angle  # major axis compass direction
    if angle > 360:
        angle = angle-360

    return angle


def rot2xy(ue, vn, projangle):
    '''
    # function [ux,vy]=rot2xy(ue,vn,projangle)
    # script to rotate east and north components into along and cross channel
    # components. By SE Suttles, Nov 26, 2006.
    #  translated to python by Mathew Biddle June 26, 2019
    # ue = east velocity
    # vn = north velocity
    # projangle is the angle on which you wish to project the positive ux
        # current vector, measured in the geographic/compass convention as
        # degrees clockwise from North
    # ux - along channel velocity
    # vy - cross channel velocity (positive to the left of the along channel
        # velocity (right-hand convention with vertical + up)
    '''
    theta = 450-projangle  # convert compass direction to geometric direction
    ux = +ue*np.cos(theta*np.pi/180) + vn*np.sin(theta*np.pi/180)
    vy = -ue*np.sin(theta*np.pi/180) + vn*np.cos(theta*np.pi/180)
    return ux, vy


def nearest_ind(items, value):
    '''
    :param items: list of values
    :param value: value to find closest
    :return: index of closest value in list
    '''
    diff = np.abs([item - value for item in items])
    return diff.argmin(0)


def time_periods():
    '''
    Provides the following time periods for subsetting:

    typical =  '2011-08-01' to '2011-08-06'
    Irene =  '2011-08-27' to '2011-08-30'
    Lee = '2011-09-07' to '2011-09-16'
    post-Lee = '2011-10-13' to '2011-10-24'

    :return:
    A dictionary of time periods for the different events
    Each event key returns a list of start and end dates.

    times[key] = [start,end]
    '''
    import datetime
    times = dict()
    fmt = '%Y-%m-%d'
    times['typical'] = [datetime.datetime.strptime('2011-08-01',fmt), datetime.datetime.strptime('2011-08-06',fmt)]
    times['Irene'] = [datetime.datetime.strptime('2011-08-27',fmt), datetime.datetime.strptime('2011-08-30',fmt)]
    times['Lee'] = [datetime.datetime.strptime('2011-09-07',fmt), datetime.datetime.strptime('2011-09-16',fmt)]
    times['post-Lee'] = [datetime.datetime.strptime('2011-10-13',fmt), datetime.datetime.strptime('2011-10-24',fmt)]
    return times


def get_transect_indexes():
    ###########################
    # Susquehanna River mouth #
    ###########################
    transect=dict()
    transect['T1'] = dict()
    transect['T1']['x'] = np.array(list(range(27, 34)))
    transect['T1']['y'] = np.array([10] * len(transect['T1']['x']))

    ###############################
    # Turkey Point to Sandy Point #
    ###############################
    transect['T2'] = dict()
    transect['T2']['x'] = np.array(list(range(40, 69)))
    transect['T2']['y'] = np.array([58] * len(transect['T2']['x']))
    return transect


def get_sed_flux_data(run,transect_indexes):
    '''
    Computes the sediment flux across transects by computing the flux sum in South and East direction.

    :param run: Select which run to grab the data for.
    Options are:
    'noveg'
    'veg'

    :param transects: a dictionary of transect indexes
    {'T1': {'x': array([ ]),
        'y': array([ ])},

     'T2': {'x': array([ ]),
        'y': array([ ])}
    }
    :return:
    transect - a dictionary of the sediment flux data for each transect
    transect = {'T1': {...}, 'T2': {...}, 'time': [...]}
    '''
    import netCDF4

    import copy

    # bring in the data
    runs_dir = '/Users/mbiddle/Documents/Personal_Documents/Graduate_School/Thesis/COAWST/COAWST_RUNS/COAWST_OUTPUT/'
    if run == 'noveg':
        direct = runs_dir+'Full_20110719T23_20111101_final_noveg'
    elif run == 'veg':
        direct = runs_dir+'Full_20110719T23_20111101_final'
    inputfile = direct+'/upper_ches_avg.nc'
    print('Reading %s %s...' % (inputfile.split("/")[-1],run))
    f = netCDF4.Dataset(inputfile, 'r')

    # Get the data we want
    ocean_time = f.variables['ocean_time'][:]
    lat = f.variables['lat_rho'][:]
    lon = f.variables['lon_rho'][:]
    mask_rho = f.variables['mask_rho'][:]

    Huon = dict()
    # sand
    Huon['sand'] = f.variables['Huon_sand_01'][:] # South

    # mud
    Huon['mud'] = f.variables['Huon_mud_01'][:] # South

    transect = copy.deepcopy(transect_indexes)
    transect['time'] = []
    for sec in ocean_time:
        transect['time'].append(
            netCDF4.num2date(sec, units=f.variables['ocean_time'].units, calendar=f.variables['ocean_time'].calendar))

    sediment_type = ['sand','mud']
    ## Iterate through each transect
    for t in transect: # each transect
        if t == 'time':
            continue
        # get index of interest
        x = transect[t]['x']
        y = transect[t]['y']
        mask_rho[x, y] = 5

        for sed in sediment_type:
            transect[t]['Huon_%s' % sed] = Huon[sed][:, :, x, y] # kg/s
            transect[t]['Huon_%s_units' % sed] = 'kg/s'

            transect[t]['Huon_sum_%s' % sed] = np.sum(transect[t]['Huon_%s' % sed], axis=(1, 2)) # kg/s
            transect[t]['Huon_sum_%s_units' % sed] = 'kg/s'

            #transect[t]['Huon_%s_integrated' % sed] = (integrate.trapz(transect[t]['Huon_sum_%s' % sed]) * 3600)/1000 # integrate per hour (tons)
            #transect[t]['Huon_%s_integrated_units' % sed] = 'tons'

    return transect