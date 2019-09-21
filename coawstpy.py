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
    times = dict()
    times['typical'] = ['2011-08-01', '2011-08-06']
    times['Irene'] = ['2011-08-27', '2011-08-30']
    times['Lee'] = ['2011-09-07', '2011-09-16']
    times['post-Lee'] = ['2011-10-13', '2011-10-24']
    return times

