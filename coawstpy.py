import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import scipy.integrate as integrate


def stick_plot(time, u, v, **kw):

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


def sediment_flux1(mud_flux_xe,mud_flux_yn,sand_flux_xe,sand_flux_yn,srho_angle,tx1,tx2):
    '''
    Calculates the sediment flux through already subsetted data.
    Flux data is provided as raw values for each cell in the transect
    :param mud_flux_xe: mud flux in the east direction
    :param mud_flux_yn: mud flux in the west direction
    :param sand_flux_xe: sand flux in the east direction
    :param sand_flux_yn: sand flux in the north direction
    :param srho_angle: depth coordinate for angle computation
    :param tx1: start time index for angle computation
    :param tx2: end time index for angle computation
    :return: dictionary of all flux values

    outdict['angle_list']
    outdict['angle']
    outdict['mud_flux_ux_rot']
    outdict['mud_flux_vy_rot']
    outdict['sand_flux_ux_rot']
    outdict['sand_flux_vy_rot']
    outdict['mag_mud']
    outdict['mag_sand']
    outdict['mag_ssc']
    outdict['SSC_ts_sum']
    outdict['cumtrapz_mud']
    outdict['cumtrapz_sand']
    outdict['cumtrapz_ssc']
    outdict['total_mud']
    outdict['total_sand']
    outdict['total_sed']
    '''
    mud_flux_ux_rot = np.empty(mud_flux_xe.shape)
    mud_flux_vy_rot = np.empty(mud_flux_yn.shape)
    sand_flux_ux_rot = np.empty(sand_flux_xe.shape)
    sand_flux_vy_rot = np.empty(sand_flux_yn.shape)

    # compute angle and rotate
    angle = maj_ax(mud_flux_xe[tx1:tx2, srho_angle, :], mud_flux_yn[tx1:tx2, srho_angle, :])
    angle_list = []
    for t in range(0, mud_flux_yn.shape[0]):  # time
        angle_list.append(maj_ax(mud_flux_xe[t, srho_angle, :], mud_flux_yn[t, srho_angle, :]))
        for xy in range(0, mud_flux_yn.shape[2]):  # cell in xy
            for z in range(0, mud_flux_yn.shape[1]): # depth
                # ux - along channel velocity
                # vy - cross channel velocity (positive to the left of the along channel
                mud_flux_ux_rot[t, z, xy], mud_flux_vy_rot[t, z, xy] = rot2xy(
                    mud_flux_xe[t, z, xy], mud_flux_yn[t, z, xy], angle)

                sand_flux_ux_rot[t, z, xy], sand_flux_vy_rot[t, z, xy] = rot2xy(
                    sand_flux_xe[t, z, xy], sand_flux_yn[t, z, xy], angle)

    # magnitude of the rotated flux
    mag_mud = np.sqrt(mud_flux_ux_rot ** 2 + mud_flux_vy_rot ** 2)
    mag_sand = np.sqrt(sand_flux_ux_rot ** 2 + sand_flux_vy_rot ** 2)
    mag_ssc = mag_mud + mag_sand  # total sediment together
    SSC_ts_sum = np.sum(mag_ssc, axis=(1, 2))  # sum over depth and transect
    cumsum_ssc = np.nancumsum(SSC_ts_sum, axis=0)  # cumulative sum over time
    cumtrapz_mud = integrate.cumtrapz(np.sum(mag_mud,axis=(1,2)), axis=0)
    cumtrapz_sand = integrate.cumtrapz(np.sum(mag_sand,axis=(1,2)), axis=0)
    cumtrapz_ssc = integrate.cumtrapz(SSC_ts_sum, axis=0)  # cumulative integral over time
    total_mud = cumtrapz_mud[-1] * 3600
    total_sand = cumtrapz_sand[-1] * 3600
    total_sed = cumtrapz_ssc[-1] * 3600  # sum is integrated for every hour, multiply by 3600 seconds = total weight

    outdict = dict()
    outdict['angle_list'] = angle_list
    outdict['angle'] = angle
    outdict['mud_flux_ux_rot'] = mud_flux_ux_rot
    outdict['mud_flux_vy_rot'] = mud_flux_vy_rot
    outdict['sand_flux_ux_rot'] = sand_flux_ux_rot
    outdict['sand_flux_vy_rot'] = sand_flux_vy_rot
    outdict['mag_mud'] = mag_mud
    outdict['mag_sand'] = mag_sand
    outdict['mag_ssc'] = mag_ssc
    outdict['SSC_ts_sum'] = SSC_ts_sum
    outdict['cumtrapz_mud'] = cumtrapz_mud
    outdict['cumtrapz_sand'] = cumtrapz_sand
    outdict['cumtrapz_ssc'] = cumtrapz_ssc
    outdict['total_mud'] = total_mud
    outdict['total_sand'] = total_sand
    outdict['total_sed'] = total_sed
    return outdict


def sediment_flux2(mud_flux_xe,mud_flux_yn,sand_flux_xe,sand_flux_yn,tx1,tx2):
    ''' #srho_angle
    Calculates the sediment flux through already subsetted data.
    Flux values are provided as a single East and North component for entire transect
    :param mud_flux_xe: mud flux in the east direction
    :param mud_flux_yn: mud flux in the west direction
    :param sand_flux_xe: sand flux in the east direction
    :param sand_flux_yn: sand flux in the north direction
    :param srho_angle: depth coordinate for angle computation
    :param tx1: start time index for angle computation
    :param tx2: end time index for angle computation
    :return: dictionary of all flux values

    outdict['angle_list']
    outdict['angle']
    outdict['mud_flux_ux_rot']
    outdict['mud_flux_vy_rot']
    outdict['sand_flux_ux_rot']
    outdict['sand_flux_vy_rot']
    outdict['mag_mud']
    outdict['mag_sand']
    outdict['mag_ssc']
    outdict['SSC_ts_sum']
    outdict['cumtrapz_mud']
    outdict['cumtrapz_sand']
    outdict['cumtrapz_ssc']
    outdict['total_mud']
    outdict['total_sand']
    outdict['total_sed']
    '''
    mud_flux_ux_rot = np.empty(mud_flux_xe.shape)
    mud_flux_vy_rot = np.empty(mud_flux_yn.shape)
    sand_flux_ux_rot = np.empty(sand_flux_xe.shape)
    sand_flux_vy_rot = np.empty(sand_flux_yn.shape)

    # compute angle and rotate
    angle = maj_ax(mud_flux_xe[tx1:tx2], mud_flux_yn[tx1:tx2])
    angle_list = []
    for t in range(0, mud_flux_yn.shape[0]):  # time
        angle_list.append(maj_ax(mud_flux_xe[t], mud_flux_yn[t]))
        mud_flux_ux_rot[t], mud_flux_vy_rot[t] = rot2xy(mud_flux_xe[t], mud_flux_yn[t], angle)
        sand_flux_ux_rot[t], sand_flux_vy_rot[t] = rot2xy(sand_flux_xe[t], sand_flux_yn[t], angle)

    # magnitude of the rotated flux
    mag_mud = np.sqrt(mud_flux_ux_rot ** 2 + mud_flux_vy_rot ** 2)
    mag_sand = np.sqrt(sand_flux_ux_rot ** 2 + sand_flux_vy_rot ** 2)
    mag_ssc = mag_mud + mag_sand  # total sediment together
    cumtrapz_ssc = integrate.cumtrapz(mag_ssc, axis=0)
    cumtrapz_mud = integrate.cumtrapz(mag_mud, axis=0)
    cumtrapz_sand = integrate.cumtrapz(mag_sand, axis=0)
    total_mud = cumtrapz_mud[-1] * 3600
    total_sand = cumtrapz_sand[-1] * 3600
    total_sed = cumtrapz_ssc[-1] * 3600  # sum is integrated for every hour, multiply by 3600 seconds = total weight

    outdict = dict()
    outdict['angle_list'] = angle_list
    outdict['angle'] = angle
    outdict['mud_flux_ux_rot'] = mud_flux_ux_rot
    outdict['mud_flux_vy_rot'] = mud_flux_vy_rot
    outdict['sand_flux_ux_rot'] = sand_flux_ux_rot
    outdict['sand_flux_vy_rot'] = sand_flux_vy_rot
    outdict['mag_mud'] = mag_mud
    outdict['mag_sand'] = mag_sand
    outdict['mag_ssc'] = mag_ssc
    outdict['cumtrapz_mud'] = cumtrapz_mud
    outdict['cumtrapz_sand'] = cumtrapz_sand
    outdict['cumtrapz_ssc'] = cumtrapz_ssc
    outdict['total_mud'] = total_mud
    outdict['total_sand'] = total_sand
    outdict['total_sed'] = total_sed
    return outdict


def sediment_flux3(xe,yn,tx1,tx2):
    ''' #srho_angle
    Calculates the sediment flux through already subsetted data.
    Flux values are provided as a single East and North component for entire transect
    :param xe: flux in the east direction
    :param yn: flux in the west direction
    :param tx1: start time index for angle computation
    :param tx2: end time index for angle computation
    :return: dictionary of all flux values

    outdict['angle_list']
    outdict['angle']
    outdict['ux_rot']
    outdict['vy_rot']
    outdict['mag']
    outdict['cumtrapz']
    outdict['total']
    '''

    ux_rot = np.empty(xe.shape)
    vy_rot = np.empty(yn.shape)

    # compute angle and rotate
    angle = maj_ax(xe[tx1:tx2], yn[tx1:tx2])
    angle_list = []
    #for t in range(0, yn.shape[0]):  # time
    #    angle_list.append(maj_ax(xe[t], yn[t]))
    #    ux_rot[t], vy_rot[t] = rot2xy(xe[t], yn[t], angle)
    ux_rot, vy_rot = rot2xy(xe, yn, angle)

    # magnitude of the rotated flux
    mag = np.sqrt(ux_rot ** 2 + vy_rot ** 2)
    cumtrapz = integrate.cumtrapz(mag, axis=0)

    outdict = dict()
    outdict['xe'] = xe
    outdict['yn'] = yn
    #outdict['angle_list'] = angle_list
    outdict['angle'] = angle
    outdict['ux_rot'] = ux_rot
    outdict['vy_rot'] = vy_rot
    outdict['mag'] = mag
    outdict['cumtrapz'] = cumtrapz
    return outdict