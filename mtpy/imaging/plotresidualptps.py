# -*- coding: utf-8 -*-
"""
PlotResidualPhaseTensorPseudoSection
=======================

    *plots the residual phase tensor for two different sets of measurments


Created on Wed Oct 16 14:56:04 2013

@author: jpeacock-pr
"""
#==============================================================================
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.colorbar as mcb
import mtpy.imaging.mtcolors as mtcl
import mtpy.imaging.mtplottools as mtpl
import mtpy.analysis.pt as mtpt
import scipy.signal as sps

#==============================================================================


class PlotResidualPTps(mtpl.MTEllipse):
    """
    This will plot residual phase tensors in a pseudo section.  The data is
    read in and stored in 2 ways, one as a list ResidualPhaseTensor object for
    each matching station and the other in a structured array with all the
    important information.  The structured array is the one that is used for
    plotting.  It is computed each time plot() is called so if it is
    manipulated it is reset.  The array is sorted by relative offset, so no
    special order of input is needed for the file names.  However, the
    station names should be verbatim between surveys, otherwise it will not
    work.

    The residual phase tensor is calculated as I-(Phi_2)^-1 (Phi_1)

    The default coloring is by the geometric mean as sqrt(Phi_min*Phi_max),
    which defines the percent change between measurements.

    There are a lot of parameters to change how the plot looks, have a look
    below if you figure looks a little funny.  The most useful will be
    xstretch, ystretch and ellipse_size

    The ellipses are normalized by the largest Phi_max of the survey.


    Arguments:
    --------------

        **fn_list1** : list of strings
                        full paths to .edi files for survey 1

        **fn_list2** : list of strings
                        full paths to .edi files for survey 2

        **rot90** : [ True | False ] True to rotate residual phase tensor
                    by 90 degrees. False to leave as is.

        **ellipse_dict** : dictionary
                          dictionary of parameters for the phase tensor
                          ellipses with keys:
                              * 'size' -> size of ellipse in points
                                         *default* is 2

                              * 'colorby' : [ 'phimin' | 'phimax' | 'skew' |
                                              'skew_seg' | 'phidet' |
                                              'ellipticity' | 'geometric_mean']

                                        - 'phimin' -> colors by minimum phase
                                        - 'phimax' -> colors by maximum phase
                                        - 'skew' -> colors by beta (skew)
                                        - 'skew_seg' -> colors by beta in
                                                       discrete segments
                                                       defined by the range
                                        - 'phidet' -> colors by determinant of
                                                     the phase tensor
                                        - 'ellipticity' -> colors by ellipticity
                                        - 'geometric_mean' -> sqrt(phimin*phimax)
                                        *default* is 'geometric_mean'

                               * 'range' : tuple (min, max, step)
                                     Need to input at least the min and max
                                     and if using 'skew_seg' to plot
                                     discrete values input step as well
                                     *default* depends on 'colorby'

                          * 'cmap' : [ 'mt_yl2rd' | 'mt_bl2yl2rd' |
                                      'mt_wh2bl' | 'mt_rd2bl' |
                                      'mt_bl2wh2rd' | 'mt_seg_bl2wh2rd' |
                                      'mt_rd2gr2bl' | 'mt_wh2or ]

                                   - 'mt_yl2rd' -> yellow to red *default*
                                   - 'mt_bl2yl2rd' -> blue to yellow to red
                                   - 'mt_wh2bl' -> white to blue
                                   - 'mt_rd2bl' -> red to blue
                                   - 'mt_bl2wh2rd' -> blue to white to red
                                   - 'mt_bl2gr2rd' -> blue to green to red
                                   - 'mt_rd2gr2bl' -> red to green to blue
                                   - 'mt_seg_bl2wh2rd' -> discrete blue to
                                                         white to red
                                   - 'mt_wh2or' -> white to orange


        **ellipse_scale** : float
                            value to which all ellipses are normalized to. So
                            if you think a the maximum of change in the
                            phase tensor is about 10 degrees, set this value
                            to 10, then any phimax that is 10 degrees will make
                            an ellipse with unit length.


        **med_filt_kernel** : tuple(station, period)
                              kernel size for the 2D median filter.
                              the first is the number of stations to smooth
                              over. The first number is the number of periods
                              to smooth over. Both should be odd.
                              *default* is None

        **xstretch** : float
                        is a factor that scales the distance from one
                        station to the next to make the plot readable.
                        *Default* is 1000
        **ystretch** : float
                        is a factor that scales the distance from one
                        period to the next to make the plot readable.
                        *Default* is 25

        **linedir** : [ 'ns' | 'ew' ]
                      predominant direction of profile line
                      * 'ns' -> North-South Line
                      * 'ew' -> East-West line
                      *Default* is 'ns'

        **station_id** : tuple or list
                        start and stop of station name indicies.
                        ex: for MT01dr station_id=(0,4) will be MT01

        **rotz** : float or np.ndarray
                   angle in degrees to rotate the data clockwise positive.
                   Can be an array of angle to individually rotate stations or
                   periods or both.
                       - If rotating each station by a constant
                         angle the array needs to have a shape of
                         (# of stations)
                        - If rotating by period needs to have shape
                           # of periods
                        - If rotating both individually shape=(ns, nf)
                  *Default* is 0

        **title** : string
                    figure title

        **dpi** : int
                  dots per inch of the resolution. *default* is 300


        **fig_num** : int
                     figure number.  *Default* is 1


        **tscale** : [ 'period' | 'frequency' ]

                     * 'period'    -> plot vertical scale in period

                     * 'frequency' -> plot vertical scale in frequency

        **cb_dict** : dictionary to control the color bar

                      * 'orientation' : [ 'vertical' | 'horizontal' ]
                                       orientation of the color bar
                                       *default* is vertical

                      * 'position' : tuple (x,y,dx,dy)
                                    - x -> lateral position of left hand corner
                                          of the color bar in figure between
                                          [0,1], 0 is left side

                                    - y -> vertical position of the bottom of
                                          the color bar in figure between
                                          [0,1], 0 is bottom side.

                                    - dx -> width of the color bar [0,1]

                                    - dy -> height of the color bar [0,1]
        **font_size** : float
                        size of the font that labels the plot, 2 will be added
                        to this number for the axis labels.

        **plot_yn** : [ 'y' | 'n' ]
                      * 'y' to plot on creating an instance

                      * 'n' to not plot on creating an instance

        **xlimits** : tuple(xmin, xmax)
                   min and max along the x-axis in relative distance of degrees
                   and multiplied by xstretch

        **ylimits** : tuple(ymin, ymax)
                   min and max period to plot, note that the scaling will be
                   done in the code.  So if you want to plot from (.1s, 100s)
                   input ylim=(.1,100)


    ==================== ======================================================
      Attributes          Description
    ==================== ======================================================
     ax                   matplotlib.axes instance for the main plot
     ax2                  matplotlib.axes instance for the color bar
     cb                   matplotlib.colors.ColorBar instance for color bar
     cb_orientation       color bar orientation ('vertical' | 'horizontal')
     cb_position          color bar position (x, y, dx, dy)
     ellipse_cmap         ellipse color map, see above for options
     ellipse_colorby      parameter to color ellipse by
     ellipse_range        (min, max, step) values to color ellipses
     ellipse_scale        value to which all ellipses are normalized to
     ellipse_size         scaling factor to make ellipses visible
     fig                  matplotlib.figure instance for the figure
     fig_dpi              dots-per-inch resolution
     fig_num              number of figure being plotted
     fig_size             size of figure in inches
     fn_list1             list of .edi file names for survey 1
     fn_list2             list of .edi file names for survey 2
     font_size            font size of axes tick label, axes labels will be
                          font_size + 2
     freq_list            list of frequencies from all .edi files
     linedir              prominent direction of profile being plotted
     med_filt_kernel      (station, frequency) kernel to apply median smoothing
                          to the data.
     mt_list1             list of mtplot.MTplot instances containing all
                          important information for each station in survey 1
     mt_list2             list of mtplot.MTplot instances containing all
                          important information for each station in survey 2
     offset_list          array of relative offsets of each station
     plot_title           title of the plot
     plot_yn              plot the pseudo section on instance creation
     residual_pt_list     list ofmtpy.pt.ResidualPhaseTensor objects
     rot90                rotates the residual phase tensors by 90 degrees
                          if set to True
     rpt_array            structured array with all the important information.
                          This is the important array from which plotting
                          occurs.
     station_font_dict    font dictionary for station labels
     station_id           index [min, max] to reaad station name
     station_list         list of stations plotted
     station_pad          padding between axis and station label
     subplot_bottom       spacing between plot and bottom of figure window
     subplot_hspace       vertical spacing between subplots
     subplot_left         spacing between plot and left of figure window
     subplot_right        spacing between plot and right of figure window
     subplot_top          spacing between plot and top of figure window
     subplot_wspace       horizontal spacing between subplots
     tscale               temporal scale of y-axis ('frequency' | 'period')
     xlimits              limits on x-axis (xmin, xmax)
     xstretch             scaling factor to stretch x offsets
     xstep                interval of station labels to plot
                          *default* is 1 to label all stations
     ylimits              limits on y-axis (ymin, ymax)
     ystep                interval of period labels to plot
                          *default* is 1 to label all powers of 10
     ystretch             scaling factor to strech axes in y direction
    ==================== ======================================================


    ======================= ===================================================
    Methods                 Description
    ======================= ===================================================
     plot                   plots the pseudo section
     redraw_plot            on call redraws the plot from scratch
     save_figure            saves figure to a file of given format
     update_plot            updates the plot while still active
     _apply_median_filter   apply a 2D median filter to the data
     _compute_residual_pt   compute residual pt and fill rpt_array and
                            fill residaul_pt_list
     _get_freq_list         get a list of all possible frequencies from .edi's
     _get_offsets           get a list of offsets of the station locations
                            and fill rpt_array['offset']
     _read_ellipse_dict     read ellipse dictionary and return a ellipse object
    ======================= ===================================================

    To get a list of .edi files that you want to plot -->
    :Example: ::

        >>> import mtpy.imaging.mtplot as mtplot
        >>> import os
        >>> edipath1 = r"/home/EDIfiles1"
        >>> edilist1 = [os.path.join(edipath1,edi) for edi in os.listdir(edipath1)
        >>> ...       if edi.find('.edi')>0]
        >>> edipath2 = r"/home/EDIfiles2"
        >>> edilist2 = [os.path.join(edipath2,edi) for edi in os.listdir(edipath2)
        >>> ...       if edi.find('.edi')>0]
        >>> # color by phimin with a range of 0-5 deg

    * If you want to plot geometric mean colored from white to orange in a
      range of 0 to 10 percent you can do it one of two ways-->

    1)
    :Example: ::

        >>> edict = {'range':(0,10), 'cmap':'mt_wh2or', \
                     'colorby':'geometric_mean', 'size':10}
        >>> pt1 = mtplot.residual_pt_ps(edilist1, edilst2, ellipse_dict=edict)

    2)
    :Example: ::

        >>> pt1 = mtplot.residual_pt_ps(edilist1, edilst2, ellipse_dict=edict,\
                                        plot_yn='n')
        >>> pt1.ellipse_colorby = 'geometric_mean'
        >>> pt1.ellipse_cmap = 'mt_wh2or'
        >>> pt1.ellipse_range = (0, 10)
        >>> pt1.ellipse_size = 10
        >>> pt1.plot()


    * If you want to save the plot as a pdf with a generic name -->
    :Example: ::
        >>> pt1.save_figure(r"/home/PTFigures", file_format='pdf', dpi=300)
        File saved to '/home/PTFigures/PTPseudoSection.pdf'

    """

    def __init__(self, fn_list1, fn_list2, **kwargs):
        assert len(fn_list1) == len(fn_list2)
        self.fn_list1 = fn_list1
        self.fn_list2 = fn_list2

        self.mt_list1 = mtpl.get_mtlist(fn_list=fn_list1)
        self.mt_list2 = mtpl.get_mtlist(fn_list=fn_list2)

        self.residual_pt_list = []
        self.freq_list = None
        self.rpt_array = None
        self.med_filt_kernel = kwargs.pop('med_filt_kernel', None)
        self.rot90 = kwargs.pop('rot90', True)

        #--> set the ellipse properties
        self._ellipse_dict = kwargs.pop('ellipse_dict',
                                        {'cmap': 'mt_yl2rd',
                                         'range': (0, 10),
                                         'colorby': 'geometric_mean'})
        self._read_ellipse_dict()
        self.ellipse_scale = kwargs.pop('ellipse_scale', 10)

        #--> set colorbar properties---------------------------------
        # set orientation to horizontal
        cb_dict = kwargs.pop('cb_dict', {})
        try:
            self.cb_orientation = cb_dict['orientation']
        except KeyError:
            self.cb_orientation = 'vertical'

        # set the position to middle outside the plot
        try:
            self.cb_position = cb_dict['position']
        except KeyError:
            self.cb_position = None

        #--> set plot properties ------------------------------
        self.fig_num = kwargs.pop('fig_num', 'residual_PT')
        self.plot_title = kwargs.pop('plot_title', None)
        self.fig_dpi = kwargs.pop('fig_dpi', 300)
        self.tscale = kwargs.pop('tscale', 'period')
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.linedir = kwargs.pop('linedir', 'ew')
        self.font_size = kwargs.pop('font_size', 7)
        self.station_id = kwargs.pop('station_id', [0, 4])
        self.ystep = kwargs.pop('ystep', 4)
        self.xstep = kwargs.pop('xstep', 1)
        self.xlimits = kwargs.pop('xlimits', None)
        self.ylimits = kwargs.pop('ylimits', None)

        #--> set spacing of plot
        self.subplot_wspace = .1
        self.subplot_hspace = .15
        self.subplot_right = .98
        self.subplot_left = .085
        self.subplot_top = .93
        self.subplot_bottom = .1

        # set the stretching in each direction
        self.xstretch = kwargs.pop('xstretch', 1000)
        self.ystretch = kwargs.pop('ystretch', 25)

        # if rotation angle is an int or float make an array the length of
        # mt_list for plotting purposes

        self._rot_z = kwargs.pop('rot_z', 0)
        if isinstance(self._rot_z, float) or isinstance(self._rot_z, int):
            self._rot_z = np.array([self._rot_z] * len(self.mt_list1))

        # if the rotation angle is an array for rotation of different
        # freq than repeat that rotation array to the len(mt_list)
        elif isinstance(self._rot_z, np.ndarray):
            if self._rot_z.shape[0] != len(self.mt_list1):
                self._rot_z = np.repeat(self._rot_z, len(self.mt_list1))

        else:
            pass

        #--> set station name properties
        station_dict = kwargs.pop('station_dict', {})
        self.station_id = station_dict.pop('id', self.station_id)
        self.station_pad = station_dict.pop('pad', .0005)
        self.station_font_dict = station_dict.pop('font_dict',
                                                  {'size': self.font_size,
                                                   'weight': 'bold'})

        #--> plot if desired ------------------------
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()

    #---need to rotate data on setting rotz
    def _set_rot_z(self, rot_z):
        """
        need to rotate data when setting z
        """

        # if rotation angle is an int or float make an array the length of
        # mt_list for plotting purposes
        if isinstance(rot_z, float) or isinstance(rot_z, int):
            self._rot_z = np.array([rot_z] * len(self.mt_list))

        # if the rotation angle is an array for rotation of different
        # freq than repeat that rotation array to the len(mt_list)
        elif isinstance(rot_z, np.ndarray):
            if rot_z.shape[0] != len(self.mt_list):
                self._rot_z = np.repeat(rot_z, len(self.mt_list))

        else:
            pass

        for ii, rpt in enumerate(self.residual_pt_list):
            rpt.rot_z = self._rot_z[ii]

    def _get_rot_z(self):
        return self._rot_z

    rot_z = property(fget=_get_rot_z, fset=_set_rot_z,
                     doc="""rotation angle(s)""")
    #-------------------------------------------------------------------

    def _get_freq_list(self):
        """
        get all possible periods to plot

        """

        freq_list = []
        for mt1 in self.mt_list1:
            freq_list.extend(mt1.freq)
        for mt2 in self.mt_list2:
            freq_list.extend(mt2.freq)

        self.freq_list = np.array(sorted(set(freq_list), reverse=True))

    #------------------------------------------------------------------
    def _compute_residual_pt(self):
        """
        compute residual phase tensor so the result is something useful to
        plot
        """
        log_path = os.path.dirname(os.path.dirname(self.fn_list1[0]))
        log_fn = os.path.join(log_path, 'Residual_PT.log')
        logfid = file(log_fn, 'w')

        self._get_freq_list()
        freq_dict = dict([(np.round(key, 5), value)
                          for value, key in enumerate(self.freq_list)])

        num_freq = self.freq_list.shape[0]
        num_station = len(self.mt_list1)

        # make a structured array to put stuff into for easier manipulation
        self.rpt_array = np.zeros(num_station,
                                  dtype=[('station', '|S10'),
                                         ('freq', (np.float, num_freq)),
                                         ('lat', np.float),
                                         ('lon', np.float),
                                         ('elev', np.float),
                                         ('offset', np.float),
                                         ('phimin', (np.float, num_freq)),
                                         ('phimax', (np.float, num_freq)),
                                         ('skew', (np.float, num_freq)),
                                         ('azimuth', (np.float, num_freq)),
                                         ('geometric_mean', (np.float, num_freq))])

        self.residual_pt_list = []
        for mm, mt1 in enumerate(self.mt_list1):
            station_find = False
            fdict1 = dict([(np.round(ff, 5), ii)
                           for ii, ff in enumerate(mt1.freq)])
            for mt2 in self.mt_list2:
                if mt2.station == mt1.station:
                    logfid.write('{0}{1}{0}\n'.format('=' * 30, mt1.station))
                    fdict2 = dict([(np.round(ff, 5), ii)
                                   for ii, ff in enumerate(mt2.freq)])

                    # need to make sure only matched frequencies are compared
                    index_1 = []
                    index_2 = []
                    for key1 in sorted(fdict1.keys()):
                        try:
                            index_2.append(fdict2[key1])
                            index_1.append(fdict1[key1])
                            logfid.write('-' * 20 + '\n')
                            logfid.write('found {0} in both\n'.format(key1))
                            logfid.write('Index 1={0}, Index 2={1}\n'.format(
                                         fdict1[key1], fdict2[key1]))
                        except KeyError:
                            'Did not find {0:.4e} Hz in {1}'.format(key1,
                                                                    mt2.fn)

                    # need to sort the index list, otherwise weird things
                    # happen
                    index_1.sort()
                    index_2.sort()

                    # create new Z objects that have similar frequencies
                    new_z1 = mtpl.mtz.Z(z_array=mt1.z[index_1],
                                        z_err_array=mt1.z_err[index_1],
                                        freq=mt1.freq[index_1])
                    new_z2 = mtpl.mtz.Z(z_array=mt2.z[index_2],
                                        z_err_array=mt2.z_err[index_2],
                                        freq=mt2.freq[index_2])

                    # make new phase tensor objects
                    pt1 = mtpt.PhaseTensor(z_object=new_z1)
                    pt2 = mtpt.PhaseTensor(z_object=new_z2)

                    # compute residual phase tensor
                    rpt = mtpt.ResidualPhaseTensor(pt1, pt2)
                    rpt.compute_residual_pt(pt1, pt2)

                    # add some attributes to residual phase tensor object
                    rpt.station = mt1.station
                    rpt.lat = mt1.lat
                    rpt.lon = mt1.lon

                    # append to list for manipulating later
                    self.residual_pt_list.append(rpt)

                    # be sure to tell the program you found the station
                    station_find = True

                    # put stuff into an array because we cannot set values of
                    # rpt, need this for filtering.
                    st_1, st_2 = self.station_id
                    self.rpt_array[mm]['station'] = mt1.station[st_1:st_2]
                    self.rpt_array[mm]['lat'] = mt1.lat
                    self.rpt_array[mm]['lon'] = mt1.lon
                    self.rpt_array[mm]['elev'] = mt1.elev
                    self.rpt_array[mm]['freq'][:] = self.freq_list

                    rpt_fdict = dict([(np.round(key, 5), value)
                                      for value, key in enumerate(rpt.freq)])

                    for f_index, freq in enumerate(rpt.freq):
                        aa = freq_dict[np.round(freq, 5)]
                        try:
                            rr = rpt_fdict[np.round(freq, 5)]
                            try:
                                self.rpt_array[mm]['phimin'][aa] = \
                                    rpt.residual_pt.phimin[0][rr]
                                self.rpt_array[mm]['phimax'][aa] = \
                                    rpt.residual_pt.phimax[0][rr]
                                self.rpt_array[mm]['skew'][aa] = \
                                    rpt.residual_pt.beta[0][rr]
                                self.rpt_array[mm]['azimuth'][aa] = \
                                    rpt.residual_pt.azimuth[0][rr]
                                self.rpt_array[mm]['geometric_mean'][aa] = \
                                    np.sqrt(abs(rpt.residual_pt.phimin[0][rr] *
                                                rpt.residual_pt.phimax[0][rr]))
                                logfid.write('Freq={0:.5f} '.format(freq))
                                logfid.write('Freq_list_index={0} '.format(
                                    np.where(self.freq_list == freq)[0][0]))
                                logfid.write('rpt_array_index={0} '.format(aa))
                                logfid.write('rpt_dict={0} '.format(rr))
                                logfid.write('rpt.freq_index={0} '.format(
                                    np.where(rpt.freq == freq)[0][0]))
                                logfid.write('Phi_max={0:2f} '.format(
                                    rpt.residual_pt.phimax[0][rr]))
                                logfid.write('Phi_min={0:2f} '.format(
                                    rpt.residual_pt.phimin[0][rr]))
                                logfid.write('Skew={0:2f} '.format(
                                    rpt.residual_pt.beta[0][rr]))
                                logfid.write('Azimuth={0:2f}\n'.format(
                                    rpt.residual_pt.azimuth[0][rr]))

                            except IndexError:
                                print('-' * 50)
                                print(mt1.station)
                                print('freq_index for 1:  {0}'.format(f_index))
                                print('freq looking for:  {0}'.format(freq))
                                print('index in big    :  {0}'.format(aa))
                                print('index in 1      :  {0} '.format(rr))
                                print('len_1 = {0}, len_2 = {1}'.format(
                                    len(mt2.freq), len(mt1.freq)))
                                print('len rpt_freq = {0}'.format(len(rpt.freq)))
                        except KeyError:
                            print('Station {0} does not have {1:.5f}Hz'.format(
                                mt1.station, freq))

                    break
                else:
                    pass
            if station_find == False:
                print('Did not find {0} from list 1 in list 2'.format(
                    mt1.station))

        # from the data get the relative offsets and sort the data by them
        self._get_offsets()

        logfid.close()
    #-------------------------------------------------------------------

    def _apply_median_filter(self, kernel=(3, 3)):
        """
        apply a median filter to the data to remove extreme outliers

        kernel is (station, frequency)

        """

        filt_phimin_arr = sps.medfilt2d(self.rpt_array['phimin'],
                                        kernel_size=kernel)
        filt_phimax_arr = sps.medfilt2d(self.rpt_array['phimax'],
                                        kernel_size=kernel)
        filt_skew_arr = sps.medfilt2d(self.rpt_array['skew'],
                                      kernel_size=kernel)
        filt_azimuth_arr = sps.medfilt2d(self.rpt_array['azimuth'],
                                         kernel_size=kernel)

        self.rpt_array['phimin'] = filt_phimin_arr
        self.rpt_array['phimax'] = filt_phimax_arr
        self.rpt_array['skew'] = filt_skew_arr
        self.rpt_array['azimuth'] = filt_azimuth_arr
        self.rpt_array['geometric_mean'] = np.sqrt(abs(filt_phimin_arr *
                                                       filt_phimax_arr))

        print('Applying Median Filter with kernel {0}'.format(kernel))

    #-------------------------------------------------------------------
    def _get_offsets(self):
        """
        get relative offsets of stations
        """

        for ii, r_arr in enumerate(self.rpt_array):
            # set the an arbitrary origin to compare distance to all other
                # stations.
            if ii == 0:
                east0 = r_arr['lon']
                north0 = r_arr['lat']
                offset = 0.0
                r_arr['offset'] = offset
            else:
                east = r_arr['lon']
                north = r_arr['lat']
                if self.linedir == 'ew':
                    if east0 < east:
                        offset = np.sqrt((east0 - east)**2 +
                                         (north0 - north)**2)
                    elif east0 > east:
                        offset = -1 * \
                            np.sqrt((east0 - east)**2 + (north0 - north)**2)
                    else:
                        offset = 0
                    r_arr['offset'] = offset
                elif self.linedir == 'ns':
                    if north0 < north:
                        offset = np.sqrt((east0 - east)**2 +
                                         (north0 - north)**2)
                    elif north0 > north:
                        offset = -1 * \
                            np.sqrt((east0 - east)**2 + (north0 - north)**2)
                    else:
                        offset = 0
                    r_arr['offset'] = offset

        # be sure to order the structured array by offset, this will make
        # sure that the median filter is spatially correct
        self.rpt_array.sort(order='offset')

    #--------------------------------------------------------------------------
    def plot(self):
        """
        plot residual phase tensor
        """

        # get residual phase tensor for plotting
        self._compute_residual_pt()

        # filter data if desired
        if self.med_filt_kernel is not None:
            self._apply_median_filter(kernel=self.med_filt_kernel)

        # set position properties for the plot
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        plt.rcParams['figure.subplot.wspace'] = self.subplot_wspace
        plt.rcParams['figure.subplot.hspace'] = self.subplot_hspace

        # make figure instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)

        # create axis instance, be sure to make aspect equal or ellipses will
        # look funny.
        self.ax = self.fig.add_subplot(1, 1, 1, aspect='equal')

        # set local parameters with shorter names
        es = self.ellipse_size
        ck = self.ellipse_colorby
        cmap = self.ellipse_cmap
        ckmin = float(self.ellipse_range[0])
        ckmax = float(self.ellipse_range[1])

        try:
            ckstep = float(self.ellipse_range[2])
        except IndexError:
            ckstep = 3
        # set the number of segments in case a segmented map is desired
        nseg = float((ckmax - ckmin) / (2 * ckstep))

        if cmap == 'mt_seg_bl2wh2rd':
            bounds = np.arange(ckmin, ckmax + ckstep, ckstep)

        # get largest ellipse
        emax = self.ellipse_scale
        #emax = self.rpt_array['phimax'].max()

        # plot phase tensor ellipses
        for ii, rpt in enumerate(self.rpt_array):

            phimax = rpt['phimax']
            phimin = rpt['phimin']
            azimuth = rpt['azimuth']

            # get the properties to color the ellipses by
            try:
                color_array = rpt[self.ellipse_colorby]
            except ValueError:
                raise NameError('{0} is not supported'.format(
                    self.ellipse_colorby))

            for jj, ff in enumerate(self.freq_list):
                if phimin[jj] == 0.0 or phimax[jj] == 0.0:
                    pass
                else:
                    # make sure the ellipses will be visable
                    eheight = phimin[jj] / emax * es
                    ewidth = phimax[jj] / emax * es

                    # create an ellipse scaled by phimin and phimax and orient
                    # the ellipse so that north is up and east is right
                    # need to add 90 to do so instead of subtracting
                    if self.rot90 == True:
                        ellipd = patches.Ellipse((rpt['offset'] * self.xstretch,
                                                  np.log10(ff) * self.ystretch),
                                                 width=ewidth,
                                                 height=eheight,
                                                 angle=azimuth[jj] - 90)
                    else:
                        ellipd = patches.Ellipse((rpt['offset'] * self.xstretch,
                                                  np.log10(ff) * self.ystretch),
                                                 width=ewidth,
                                                 height=eheight,
                                                 angle=azimuth[jj])

                    # get ellipse color
                    if cmap.find('seg') > 0:
                        ellipd.set_facecolor(mtcl.get_plot_color(color_array[jj],
                                                                 self.ellipse_colorby,
                                                                 cmap,
                                                                 ckmin,
                                                                 ckmax,
                                                                 bounds=bounds))
                    else:
                        ellipd.set_facecolor(mtcl.get_plot_color(color_array[jj],
                                                                 self.ellipse_colorby,
                                                                 cmap,
                                                                 ckmin,
                                                                 ckmax))

                    # == =add the ellipse to the plot == ========
                    self.ax.add_artist(ellipd)

        #--> Set plot parameters
        # need to sort the offsets and station labels so they plot correctly
        sdtype = [('offset', np.float), ('station', '|S10')]
        slist = np.array([(oo, ss) for oo, ss in zip(self.rpt_array['offset'],
                                                     self.rpt_array['station'])], dtype=sdtype)

        offset_sort = np.sort(slist, order='offset')

        self.offset_list = offset_sort['offset']
        self.station_list = offset_sort['station']

        # min and max frequency of the plot
        pmin = int(np.floor(np.log10(self.freq_list.min())))
        pmax = int(np.ceil(np.log10(self.freq_list.max())))

        # set y-axis major ticks to be on each power of 10
        self.ax.yaxis.set_ticks(np.arange(pmin * self.ystretch,
                                          (pmax + 1) * self.ystretch,
                                          self.ystretch))
        # set y-ticklabels to coincide with the desired label
        if self.tscale == 'period':
            # make tick labels that will represent period
            yticklabels = [mtpl.labeldict[-ii]
                           for ii in range(pmin, pmax + 1, 1)]
            self.ax.set_ylabel('Period (s)',
                               fontsize=self.font_size + 2,
                               fontweight='bold')

        elif self.tscale == 'frequency':
            yticklabels = [mtpl.labeldict[ii]
                           for ii in range(pmin, pmax + 1, 1)]
            self.ax.set_ylabel('Frequency (Hz)',
                               fontsize=self.font_size + 2,
                               fontweight='bold')
        #--> set y-limits
        if self.ylimits is None:
            self.ax.set_ylim(pmin * self.ystretch, pmax * self.ystretch)
        else:
            pmin = np.log10(self.ylimits[0]) * self.ystretch
            pmax = np.log10(self.ylimits[1]) * self.ystretch
            self.ax.set_ylim(pmin, pmax)

        #--> set y-axis tick labels
        self.ax.set_yticklabels(yticklabels)

        #--> set x-axis label
        self.ax.set_xlabel('Station',
                           fontsize=self.font_size + 2,
                           fontweight='bold')

        # set x-axis ticks
        self.ax.set_xticks(self.offset_list * self.xstretch)

        # set x-axis tick labels as station names
        xticklabels = self.station_list
        if self.xstep != 1:
            xticklabels = np.zeros(len(self.station_list),
                                   dtype=self.station_list.dtype)
            for xx in range(0, len(self.station_list), self.xstep):
                xticklabels[xx] = self.station_list[xx]
        self.ax.set_xticklabels(xticklabels)

        #--> set x-limits
        if self.xlimits is None:
            self.ax.set_xlim(self.offset_list.min() * self.xstretch - es * 2,
                             self.offset_list.max() * self.xstretch + es * 2)
        else:
            self.ax.set_xlim(self.xlimits)

        #--> set title of the plot
        if self.plot_title is None:
            pass
        else:
            self.ax.set_title(self.plot_title, fontsize=self.font_size + 2)

        # put a grid on the plot
        self.ax.grid(alpha=.25, which='both', color=(.25, .25, .25))

        #==> make a colorbar with appropriate colors
        if self.cb_position is None:
            self.ax2, kw = mcb.make_axes(self.ax,
                                         orientation=self.cb_orientation,
                                         shrink=.35)
        else:
            self.ax2 = self.fig.add_axes(self.cb_position)

        if cmap == 'mt_seg_bl2wh2rd':
            # make a color list
            self.clist = [(cc, cc, 1)
                          for cc in np.arange(0, 1 + 1. / (nseg), 1. / (nseg))] +\
                [(1, cc, cc)
                 for cc in np.arange(1, -1. / (nseg), -1. / (nseg))]

            # make segmented colormap
            mt_seg_bl2wh2rd = colors.ListedColormap(self.clist)

            # make bounds so that the middle is white
            bounds = np.arange(ckmin - ckstep, ckmax + 2 * ckstep, ckstep)

            # normalize the colors
            norms = colors.BoundaryNorm(bounds, mt_seg_bl2wh2rd.N)

            # make the colorbar
            self.cb = mcb.ColorbarBase(self.ax2,
                                       cmap=mt_seg_bl2wh2rd,
                                       norm=norms,
                                       orientation=self.cb_orientation,
                                       ticks=bounds[1:-1])
        else:
            self.cb = mcb.ColorbarBase(self.ax2,
                                       cmap=mtcl.cmapdict[cmap],
                                       norm=colors.Normalize(vmin=ckmin,
                                                             vmax=ckmax),
                                       orientation=self.cb_orientation)

        # label the color bar accordingly
        self.cb.set_label(mtpl.ckdict[ck],
                          fontdict={'size': self.font_size, 'weight': 'bold'})

        # place the label in the correct location
        if self.cb_orientation == 'horizontal':
            self.cb.ax.xaxis.set_label_position('top')
            self.cb.ax.xaxis.set_label_coords(.5, 1.3)

        elif self.cb_orientation == 'vertical':
            self.cb.ax.yaxis.set_label_position('right')
            self.cb.ax.yaxis.set_label_coords(1.5, .5)
            self.cb.ax.yaxis.tick_left()
            self.cb.ax.tick_params(axis='y', direction='in')

        #--> add reference ellipse
        ref_ellip = patches.Ellipse((0, .0),
                                    width=es,
                                    height=es,
                                    angle=0)
        ref_ellip.set_facecolor((0, 0, 0))
        ref_ax_loc = list(self.ax2.get_position().bounds)
        ref_ax_loc[0] *= .95
        ref_ax_loc[1] -= .17
        ref_ax_loc[2] = .1
        ref_ax_loc[3] = .1
        self.ref_ax = self.fig.add_axes(ref_ax_loc, aspect='equal')
        self.ref_ax.add_artist(ref_ellip)
        self.ref_ax.set_xlim(-es / 2. * 1.05, es / 2. * 1.05)
        self.ref_ax.set_ylim(-es / 2. * 1.05, es / 2. * 1.05)
        plt.setp(self.ref_ax.xaxis.get_ticklabels(), visible=False)
        plt.setp(self.ref_ax.yaxis.get_ticklabels(), visible=False)
        self.ref_ax.set_title(r'$\Delta \Phi$ = 1')

        # put the grid lines behind
#        [line.set_zorder(10000) for line in self.ax.lines]
        self.ax.set_axisbelow(True)

        plt.show()

    #-----------------------------------------------------------------------
    def save_figure(self, save_fn, file_format='pdf',
                    orientation='portrait', fig_dpi=None, close_plot='y'):
        """
        save_plot will save the figure to save_fn.

        Arguments:
        -----------

            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as
                            save_fn/station_name_ResPhase.file_format

                          * full path -> file will be save to the given
                            path.  If you use this option then the format
                            will be assumed to be provided by the path

            **file_format** : [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps...

            **orientation** : [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait

            **fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the dpi will be that at
                          which the figure was made.  I don't think that
                          it can be larger than dpi of the figure.

            **close_plot** : [ y | n ]
                             * 'y' will close the plot after saving.
                             * 'n' will leave plot open

        :Example: ::

            >>> ptr.save_figure()

        """

        if fig_dpi is None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')
            plt.clf()
            plt.close(self.fig)

        else:
            if not os.path.exists(save_fn):
                os.mkdir(save_fn)
            if not os.path.exists(os.path.join(save_fn, 'RPT_PS')):
                os.mkdir(os.path.join(save_fn, 'RPT_PS'))
                save_fn = os.path.join(save_fn, 'RPT_PS')

            save_fn = os.path.join(save_fn, 'RPT_PS_{0}.{1}'.format(
                                   self.ellipse_colorby, file_format))
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        if close_plot == 'y':
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print('Saved figure to: ' + self.fig_fn)

    #-----------------------------------------------------------------------
    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.

        Use this if you change an of the .fig or axes properties

        :Example: ::

            >>> # to change the grid lines to only be on the major ticks
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotResPhase(r'/home/MT/mt01.edi')
            >>> [ax.grid(True, which='major') for ax in [p1.axr,p1.axp]]
            >>> p1.update_plot()

        """

        self.fig.canvas.draw()

   #-------------------------------------------------------------------------
    def redraw_plot(self):
        """
        use this function if you updated some attributes and want to re-plot.

        :Example: ::

            >>> # change the color and marker of the xy components
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotResPhase(r'/home/MT/mt01.edi')
            >>> p1.xy_color = (.5,.5,.9)
            >>> p1.xy_marker = '*'
            >>> p1.redraw_plot()
        """

        plt.close(self.fig)
        self.plot()
