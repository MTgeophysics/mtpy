# -*- coding: utf-8 -*-
"""
Created on Thu May 30 18:10:55 2013

@author: jpeacock-pr
"""

# ==============================================================================

import os

import matplotlib.colorbar as mcb
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

import mtpy.imaging.mtcolors as mtcl
import mtpy.imaging.mtplottools as mtpl


# ==============================================================================

class PlotPhaseTensorPseudoSection(mtpl.PlotSettings):
    """
    PlotPhaseTensorPseudoSection will plot the phase tensor ellipses in a
    pseudo section format


    Arguments:
    ------------

        **fn_list** : list of strings
                          full paths to .edi files to plot

        **z_object** : class mtpy.core.z.Z
                      object of mtpy.core.z.  If this is input be sure the
                      attribute z.frequency is filled.  *default* is None

        **mt_object** : class mtpy.imaging.mtplot.MTplot
                        object of mtpy.imaging.mtplot.MTplot
                        *default* is None

        **pt_object** : class mtpy.analysis.pt
                        phase tensor object of mtpy.analysis.pt.  If this is
                        input then the ._mt attribute is set to None cause
                        at the moment cannot tranform the phase tensor to z
                        *default* is None

        **ellipse_dict** : dictionary
                          dictionary of parameters for the phase tensor
                          ellipses with keys:
                              * 'size' -> size of ellipse in points
                                         *default* is 2

                              * 'colorby' : [ 'phimin' | 'phimax' | 'skew' |
                                              'skew_seg' | 'phidet' |
                                              'ellipticity' ]

                                        - 'phimin' -> colors by minimum phase
                                        - 'phimax' -> colors by maximum phase
                                        - 'skew' -> colors by skew
                                        - 'skew_seg' -> colors by skew in
                                                       discrete segments
                                                       defined by the range
                                        - 'normalized_skew' -> colors by skew
                                                see [Booker, 2014]
                                        - 'normalized_skew_seg' -> colors by
                                                       normalized skew in
                                                       discrete segments
                                                       defined by the range
                                        - 'phidet' -> colors by determinant of
                                                     the phase tensor
                                        - 'ellipticity' -> colors by ellipticity
                                        *default* is 'phimin'

                               * 'range' : tuple (min, max, step)
                                     Need to input at least the min and max
                                     and if using 'skew_seg' to plot
                                     discrete values input step as well
                                     *default* depends on 'colorby'

                          * 'cmap' : [ 'mt_yl2rd' | 'mt_bl2yl2rd' |
                                      'mt_wh2bl' | 'mt_rd2bl' |
                                      'mt_bl2wh2rd' | 'mt_seg_bl2wh2rd' |
                                      'mt_rd2gr2bl' ]

                                   - 'mt_yl2rd' -> yellow to red
                                   - 'mt_bl2yl2rd' -> blue to yellow to red
                                   - 'mt_wh2bl' -> white to blue
                                   - 'mt_rd2bl' -> red to blue
                                   - 'mt_bl2wh2rd' -> blue to white to red
                                   - 'mt_bl2gr2rd' -> blue to green to red
                                   - 'mt_rd2gr2bl' -> red to green to blue
                                   - 'mt_seg_bl2wh2rd' -> discrete blue to
                                                         white to red



        **stretch** : float or tuple (xstretch, ystretch)
                        is a factor that scales the distance from one
                        station to the next to make the plot readable.
                        *Default* is 200

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


        **fignum** : int
                     figure number.  *Default* is 1

        **plot_tipper** : [ 'yri' | 'yr' | 'yi' | 'n' ]
                        * 'yri' to plot induction both real and imaginary
                           induction arrows

                        * 'yr' to plot just the real induction arrows

                        * 'yi' to plot the imaginary induction arrows

                        * 'n' to not plot them

                        *Default* is 'n'

                        **Note: convention is to point towards a conductor but
                        can be changed in arrow_dict['direction']**

        **arrow_dict** : dictionary for arrow properties
                        * 'size' : float
                                  multiplier to scale the arrow. *default* is 5
                        * 'head_length' : float
                                         length of the arrow head *default* is
                                         1.5
                        * 'head_width' : float
                                        width of the arrow head *default* is
                                        1.5
                        * 'lw' : float
                                line width of the arrow *default* is .5

                        * 'color' : tuple (real, imaginary)
                                   color of the arrows for real and imaginary

                        * 'threshold': float
                                      threshold of which any arrow larger than
                                      this number will not be plotted, helps
                                      clean up if the data is not good.
                                      *default* is 1, note this is before
                                      scaling by 'size'

                        * 'direction : [ 0 | 1 ]
                                     -0 for arrows to point toward a conductor
                                     -1 for arrow to point away from conductor

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

        **xlim** : tuple(xmin, xmax)
                   min and max along the x-axis in relative distance of degrees
                   and multiplied by xstretch

        **ylim** : tuple(ymin, ymax)
                   min and max period to plot, note that the scaling will be
                   done in the code.  So if you want to plot from (.1s, 100s)
                   input ylim=(.1,100)

    To get a list of .edi files that you want to plot -->
    :Example: ::

        >>> import mtpy.imaging.mtplot as mtplot
        >>> import os
        >>> edipath = r"/home/EDIfiles"
        >>> edilist = [os.path.join(edipath,edi) for edi in os.listdir(edipath)
        >>> ...       if edi.find('.edi')>0]

    * If you want to plot minimum phase colored from blue to red in a range of
     20 to 70 degrees you can do it one of two ways-->

    1)
    :Example: ::

        >>> edict = {'range':(20,70), 'cmap':'mt_bl2gr2rd','colorby':'phimin'}
        >>> pt1 = mtplot.plot_pt_pseudosection(fn_list=edilist,
                                               ellipse_dict=edict)

    2)
    :Example: ::

        >>> pt1 = mtplot.plot_pt_pseudosection(fn_list=edilist, plot_yn='n')
        >>> pt1.ellipse_colorby = 'phimin'
        >>> pt1.ellipse_cmap = 'mt_bl2gr2rd'
        >>> pt1.ellipse_range = (20,70)
        >>> pt1.plot()

    * If you want to add real induction arrows that are scaled by 10 and point
     away from a conductor -->
    :Example: ::

        >>> pt1.plot_tipper = 'yr'
        >>> pt1.arrow_size = 10
        >>> pt1.arrow_direction = -1
        >>> pt1.redraw_plot()

    * If you want to save the plot as a pdf with a generic name -->
    :Example: ::
        >>> pt1.save_figure(r"/home/PTFigures", file_format='pdf', dpi=300)
        File saved to '/home/PTFigures/PTPseudoSection.pdf'

    Attributes:
    -------------
        -arrow_color_imag     color of imaginary induction arrow
        -arrow_color_real     color of real induction arrow
        -arrow_direction      convention of arrows pointing to or away from
                              conductors, see above.
        -arrow_head_length    length of arrow head in relative points
        -arrow_head_width     width of arrow head in relative points
        -arrow_lw             line width of arrows
        -arrow_size           scaling factor to multiple arrows by to be visible
        -arrow_threshold      threshold for plotting arrows, anything above
                              this number will not be plotted.

        -ax                   matplotlib.axes instance for the main plot
        -ax2                  matplotlib.axes instance for the color bar
        -cb                   matplotlib.colors.ColorBar instance for color bar
        -cb_orientation       color bar orientation ('vertical' | 'horizontal')
        -cb_position          color bar position (x, y, dx, dy)

        -dpi                  dots-per-inch resolution

        -ellipse_cmap         ellipse color map, see above for options
        -ellipse_colorby      parameter to color ellipse by
        -ellipse_range        (min, max, step) values to color ellipses
        -ellipse_size         scaling factor to make ellipses visible

        -fig                  matplotlib.figure instance for the figure
        -fignum               number of figure being plotted
        -figsize              size of figure in inches
        -font_size            font size of axes tick label, axes labels will be
                              font_size + 2

        -linedir              prominent direction of profile being plotted

        -mt_list               list of mtplot.MTplot instances containing all
                              the important information for each station
        -offsetlist            array of relative offsets of each station

        -plot_tipper          string to inform program to plot induction arrows
        -plot_yn              plot the pseudo section on instance creation

        -rot_z                rotates the data by this angle assuming North is
                              0 and angle measures clockwise

        -station_id            index [min, max] to reaad station name
        -stationlist           list of stations plotted
        -title                title of figure
        -tscale               temporal scale of y-axis ('frequency' | 'period')

        -xlimits              limits on x-axis (xmin, xmax)
        -xstretch             scaling factor to stretch x offsets

        -ylimits              limits on y-axis (ymin, ymax)
        -ystep                step to set major ticks on y-axis
        -ystretch             scaling factor to strech axes in y direction

    Methods:
    ----------

        -plot                 plots the pseudo section
        -redraw_plot          on call redraws the plot from scratch
        -save_figure          saves figure to a file of given format
        -update_plot          updates the plot while still active
        -export_pt_params_to_file       writes parameters of the phase tensor and tipper
                              to text files.

    """

    def __init__(self, **kwargs):
        super(PlotPhaseTensorPseudoSection, self).__init__()
        mtpl.PlotSettings.__init__(self)
        self._rotation_angle = 0

        fn_list = kwargs.pop('fn_list', None)
        z_object_list = kwargs.pop('z_object_list', None)
        tipper_object_list = kwargs.pop('tipper_object_list', None)
        mt_object_list = kwargs.pop('mt_object_list', None)
        res_object_list = kwargs.pop('res_object_list', None)

        # ----set attributes for the class-------------------------
        self.mt_list = mtpl.get_mtlist(fn_list=fn_list,
                                       res_object_list=res_object_list,
                                       z_object_list=z_object_list,
                                       tipper_object_list=tipper_object_list,
                                       mt_object_list=mt_object_list)

        # --> set the ellipse properties
        self._ellipse_dict = kwargs.pop('ellipse_dict', {})
        self._read_ellipse_dict(self._ellipse_dict)

        # --> set colorbar properties
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

        # set the stretching in each direction
        stretch = kwargs.pop('stretch', (200, 25))
        if isinstance(stretch, float) or isinstance(stretch, int):
            self.xstretch = stretch
            self.ystretch = stretch
        else:
            self.xstretch = stretch[0]
            self.ystretch = stretch[1]

        # --> set plot properties
        self.fig_num = kwargs.pop('fig_num', 1)
        self.plot_num = kwargs.pop('plot_num', 1)
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
        self.scale_arrow = kwargs.pop('scale_arrow', False)
        self.scale_arrow_dict = kwargs.pop('scale_arrow_dict', {})

        if 'size' not in list(self.scale_arrow_dict.keys()):
            self.scale_arrow_dict['size'] = 1.
        if 'text_offset_y' not in list(self.scale_arrow_dict.keys()):
            self.scale_arrow_dict['text_offset_y'] = 0.

        self._rot_z = kwargs.pop('rot_z', 0)
        if isinstance(self._rot_z, float) or isinstance(self._rot_z, int):
            self._rot_z = np.array([self._rot_z] * len(self.mt_list))

        # if the rotation angle is an array for rotation of different
        # freq than repeat that rotation array to the len(mt_list)
        elif isinstance(self._rot_z, np.ndarray):
            if self._rot_z.shape[0] != len(self.mt_list):
                self._rot_z = np.repeat(self._rot_z, len(self.mt_list))

        else:
            pass

        # --> set induction arrow properties -------------------------------
        self.plot_tipper = kwargs.pop('plot_tipper', 'n')

        self._arrow_dict = kwargs.pop('arrow_dict', {})
        self._read_arrow_dict(self._arrow_dict)
        
        self.subplot_left = .10
        self.subplot_right = .90
        self.subplot_bottom = .2
        self.subplot_top = 0.9
        self.subplot_wspace = .05
        self.subplot_hspace = .05
        
        for key, value in kwargs.items():
            setattr(self, key, value)

    #---need to rotate data on setting rotz
    @property
    def rotation_angle(self):
        return self._rotation_angle
    
    @rotation_angle.setter
    def rotation_angle(self, value):
        """
        only a single value is allowed
        """
        for ii, mt in enumerate(self.mt_list):
            # JP: need to set the rotation angle negative for plotting
            # I think its because the way polar plots work by measuring 
            # counter clockwise
            mt.rotation_angle = value
            
        self._rotation_angle = value
            

    def plot(self, show=True):
        """
        plots the phase tensor pseudo section.  See class doc string for
        more details.
        """

        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        plt.rcParams['figure.subplot.wspace'] = self.subplot_wspace
        plt.rcParams['figure.subplot.hspace'] = self.subplot_hspace

        # create a plot instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        self.ax = self.fig.add_subplot(1, 1, 1, aspect='equal')

        # FZ: control tick rotation=30 not that good
        plt.xticks(rotation='vertical')

        # create empty lists to put things into
        self.stationlist = []
        self.offsetlist = []
        minlist = []
        maxlist = []
        plot_periodlist = None

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

        nseg = float((ckmax - ckmin) / (2 * ckstep))

        if cmap == 'mt_seg_bl2wh2rd':
            bounds = np.arange(ckmin, ckmax + ckstep, ckstep)
        # plot phase tensor ellipses
        for ii, mt in enumerate(self.mt_list):
            self.stationlist.append(
                mt.station[self.station_id[0]:self.station_id[1]])

            # set the an arbitrary origin to compare distance to all other
            # stations.
            if ii == 0:
                east0 = mt.lon
                north0 = mt.lat
                offset = 0.0
            else:
                east = mt.lon
                north = mt.lat
                if self.linedir == 'ew':
                    if east0 < east:
                        offset = np.sqrt((east0 - east) **
                                         2 + (north0 - north) ** 2)
                    elif east0 > east:
                        offset = -1 * \
                            np.sqrt((east0 - east) ** 2 +
                                    (north0 - north) ** 2)
                    else:
                        offset = 0
                elif self.linedir == 'ns':
                    if north0 < north:
                        offset = np.sqrt((east0 - east) **
                                         2 + (north0 - north) ** 2)
                    elif north0 > north:
                        offset = -1 * \
                            np.sqrt((east0 - east) ** 2 +
                                    (north0 - north) ** 2)
                    else:
                        offset = 0

            self.offsetlist.append(offset)

            # get phase tensor elements and flip so the top is small
            # periods/high frequency
            pt = mt.pt

            periodlist = mt.period[::-1]
            phimax = pt.phimax[::-1]
            phimin = pt.phimin[::-1]
            azimuth = pt.azimuth[::-1]

            # if there are induction arrows, flip them as pt
            if self.plot_tipper.find('y') == 0:
                tip = mt.Tipper
                if tip.mag_real is not None:
                    tmr = tip.mag_real[::-1]
                    tmi = tip.mag_imag[::-1]
                    tar = tip.angle_real[::-1]
                    tai = tip.angle_imag[::-1]
                else:
                    tmr = np.zeros(len(mt.period))
                    tmi = np.zeros(len(mt.period))
                    tar = np.zeros(len(mt.period))
                    tai = np.zeros(len(mt.period))

                aheight = self.arrow_head_length
                awidth = self.arrow_head_width
                alw = self.arrow_lw

            # get the properties to color the ellipses by
            if self.ellipse_colorby == 'phimin':
                colorarray = pt.phimin[::-1]

            elif self.ellipse_colorby == 'phimax':
                colorarray = pt.phimin[::-1]

            elif self.ellipse_colorby == 'phidet':
                colorarray = np.sqrt(abs(pt.det[::-1])) * (180 / np.pi)

            elif self.ellipse_colorby == 'skew' or \
                    self.ellipse_colorby == 'skew_seg':
                colorarray = pt.beta[::-1]

            elif self.ellipse_colorby == 'normalized_skew' or \
                    self.ellipse_colorby == 'normalized_skew_seg':
                colorarray = 2 * pt.beta[::-1]

            elif self.ellipse_colorby == 'ellipticity':
                colorarray = pt.ellipticity[::-1]
            elif self.ellipse_colorby in ['strike', 'azimuth']:
                colorarray = pt.azimuth[::-1] % 180
            else:
                raise NameError(self.ellipse_colorby + ' is not supported')

            # get the number of periods
            n = len(periodlist)

            if ii == 0:
                plot_periodlist = periodlist

            else:
                if n > len(plot_periodlist):
                    plot_periodlist = periodlist

            # get min and max of the color array for scaling later
            minlist.append(min(colorarray))
            maxlist.append(max(colorarray))

            for jj, ff in enumerate(periodlist):

                # make sure the ellipses will be visable
                eheight = phimin[jj] / phimax[jj] * es
                ewidth = phimax[jj] / phimax[jj] * es

                # create an ellipse scaled by phimin and phimax and orient
                # the ellipse so that north is up and east is right
                # need to add 90 to do so instead of subtracting
                ellipd = patches.Ellipse((offset * self.xstretch,
                                          np.log10(ff) * self.ystretch),
                                         width=ewidth,
                                         height=eheight,
                                         edgecolor='k',
                                         lw=0.5,
                                         angle=azimuth[jj] + 90)

                # get ellipse color
                if cmap.find('seg') > 0:
                    ellipd.set_facecolor(mtcl.get_plot_color(colorarray[jj],
                                                             self.ellipse_colorby,
                                                             cmap,
                                                             ckmin,
                                                             ckmax,
                                                             bounds=bounds))
                else:
                    ellipd.set_facecolor(mtcl.get_plot_color(colorarray[jj],
                                                             self.ellipse_colorby,
                                                             cmap,
                                                             ckmin,
                                                             ckmax))

                # == =add the ellipse to the plot == ========
                self.ax.add_artist(ellipd)

                # --------- Add induction arrows if desired -------------------
                if self.plot_tipper.find('y') == 0:

                    # --> plot real tipper
                    if self.plot_tipper == 'yri' or self.plot_tipper == 'yr':
                        txr = tmr[jj] * np.sin(tar[jj] * np.pi / 180 +
                                               np.pi * self.arrow_direction) * \
                            self.arrow_size
                        tyr = -tmr[jj] * np.cos(tar[jj] * np.pi / 180 +
                                                np.pi * self.arrow_direction) * \
                            self.arrow_size

                        maxlength = np.sqrt((txr / self.arrow_size) ** 2 +
                                            (tyr / self.arrow_size) ** 2)

                        if maxlength > self.arrow_threshold:
                            pass
                        elif ((txr > 0) or (tyr>0)):
#                        else:
                            self.ax.arrow(offset * self.xstretch,
                                          np.log10(ff) * self.ystretch,
                                          txr,
                                          tyr,
                                          lw=alw,
                                          facecolor=self.arrow_color_real,
                                          edgecolor=self.arrow_color_real,
                                          length_includes_head=False,
                                          head_width=awidth,
                                          head_length=aheight)

                    # --> plot imaginary tipper
                    if self.plot_tipper == 'yri' or self.plot_tipper == 'yi':
                        txi = tmi[jj] * np.sin(tai[jj] * np.pi / 180 +
                                               np.pi * self.arrow_direction) * \
                            self.arrow_size
                        tyi = -tmi[jj] * np.cos(tai[jj] * np.pi / 180 +
                                                np.pi * self.arrow_direction) * \
                            self.arrow_size

                        maxlength = np.sqrt((txi / self.arrow_size) ** 2 +
                                            (tyi / self.arrow_size) ** 2)
                        if maxlength > self.arrow_threshold:
                            pass
#                        else:
                        elif ((txr > 0) or (tyr>0)):
                            self.ax.arrow(offset * self.xstretch,
                                          np.log10(ff) * self.ystretch,
                                          txi,
                                          tyi,
                                          lw=alw,
                                          facecolor=self.arrow_color_imag,
                                          edgecolor=self.arrow_color_imag,
                                          length_includes_head=False,
                                          head_width=awidth,
                                          head_length=aheight)

        # --> Set plot parameters
        self._plot_periodlist = plot_periodlist
        n = len(plot_periodlist)

        # calculate minimum period and maximum period with a stretch factor
        #        pmin = np.log10(plot_periodlist.min())*self.ystretch
        #        pmax = np.log10(plot_periodlist.max())*self.ystretch

        pmin = int(np.floor(np.log10(plot_periodlist.min())))
        pmax = int(np.ceil(np.log10(plot_periodlist.max())))

        # need to sort the offsets and station labels so they plot correctly
        sdtype = [('offset', np.float), ('station', 'U10')]
        slist = np.array([(oo, ss) for oo, ss in zip(self.offsetlist,
                                                     self.stationlist)], dtype=sdtype)
        offset_sort = np.sort(slist, order='offset')

        self.offsetlist = offset_sort['offset']
        self.stationlist = offset_sort['station']
        #        if self.offsetlist[0] > 0:
        #            print 'rotating'
        #            print self.stationlist
        #            self.stationlist = self.stationlist[::-1]

        # set y-ticklabels
        if self.tscale == 'period':
            yticklabels = [mtpl.labeldict[ii]
                           for ii in range(pmin, pmax + 1, 1)]
            #            yticklabels = ['{0:>4}'.format('{0: .1e}'.format(plot_period_list[ll]))
            #                            for ll in np.arange(0, n, self.ystep)]+\
            #                        ['{0:>4}'.format('{0: .1e}'.format(plot_period_list[-1]))]

            self.ax.set_ylabel('Period (s)',
                               fontsize=self.font_size + 2,
                               fontweight='bold')

        elif self.tscale == 'frequency':
            #            yticklabels = ['{0:>4}'.format('{0: .1e}'.format(1./plot_period_list[ll]))
            #                            for ll in np.arange(0, n, self.ystep)]+\
            #                            ['{0:>4}'.format('{0: .1e}'.format(1./plot_period_list[-1]))]
            #
            yticklabels = [mtpl.labeldict[-ii]
                           for ii in range(pmin, pmax + 1, 1)]
            self.ax.set_ylabel('Frequency (Hz)',
                               fontsize=self.font_size + 2,
                               fontweight='bold')
        # set x-axis label
        self.ax.set_xlabel('Station',
                           fontsize=self.font_size + 2,
                           fontweight='bold')

        # --> set tick locations and labels
        # set y-axis major ticks
        #        self.ax.yaxis.set_ticks([np.log10(plot_periodlist[ll])*self.ystretch
        #                             for ll in np.arange(0, n, self.ystep)]):
        self.ax.yaxis.set_ticks(np.arange(pmin * self.ystretch,
                                          (pmax + 1) * self.ystretch,
                                          self.ystretch))

        # set y-axis minor ticks
        #        self.ax.yaxis.set_ticks([np.log10(plot_periodlist[ll])*self.ystretch
        #                             for ll in np.arange(0, n, 1)],minor=True)
        # set y-axis tick labels
        self.ax.set_yticklabels(yticklabels)

        # set x-axis ticks
        self.ax.set_xticks(self.offsetlist * self.xstretch)

        # set x-axis tick labels as station names
        xticklabels = self.stationlist
        if self.xstep != 1:
            xticklabels = np.zeros(len(self.stationlist),
                                   dtype=self.stationlist.dtype)
            for xx in range(0, len(self.stationlist), self.xstep):
                xticklabels[xx] = self.stationlist[xx]
        self.ax.set_xticklabels(xticklabels)

        # --> set x-limits
        if self.xlimits is None:
            self.ax.set_xlim(self.offsetlist.min() * self.xstretch - es * 2,
                             self.offsetlist.max() * self.xstretch + es * 2)
        else:
            self.ax.set_xlim(self.xlimits)

        # --> set y-limits
        if self.ylimits is None:
            #            self.ax.set_ylim(pmax+es*2, pmin-es*2)
            self.ax.set_ylim(pmax * self.ystretch, pmin * self.ystretch)
        else:
            pmin = np.log10(self.ylimits[0]) * self.ystretch
            pmax = np.log10(self.ylimits[1]) * self.ystretch
            self.ax.set_ylim(pmax + es * 2, pmin - es * 2)
        #            self.ax.set_ylim(pmax, pmin)

        # --> set title of the plot
        if self.plot_title is None:
            pass
        else:
            self.ax.set_title(self.plot_title, fontsize=self.font_size + 2)

        # make a legend for the induction arrows
        if self.plot_tipper.find('y') == 0:
            if self.plot_tipper == 'yri':
                treal = self.ax.plot(np.arange(10) * .000005,
                                     np.arange(10) * .00005,
                                     color=self.arrow_color_real)
                timag = self.ax.plot(np.arange(10) * .000005,
                                     np.arange(10) * .00005,
                                     color=self.arrow_color_imag)
                self.ax.legend([treal[0], timag[0]],
                               ['Tipper_real', 'Tipper_imag'],
                               loc='lower right',
                               prop={
                    'size': self.font_size - 1,
                    'weight': 'bold'},
                    ncol=2,
                    markerscale=.5,
                    borderaxespad=.005,
                    borderpad=.25)

            elif self.plot_tipper == 'yr':
                treal = self.ax.plot(np.arange(10) * .000005,
                                     np.arange(10) * .00005,
                                     color=self.arrow_color_real)
                self.ax.legend([treal[0]],
                               ['Tipper_real'],
                               loc='lower right',
                               prop={
                    'size': self.font_size - 1,
                    'weight': 'bold'},
                    ncol=2,
                    markerscale=.5,
                    borderaxespad=.005,
                    borderpad=.25)

            elif self.plot_tipper == 'yi':
                timag = self.ax.plot(np.arange(10) * .000005,
                                     np.arange(10) * .00005,
                                     color=self.arrow_color_imag)
                self.ax.legend([timag[0]],
                               ['Tipper_imag'],
                               loc='lower right',
                               prop={
                    'size': self.font_size - 1,
                    'weight': 'bold'},
                    ncol=2,
                    markerscale=.5,
                    borderaxespad=.005,
                    borderpad=.25)

            # make a scale arrow
            if self.scale_arrow:
                print((
                    np.log10(
                        self.ylimits[1] - self.scale_arrow_dict['text_offset_y'])) * self.ystretch)
                txrl = self.scale_arrow_dict['size']
                self.ax.arrow(min(self.offsetlist) * self.xstretch,
                              np.log10(self.ylimits[1]) * self.ystretch,
                              txrl * self.arrow_size,
                              0.,
                              lw=alw,
                              facecolor=self.arrow_color_real,
                              edgecolor=self.arrow_color_real,
                              length_includes_head=False,
                              head_width=awidth,
                              head_length=aheight)
                self.ax.text(min(self.offsetlist) * self.xstretch,
                             (np.log10(self.ylimits[
                              1] - self.scale_arrow_dict['text_offset_y'])) * self.ystretch,
                             '|T| = %3.1f' % txrl)

        # put a grid on the plot
        self.ax.grid(alpha=.25, which='both', color=(.25, .25, .25))

        # print out the min an max of the parameter plotted
        print('-' * 25)
        print(ck + ' min = {0:.2f}'.format(min(minlist)))
        print(ck + ' max = {0:.2f}'.format(max(maxlist)))
        print('-' * 25)

        # ==> make a colorbar with appropriate colors
        if self.cb_position is None:
            self.ax2, kw = mcb.make_axes(self.ax,
                                         orientation=self.cb_orientation,
                                         shrink=.35)
        else:
            self.ax2 = self.fig.add_axes(self.cb_position)

        if cmap == 'mt_seg_bl2wh2rd':
            # make a color list
            self.clist = [(cc, cc, 1)
                          for cc in np.arange(0, 1 + 1. / (nseg), 1. / (nseg))] + \
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

        # --> add reference ellipse
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
        self.ref_ax.set_title(r'$\Phi$ = 1')

        # put the grid lines behind
        #        [line.set_zorder(10000) for line in self.ax.lines]
        self.ax.set_axisbelow(True)

        if show:
            plt.show()

    def writeTextFiles(self, save_path=None, ptol=0.10):
        """
        This will write text files for all the phase tensor parameters
        """

        if save_path is None:
            try:
                svpath = os.path.dirname(self.mt_list[0].fn)
            except TypeError:
                raise IOError('Need to input save_path, could not find a path')
        else:
            svpath = save_path

        # check to see if plot has been run if not run it
        try:
            plist = self._plot_periodlist

        except AttributeError:
            self.plot()
            plist = self._plot_periodlist

        if plist[0] > plist[-1]:
            plist = plist[::-1]

        if self.tscale == 'frequency':
            plist = 1. / plist

        # match station list with mt list
        slist = [mt for ss in self.stationlist for mt in self.mt_list
                 if os.path.basename(mt.fn).find(ss) >= 0]

        ns = len(slist) + 1
        nt = len(plist) + 1

        # set some empty lists to put things into
        sklist = np.zeros((nt, ns), dtype='|S8')
        phiminlist = np.zeros((nt, ns), dtype='|S8')
        phimaxlist = np.zeros((nt, ns), dtype='|S8')
        elliplist = np.zeros((nt, ns), dtype='|S8')
        azimlist = np.zeros((nt, ns), dtype='|S8')
        tiplistr = np.zeros((nt, ns), dtype='|S8')
        tiplisti = np.zeros((nt, ns), dtype='|S8')
        tiplistraz = np.zeros((nt, ns), dtype='|S8')
        tiplistiaz = np.zeros((nt, ns), dtype='|S8')

        sklist[0, 0] = '{0:>8} '.format(self.tscale)
        phiminlist[0, 0] = '{0:>8} '.format(self.tscale)
        phimaxlist[0, 0] = '{0:>8} '.format(self.tscale)
        elliplist[0, 0] = '{0:>8} '.format(self.tscale)
        azimlist[0, 0] = '{0:>8} '.format(self.tscale)
        tiplistr[0, 0] = '{0:>8} '.format(self.tscale)
        tiplistraz[0, 0] = '{0:>8} '.format(self.tscale)
        tiplisti[0, 0] = '{0:>8} '.format(self.tscale)
        tiplistiaz[0, 0] = '{0:>8} '.format(self.tscale)

        # get the period as the first column
        for tt, t1 in enumerate(plist, 1):
            sklist[tt, 0] = t1
            phiminlist[tt, 0] = t1
            phimaxlist[tt, 0] = t1
            elliplist[tt, 0] = t1
            azimlist[tt, 0] = t1
            tiplistr[tt, 0] = t1
            tiplistraz[tt, 0] = t1
            tiplisti[tt, 0] = t1
            tiplistiaz[tt, 0] = t1

        # fill out the rest of the values
        for kk, mt in enumerate(slist, 1):

            pt = mt.get_PhaseTensor()
            tip = mt.Tipper

            if self.tscale == 'period':
                tlist = mt.period

            elif self.tscale == 'frequency':
                tlist = mt.frequency

            try:
                stationstr = '{0:^8}'.format(mt.station[self.station_id[0]:
                                                        self.station_id[1]])
            except AttributeError:
                stationstr = '{0:^8}'.format(mt.station)

            # -->  get station name as header in each file
            sklist[0, kk] = stationstr
            phiminlist[0, kk] = stationstr
            phimaxlist[0, kk] = stationstr
            elliplist[0, kk] = stationstr
            azimlist[0, kk] = stationstr
            tiplistr[0, kk] = stationstr
            tiplistraz[0, kk] = stationstr
            tiplisti[0, kk] = stationstr
            tiplistiaz[0, kk] = stationstr

            # If the all periods match for the station and the plotting period
            if tlist.all() == plist.all():
                if pt.pt is not None:
                    sklist[1:, kk] = pt.beta[0]
                    phiminlist[1:, kk] = pt.phimin[0]
                    phimaxlist[1:, kk] = pt.phimax[0]
                    elliplist[1:, kk] = pt.ellipticity[0]
                    azimlist[1:, kk] = pt.azimuth[0]
                if tip.mag_real is not None:
                    tiplistr[1:, kk] = tip.mag_real
                    tiplistraz[1:, kk] = tip.angle_real
                    tiplisti[1:, kk] = tip.mag_imag
                    tiplistiaz[1:, kk] = tip.angle_imag

            # otherwise search the period list to find a cooresponding period
            else:
                for mm, t1 in enumerate(plist):
                    # check to see if the periods match or are at least close in
                    # case there are frequency missing
                    t1_yn = False
                    if t1 == tlist[mm]:
                        t1_yn = True
                    elif tlist[mm] > t1 * (1 - ptol) and tlist[mm] < t1 * (1 + ptol):
                        t1_yn = True

                    if t1_yn == True:
                        # add on the value to the present row
                        if pt.beta[0] is not None:
                            sklist[mm + 1, kk] = pt.beta[0][mm]
                            phiminlist[mm + 1, kk] = pt.phimin[0][mm]
                            phimaxlist[mm + 1, kk] = pt.phimax[0][mm]
                            elliplist[mm + 1, kk] = pt.ellipticity[0][mm]
                            azimlist[mm + 1, kk] = pt.azimuth[0][mm]

                        # add on the value to the present row
                        if tip.mag_real is not None:
                            tiplistr[mm + 1, kk] = tip.mag_real[mm]
                            tiplistraz[mm + 1, kk] = tip.angle_real[mm]
                            tiplisti[mm + 1, kk] = tip.mag_imag[mm]
                            tiplistiaz[mm + 1, kk] = tip.angle_imag[mm]

                    elif t1_yn == False:
                        for ff, t2 in enumerate(tlist):
                            if t2 > t1 * (1 - ptol) and t2 < t1 * (1 + ptol):
                                # add on the value to the present row
                                if pt.beta[0] is not None:
                                    sklist[mm + 1, kk] = pt.beta[0][ff]
                                    phiminlist[mm + 1, kk] = pt.phimin[0][ff]
                                    phimaxlist[mm + 1, kk] = pt.phimax[0][ff]
                                    elliplist[
                                        mm + 1, kk] = pt.ellipticity[0][ff]
                                    azimlist[mm + 1, kk] = pt.azimuth[0][ff]

                                # add on the value to the present row
                                if tip.mag_real is not None:
                                    tiplistr[mm + 1, kk] = tip.mag_real[ff]
                                    tiplistraz[mm + 1, kk] = tip.angle_real[ff]
                                    tiplisti[mm + 1, kk] = tip.mag_imag[ff]
                                    tiplistiaz[mm + 1, kk] = tip.angle_imag[ff]
                                t1_yn = True
                                break
                            else:
                                t1_yn = False

        # write the arrays into lines properly formatted
        t1_kwargs = {'spacing': '{0:^8} ', 'value_format': '{0:.2e}',
                     'append': False, 'add': False}
        t2_kwargs = {'spacing': '{0:^8}', 'value_format': '{0: .2f}',
                     'append': False, 'add': False}
        # create empty lists to put the concatenated strings into
        sklines = []
        phiminlines = []
        phimaxlines = []
        elliplines = []
        azimlines = []
        tprlines = []
        tprazlines = []
        tpilines = []
        tpiazlines = []

        # if there are any blank strings set them as 0
        sklist[np.where(sklist == '')] = '0.0'
        phiminlist[np.where(phiminlist == '')] = '0.0'
        phimaxlist[np.where(phimaxlist == '')] = '0.0'
        elliplist[np.where(elliplist == '')] = '0.0'
        azimlist[np.where(azimlist == '')] = '0.0'
        tiplistr[np.where(tiplistr == '')] = '0.0'
        tiplistraz[np.where(tiplistraz == '')] = '0.0'
        tiplisti[np.where(tiplisti == '')] = '0.0'
        tiplistiaz[np.where(tiplistiaz == '')] = '0.0'

        for tt in range(nt):
            if tt == 0:
                skline = sklist[tt, 0] + ' '
                pminline = phiminlist[tt, 0] + ' '
                pmaxline = phimaxlist[tt, 0] + ' '
                elliline = elliplist[tt, 0] + ' '
                azline = azimlist[tt, 0] + ' '
                tprline = tiplistr[tt, 0] + ' '
                tprazline = tiplistraz[tt, 0] + ' '
                tpiline = tiplisti[tt, 0] + ' '
                tpiazline = tiplistiaz[tt, 0] + ' '
                for ss in range(1, ns):
                    skline += sklist[tt, ss]
                    pminline += phiminlist[tt, ss]
                    pmaxline += phimaxlist[tt, ss]
                    elliline += elliplist[tt, ss]
                    azline += azimlist[tt, ss]
                    tprline += tiplistr[tt, ss]
                    tprazline += tiplistraz[tt, ss]
                    tpiline += tiplisti[tt, ss]
                    tpiazline += tiplistiaz[tt, ss]
            else:
                # get period or frequency
                skline = mtpl.make_value_str(float(sklist[tt, 0]),
                                             **t1_kwargs)
                pminline = mtpl.make_value_str(float(phiminlist[tt, 0]),
                                               **t1_kwargs)
                pmaxline = mtpl.make_value_str(float(phimaxlist[tt, 0]),
                                               **t1_kwargs)
                elliline = mtpl.make_value_str(float(elliplist[tt, 0]),
                                               **t1_kwargs)
                azline = mtpl.make_value_str(float(azimlist[tt, 0]),
                                             **t1_kwargs)
                tprline = mtpl.make_value_str(float(tiplistr[tt, 0]),
                                              **t1_kwargs)
                tprazline = mtpl.make_value_str(float(tiplistraz[tt, 0]),
                                                **t1_kwargs)
                tpiline = mtpl.make_value_str(float(tiplisti[tt, 0]),
                                              **t1_kwargs)
                tpiazline = mtpl.make_value_str(float(tiplistiaz[tt, 0]),
                                                **t1_kwargs)

                # get parameter values
                for ss in range(1, ns):
                    skline += mtpl.make_value_str(float(sklist[tt, ss]),
                                                  **t2_kwargs)
                    pminline += mtpl.make_value_str(float(phiminlist[tt, ss]),
                                                    **t2_kwargs)
                    pmaxline += mtpl.make_value_str(float(phimaxlist[tt, ss]),
                                                    **t2_kwargs)
                    elliline += mtpl.make_value_str(float(elliplist[tt, ss]),
                                                    **t2_kwargs)
                    azline += mtpl.make_value_str(float(azimlist[tt, ss]),
                                                  **t2_kwargs)
                    tprline += mtpl.make_value_str(float(tiplistr[tt, ss]),
                                                   **t2_kwargs)
                    tprazline += mtpl.make_value_str(float(tiplistraz[tt, ss]),
                                                     **t2_kwargs)
                    tpiline += mtpl.make_value_str(float(tiplisti[tt, ss]),
                                                   **t2_kwargs)
                    tpiazline += mtpl.make_value_str(float(tiplistiaz[tt, ss]),
                                                     **t2_kwargs)

            # be sure to end the line after each period
            sklines.append(skline + '\n')
            phiminlines.append(pminline + '\n')
            phimaxlines.append(pmaxline + '\n')
            elliplines.append(elliline + '\n')
            azimlines.append(azline + '\n')
            tprlines.append(tprline + '\n')
            tprazlines.append(tprazline + '\n')
            tpilines.append(tpiline + '\n')
            tpiazlines.append(tpiazline + '\n')

        # write files
        skfid = file(os.path.join(svpath, 'PseudoSection.skew'), 'w')
        skfid.writelines(sklines)
        skfid.close()

        phiminfid = file(os.path.join(svpath, 'PseudoSection.phimin'), 'w')
        phiminfid.writelines(phiminlines)
        phiminfid.close()

        phimaxfid = file(os.path.join(svpath, 'PseudoSection.phimax'),
                         'w')
        phimaxfid.writelines(phimaxlines)
        phimaxfid.close()

        ellipfid = file(os.path.join(svpath, 'PseudoSection.ellipticity'),
                        'w')
        ellipfid.writelines(elliplines)
        ellipfid.close()

        azfid = file(os.path.join(svpath, 'PseudoSection.azimuth'),
                     'w')
        azfid.writelines(azimlines)
        azfid.close()

        tprfid = file(os.path.join(svpath, 'PseudoSection.tipper_mag_real'),
                      'w')
        tprfid.writelines(tprlines)
        tprfid.close()

        tprazfid = file(os.path.join(svpath, 'PseudoSection.tipper_ang_real'),
                        'w')
        tprazfid.writelines(tprazlines)
        tprazfid.close()

        tpifid = file(os.path.join(svpath, 'PseudoSection.tipper_mag_imag'),
                      'w')
        tpifid.writelines(tpilines)
        tpifid.close()

        tpiazfid = file(os.path.join(svpath, 'PseudoSection.tipper_ang_imag'),
                        'w')
        tpiazfid.writelines(tpiazlines)
        tpiazfid.close()

    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.

        Use this if you change an of the .fig or axes properties

        :Example: ::

            >>> # to change the grid lines to be on the major ticks and gray
            >>> pt1.ax.grid(True, which='major', color=(.5,.5,.5))
            >>> pt1.update_plot()

        """

        self.fig.canvas.draw()

    def redraw_plot(self):
        """
        use this function if you updated some attributes and want to re-plot.

        :Example: ::

            >>> # change ellipse size and color map to be segmented for skew
            >>> pt1.ellipse_size = 5
            >>> pt1.ellipse_colorby = 'beta_seg'
            >>> pt1.ellipse_cmap = 'mt_seg_bl2wh2rd'
            >>> pt1.ellipse_range = (-9, 9, 3)
            >>> pt1.redraw_plot()
        """

        self.fig.clf()
        self.plot()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return "Plots pseudo section of phase tensor ellipses"

    def save_figure(self, save_fn, file_format='png', orientation='portrait',
                    fig_dpi=None, close_plot='y'):
        """
        save_plot will save the figure to save_fn.

        Arguments:
        -----------

            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as
                            save_fn/station_name_PTPseudoSection.file_format

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
            >>> # save plot as a jpg
            >>> pt1.save_plot(r'/home/MT/figures', file_format='jpg')

        """

        if fig_dpi is None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) is False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        else:
            save_fn = os.path.join(save_fn, '_PTPseudoSection.' +
                                   file_format)
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        if close_plot == 'y':
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print('Saved figure to: ' + self.fig_fn)

    def save_figure2(self, save_fn, file_format='jpg', orientation='portrait', fig_dpi=None, close_plot='y'):
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

            **file_format** : [ jpg | png | pdf | eps | svg ]
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
        """

        if fig_dpi is None:
            fig_dpi=self.fig_dpi

        # FZ: fixed the following logic
        if os.path.isdir(save_fn):  # FZ: assume save-fn is a directory
            if not os.path.exists(save_fn):
                os.mkdir(save_fn)

            # make a file name
            fname='PT_Pseudo_Section_DPI%s_%s.%s' % (str(self.fig_dpi), self.ellipse_colorby,  file_format)
            path2savefile=os.path.join(save_fn, fname)
            self.fig.savefig(path2savefile, dpi=fig_dpi, format=file_format, orientation=orientation,
                             bbox_inches='tight')
        else:  # FZ: assume save-fn is a path2file= "path2/afile.fmt"
            file_format=save_fn.split('.')[-1]
            if file_format is None or file_format not in ['png', 'jpg']:
                print(("Error: output file name is not correctly provided:", save_fn))
                raise Exception("output file name is not correctly provided!!!")

            path2savefile=save_fn
            self.fig.savefig(path2savefile, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')
            # plt.clf()
            # plt.close(self.fig)

        if close_plot == 'y':
            plt.clf()
            plt.close(self.fig)
        else:
            pass

        self.fig_fn=path2savefile
        #logger.debug('Saved figure to: %s', self.fig_fn)
        print(('Saved figure to: ', self.fig_fn))

        return self.fig_fn
