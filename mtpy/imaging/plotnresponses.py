# -*- coding: utf-8 -*-
"""
plots multiple MT responses simultaneously

Created on Thu May 30 17:02:39 2013
@author: jpeacock-pr

YG: the code there is massey, todo may need to rewrite it sometime

"""

# ============================================================================

import matplotlib.colorbar as mcb
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

import mtpy.imaging.mtcolors as mtcl
import mtpy.imaging.mtplottools as mtpl
from mtpy.analysis.pt import PhaseTensor
from mtpy.analysis.zinvariants import Zinvariants


# reload(mtpl)

# ============================================================================
from mtpy.core.mt import MT


class PlotMultipleResponses(mtpl.PlotSettings):
    """
    plots multiple MT responses simultaneously either in single plots or in
    one plot of sub-figures or in a single plot with subfigures for each
    component.

    expecting only one type of input --> can be:
        **fn_list** : list of filenames to plot

         **z_object_list** : list of mtpy.core.z.Z objects

         **res_object_list** : list of mtpy.imaging.mtplot.ResPhase objects

         **tipper_object_list** : list of mtpy.imaging.mtplot.Tipper objects

         **mt_object_list** : list of mtpy.imaging.mtplot.MTplot objects


    Arguments:
    ----------
        **fn_list** : list of filenames to plot
                     ie. [fn_1, fn_2, ...], *default* is None

         **z_object_list** : list of mtpy.core.z.Z objects
                            *default* is None

         **res_object_list** : list of mtpy.imaging.mtplot.ResPhase objects
                              *default* is None

         **tipper_object_list** : list of mtpy.imaging.mtplot.Tipper objects
                                 *default* is None

         **mt_object_list** : list of mtpy.imaging.mtplot.MTplot objects
                             *default* is None

        **fig_num** : int
                     figure number
                     *default* is 1
        **fig_size** : [width, height] of figure size in inches

        **rot_z** : float or np.ndarray
                   rotation angle of impedance tensor (deg or radians),
                   *Note* : rotaion is clockwise positive
                   *default* is 0
                   Can input so each station is rotated at a constant angle or
                   each period is rotated differently, or both.

        **plot_num** : [ 1 | 2 | 3 ]
                        * 1 for just Ex/By and Ey/Bx *default*
                        * 2 for all 4 components
                        * 3 for off diagonal plus the determinant

        **plot_style** : [ '1' | 'all' | 'compare' ]
                        determines the plotting style:
                            * '1' for plotting each station in a different
                                  figure. *default*

                            * 'all' for plotting each station in a subplot
                                    all in the same figure

                            * 'compare' for comparing the responses all in
                                        one plot.  Here the responses are
                                        colored from dark to light.  This
                                        plot can get messy if too many stations
                                        are plotted.


        **plot_title** : string
                    title of plot
                    *default* is station name

        **plot_tipper** : [ 'yri' | 'yr' | 'yi' | 'n' ]
                          Plots the tipper in a bottom pannel
                          * 'yri'  --> plots the real and imaginar parts
                          * 'yr'   --> plots just the real part
                          * 'yi'   --> plots just the imaginary part

                          **Note:** the convention is to point towards a
                          conductor.  Can change this by setting the
                          parameter arrow_direction = 1.

        **plot_strike** : [ 'y' | 1 | 2 | 3 | 'n' ]
                          Plots the strike angle from different parameters:
                              * 'y'  --> plots strike angle determined from
                                         the invariants of Weaver et al. [2000]
                                         and the phase tensor of
                                         Caldwell et al. [2004], if Tipper is
                                         plotted the strike of the tipper is
                                         also plotted.

                               * 1  --> plots strike angle determined from
                                        the invariants of Weaver et al. [2000]
                               * 2  --> plots strike angle determined from
                                        the phase tensor of
                                        Caldwell et al. [2004]
                               * 3  --> plots strike angle determined from
                                        the tipper
                               * 'n' --> doesn't plot the strike, *default*

        **plot_skew** : [ 'y' | 'n' ]
                       string for plotting skew angle.  This is plotted in
                       the same plot as strike angle at the moment.
                           * 'y' for plotting the skew
                           * 'n' for not plotting skew *default*

        **fig_dpi** : int
                 dots-per-inch resolution, *default* is 300


        :Example: ::

            >>> import mtpy.imaging.mtplottools as mtplot
            >>> import os
            >>> edipath = r"/home/Edifiles"
            >>> edilist = [os.path.join(edipath,edi)
            >>> ...       for edi in os.listdir(edipath)
            >>> ...       if edi.find('.edi')>0]
            >>> plot each station in a subplot all in one figure with tipper
            >>> rp1 = mtplot.PlotMultipleResPhase(fn_list=edilist, plotnum=1,
            >>> ...                                plot_tipper='yr',
            >>> ...                                plot_style='all')


    Attributes:
    -----------
        -mt_list         list of mtplot.MTplot objects made from inputs
        -fignum         figure number for plotting
        -fig_size       figure size in inches [width, height]
        -plotnum        plot type, see arguments for details
        -title          title of the plot, *default* is station name
        -dpi            Dots-per-inch resolution of plot, *default* is 300
        -rotz           Rotate impedance tensor by this angle (deg) assuming
                        that North is 0 and angle is positive clockwise

        -plot_tipper    string to tell the program to plot tipper arrows or
                        not, see accepted values above in arguments

        -plot_strike    string or integer telling the program to plot the
                        strike angle, see values above in arguments

        -plot_skew      string to tell the program to plot skew angle.
                        The skew is plotted in the same subplot as the strike
                        angle at the moment


        -period          period array cooresponding to the impedance tensor
        -font_size       size of font for the axis ticklabels, note that the
                         axis labels will be font_size+2

        -axr             matplotlib.axes object for the xy,yx resistivity plot.
        -axp             matplotlib.axes object for the xy,yx phase plot
        -axt             matplotlib.axes object for the tipper plot
        -ax2r            matplotlib.axes object for the xx,yy resistivity plot
        -ax2p            matplotlib.axes object for the xx,yy phase plot
        -axs             matplotlib.axes object for the strike plot
        -axs2            matplotlib.axes object for the skew plot

        ..

             **Note:** that from these axes object you have control of the
             plot.  You can do this by changing any parameter in the
             axes object and then calling update_plot()

        -erxyr          class matplotlib.container.ErrorbarContainer for
                        xy apparent resistivity.
        -erxyp          class matplotlib.container.ErrorbarContainer for
                        xy.
        -eryxr          class matplotlib.container.ErrorbarContainer for
                        yx apparent resistivity.
        -eryxp          class matplotlib.container.ErrorbarContainer for
                        yx phase.

        ..

            **Note:** that from these line objects you can manipulate the
            error bar properties and then call update_plot()

        -xy_ls           line style for xy and xx components, *default* is None
        -yx_ls           line style for yx and yy components, *default* is None
        -det_ls          line style for determinant, *default* is None

        -xy_marker       marker for xy and xx, *default* is squares
        -yx_marker       marker for yx and yy, *default* is circles
        -det_marker      marker for determinant, *default* is diamonds

        -xy_color        marker color for xy and xx, *default* is blue
        -yx_color        marker color for yx and yy, *default* is red
        -det_color       marker color for determinant, *default* is green

        -xy_mfc          marker face color for xy and xx, *default* is None
        -yx_mfc          marker face color for yx and yy, *default* is None
        -det_mfc         marker face color for determinant, *default* is None

        -skew_marker     marker for skew angle, *default* is 'd'
        -skew_color      color for skew angle, *default* is 'orange'

        -strike_inv_marker  marker for strike angle determined by invariants
                            *default* is '^'
        -strike_inv_color   color for strike angle determined by invaraiants
                            *default* is (.2, .2, .7)
        -strike_pt_marker  marker for strike angle determined by pt,
                           *default* is'v'
        -strike_pt_color   color for strike angle determined by pt
                           *default* is (.7, .2, .2)

        -strike_tip_marker  marker for strike angle determined by tipper
                            *default* is '>'
        -strike_tip_color   color for strike angle determined by tipper
                            *default* is (.2, .7, .2)

        -marker_size     size of marker in relative dimenstions, *default* is 2
        -marker_lw       line width of marker, *default* is 100./dpi
        ..

         *For more on line and marker styles see matplotlib.lines.Line2D*

        -arrow_lw          line width of the arrow, *default* is 0.75
        -arrow_head_width  head width of the arrow, *default* is 0 for no arrow
                           head.  Haven't found a good way to scale the arrow
                           heads in a log scale.

        -arrow_head_height  head width of the arrow, *default* is 0 for no arrow
                            head.  Haven't found a good way to scale the arrow
                            heads in a log scale.

        -arrow_color_real  color of the real arrows, *default* is black
        -arrow_color_imag  color of the imaginary arrows, *default* is blue

        -arrow_direction   0 for pointing towards a conductor and -1 for
                           pointing away from a conductor.


        -xlimits        limits on the x-limits (period), *default* is None
                        which will estimate the min and max from the data,
                        setting the min as the floor(min(period)) and the max
                        as ceil(max(period)).  Input in linear scale if you
                        want to change the period limits, ie. (.1,1000)

        -res_limits     limits on the resistivity, *default* is None, which
                        will estimate the min and max from the data, rounding
                        to the lowest and highest increments to the power of 10
                        Input in linear scale if you want to change them,
                        ie. (1,10000). Note this only sets the xy and yx
                        components, not the xx and yy.

        -phase_limits   limits on the phase, *default* is (0,90) but will
                        adapt to the data if there is phase above 90 or below
                        0.  Input in degrees.  Note this only changes the xy
                        and yx components.

        -tipper_limits  limits of the y-axis, *default* is (-1,1)

    """

    def __init__(self, **kwargs):
        """
        Initialize parameters
        """

        super(PlotMultipleResponses, self).__init__()

        fn_list = kwargs.pop('fn_list', None)
        z_object_list = kwargs.pop('z_object_list', None)
        tipper_object_list = kwargs.pop('tipper_object_list', None)
        mt_object_list = kwargs.pop('mt_object_list', None)

        # --> get the inputs into a list of mt objects
        self.mt_list = mtpl.get_mtlist(fn_list=fn_list,
                                       z_object_list=z_object_list,
                                       tipper_object_list=tipper_object_list,
                                       mt_object_list=mt_object_list)
        self.fig_num = kwargs.pop('fig_num', self.fig_num)

        # set some of the properties as attributes much to Lars' discontent
        self.plot_num = kwargs.pop('plot_num', 1)
        self.plot_style = kwargs.pop('plot_style', '1')
        self.plot_title = kwargs.pop('plot_title', None)

        # if rotation angle is an int or float make an array the length of
        # mt_list for plotting purposes
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

        self._set_rot_z(self._rot_z)

        # set plot limits
        self.xlimits = kwargs.pop('xlimits', None)
        self.res_limits = kwargs.pop('res_limits', None)
        self.phase_limits = kwargs.pop('phase_limits', None)
        self.tipper_limits = kwargs.pop('tipper_limits', None)
        self.strike_limits = kwargs.pop('strike_limits', None)
        self.skew_limits = kwargs.pop('skew_limits', (-9, 9))
        self.pt_limits = kwargs.pop('pt_limits', None)

        # set font parameters
        self.font_size = kwargs.pop('font_size', 7)

        # set plot tipper or not
        self._plot_tipper = kwargs.pop('plot_tipper', 'n')

        # plot strike angle or not
        self._plot_strike = kwargs.pop('plot_strike', 'n')

        # plot skew angle
        self._plot_skew = kwargs.pop('plot_skew', 'n')

        # plot phase tensor ellipses
        self._plot_pt = kwargs.pop('plot_pt', 'n')

        # order of plots
        self.plot_order = kwargs.pop('plot_order',
                                     ['tip', 'pt', 'strike', 'skew'])

        self.plot_dict = dict([(kk, vv) for kk, vv in zip(['tip', 'pt',
                                                           'strike', 'skew'],
                                                          [self._plot_tipper,
                                                           self._plot_pt,
                                                           self._plot_strike,
                                                           self._plot_skew])])

        # set arrow properties
        self.arrow_head_length = 0.03
        self.arrow_head_width = 0.03
        self.arrow_lw = .5

        # ellipse_properties
        self.ellipse_size = 0.25
        self.ellipse_spacing = kwargs.pop('ellipse_spacing', 1)
        if self.ellipse_size == 2 and self.ellipse_spacing == 1:
            self.ellipse_size = 0.25

        # --> set text box parameters
        self.text_location = kwargs.pop('text_location', None)
        self.text_xpad = kwargs.pop('text_xpad', 1.35)
        self.text_ypad = kwargs.pop('text_ypad', .75)
        self.text_size = kwargs.pop('text_size', 7)
        self.text_weight = kwargs.pop('text_weight', 'bold')

        self.plot_yn = kwargs.pop('plot_yn', 'y')

        # plot on initializing
        if self.plot_yn == 'y':
            self.plot()

    # ---rotate data on setting rot_z
    def _set_rot_z(self, rot_z):
        """
        need to rotate data when setting z
        """

        # if rotation angle is an int or float make an array the length of
        # mt_list for plotting purposes
        if isinstance(rot_z, float) or isinstance(rot_z, int):
            self._rot_z += np.array([rot_z] * len(self.mt_list))

        # if the rotation angle is an array for rotation of different
        # freq than repeat that rotation array to the len(mt_list)
        elif isinstance(rot_z, np.ndarray):
            if rot_z.shape[0] != len(self.mt_list):
                self._rot_z += np.repeat(rot_z, len(self.mt_list))

        else:
            pass

        for ii, mt in enumerate(self.mt_list):
            mt.rot_z = self._rot_z[ii]

    def _get_rot_z(self):
        return self._rot_z

    rot_z = property(fget=_get_rot_z, fset=_set_rot_z,
                     doc="""rotation angle(s)""")

    # --> on setting plot_ make sure to update the order and list
    def _set_plot_tipper(self, plot_tipper):
        """
        If plotting tipper make arrow attributes

        """

        self._plot_tipper = plot_tipper

        self.plot_dict['tip'] = self._plot_tipper

    def _get_plot_tipper(self):
        return self._plot_tipper

    plot_tipper = property(fget=_get_plot_tipper, fset=_set_plot_tipper,
                           doc="""string to plot tipper""")

    def _set_plot_pt(self, plot_pt):
        """
        If plotting tipper make arrow attributes

        """

        self._plot_pt = plot_pt

        self.plot_dict['pt'] = self._plot_pt

    def _get_plot_pt(self):
        return self._plot_pt

    plot_pt = property(fget=_get_plot_pt, fset=_set_plot_pt,
                       doc="""string to plot phase tensor ellipses""")

    def _set_plot_strike(self, plot_strike):
        """
        change plot_dict when changing plot_strike

        """

        self._plot_strike = plot_strike

        self.plot_dict['strike'] = self._plot_strike

    def _get_plot_strike(self):
        return self._plot_strike

    plot_strike = property(fget=_get_plot_strike, fset=_set_plot_strike,
                           doc="""string to plot strike""")

    def _set_plot_skew(self, plot_skew):
        """
        change plot_dict when changing plot_strike

        """

        self._plot_skew = plot_skew

        self.plot_dict['skew'] = self._plot_skew

    def _get_plot_skew(self):
        return self._plot_skew

    plot_skew = property(fget=_get_plot_skew, fset=_set_plot_skew,
                         doc="""string to plot skew""")

    # ---plot the resistivity and phase
    def plot(self, show=True):
        """
        plot the apparent resistivity and phase
        """
        # create a dictionary for the number of subplots needed
        pdict = {'res': 0,
                 'phase': 1}
        # start the index at 2 because resistivity and phase is permanent
        # for now
        index = 2
        for key in self.plot_order:
            if self.plot_dict[key].find('y') == 0:
                pdict[key] = index
                index += 1

        # get number of rows needed
        nrows = index

        # set height ratios of the subplots
        hr = [2, 1.5] + [1] * (len(list(pdict.keys())) - 2)

        #        if self.plot_style == '1':
        #            self.plotlist = []
        #
        #            #--> plot from edi's if given, don't need to rotate because
        #            #    data has already been rotated by the funcion _set_rot_z
        ##            if self.fig_size is None:
        ##                self.fig_size = [6, 6]
        #            for ii, mt in enumerate(self.mt_list, 1):
        #                p1 = plotresponse(mt_object=mt,
        #                                  fig_num=ii,
        #                                  fig_size=self.fig_size,
        #                                  plot_num=self.plot_num,
        #                                  fig_dpi=self.fig_dpi,
        #                                  plot_yn='n',
        #                                  plot_tipper=self._plot_tipper,
        #                                  plot_strike=self._plot_strike,
        #                                  plot_skew=self._plot_skew,
        #                                  plot_pt=self._plot_pt)
        #
        #                #make sure all the properties are set to match the users
        #                #line style between points
        #                p1.xy_ls = self.xy_ls
        #                p1.yx_ls = self.yx_ls
        #                p1.det_ls = self.det_ls
        #
        #                #outline color
        #                p1.xy_color = self.xy_color
        #                p1.yx_color = self.yx_color
        #                p1.det_color = self.det_color
        #
        #                #face color
        #                p1.xy_mfc = self.xy_mfc
        #                p1.yx_mfc = self.yx_mfc
        #                p1.det_mfc = self.det_mfc
        #
        #                #maker
        #                p1.xy_marker = self.xy_marker
        #                p1.yx_marker = self.yx_marker
        #                p1.det_marker = self.det_marker
        #
        #                #size
        #                p1.marker_size = 2
        #
        #                #set plot limits
        #                p1.xlimits = self.xlimits
        #                p1.res_limits = self.res_limits
        #                p1.phase_limits = self.phase_limits
        #
        #                #set font parameters
        #                p1.font_size = self.font_size
        #
        #                #set arrow properties
        #                p1.arrow_lw = self.arrow_lw
        #                p1.arrow_head_width = self.arrow_head_width
        #                p1.arrow_head_length = self.arrow_head_length
        #                p1.arrow_color_real = self.arrow_color_real
        #                p1.arrow_color_imag = self.arrow_color_imag
        #                p1.arrow_direction = self.arrow_direction
        #                p1.tipper_limits = self.tipper_limits
        #
        #                #skew properties
        #                p1.skew_color = self.skew_color
        #                p1.skew_marker = self.skew_marker
        #
        #                #strike properties
        #                p1.strike_inv_marker = self.strike_inv_marker
        #                p1.strike_inv_color = self.strike_inv_color
        #
        #                p1.strike_pt_marker = self.strike_pt_marker
        #                p1.strike_pt_color = self.strike_pt_color
        #
        #                p1.strike_tip_marker = self.strike_tip_marker
        #                p1.strike_tip_color = self.strike_tip_color
        #
        #                #--> plot the apparent resistivity and phase
        #                self.plotlist.append(p1)
        #
        #                p1.plot()
        #

        # -----Plot All in one figure with each plot as a subfigure------------
        if self.plot_style == 'all':

            stlist = []
            stlabel = []
            st_maxlist = []
            st_minlist = []

            ns = len(self.mt_list)

            # set some parameters of the figure and subplot spacing
            plt.rcParams['font.size'] = self.font_size
            if self.plot_skew == 'y':
                plt.rcParams['figure.subplot.right'] = .94
            else:
                plt.rcParams['figure.subplot.right'] = .98
            plt.rcParams['figure.subplot.bottom'] = .1
            plt.rcParams['figure.subplot.top'] = .93

            # set the font properties for the axis labels
            fontdict = {'size': self.font_size + 2, 'weight': 'bold'}

            # set figure size according to what the plot will be.
            if self.fig_size is None:
                if self.plot_num == 1 or self.plot_num == 3:
                    self.fig_size = [ns * 4, 6]

                elif self.plot_num == 2:
                    self.fig_size = [ns * 8, 6]

            # make a figure instance
            self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)

            # make subplots as columns for all stations that need to be plotted
            gs0 = gridspec.GridSpec(1, ns)

            # space out the subplots
            gs0.update(hspace=.025, wspace=.025, left=.085)

            labelcoords = (-0.145, 0.5)
            axr = None
            axt = None
            axpt = None
            axst = None
            axsk = None

            for ii, mt in enumerate(self.mt_list):  # type: int, MT
                # get the reistivity and phase object
                rp = mt.Z

                # set x-axis limits from short period to long period
                if self.xlimits is None:
                    self.xlimits = (10 ** (np.floor(np.log10(mt.period[0]))),
                                    10 ** (np.ceil(np.log10((mt.period[-1])))))
                # if self.phase_limits is None:
                #     pass

                if self.res_limits is None:
                    self.res_limits = (10 ** (np.floor(
                        np.log10(min([mt.Z.res_xy.min(),
                                      mt.Z.res_yx.min()])))),
                                       10 ** (np.ceil(
                                           np.log10(max([mt.Z.res_xy.max(),
                                                         mt.Z.res_yx.max()])))))

                # create a grid to place the figures into, set to have 2 rows
                # and 2 columns to put any of the 4 components.  Make the phase
                # plot slightly shorter than the apparent resistivity plot and
                # have the two close to eachother vertically.
                gs = gridspec.GridSpecFromSubplotSpec(nrows, 2,
                                                      subplot_spec=gs0[ii],
                                                      height_ratios=hr,
                                                      hspace=0.05,
                                                      wspace=.0125)

                # --> create the axes instances for xy, yx
                if self.plot_num == 1 or self.plot_num == 3:
                    # apparent resistivity axis
                    axr = self.fig.add_subplot(gs[0, :], sharex=axr)

                    # phase axis that shares period axis with resistivity
                    axp = self.fig.add_subplot(gs[1, :], sharex=axr)# --> make figure for xy,yx components
                    # space out the subplots
                    # gs.update(hspace=.05, wspace=.02, left=.1)

                # --> make figure for all 4 components
                elif self.plot_num == 2:
                    # --> create the axes instances
                    # apparent resistivity axis
                    axr = self.fig.add_subplot(gs[0, 0], sharex=axr)

                    # phase axis that shares period axis with resistivity
                    axp = self.fig.add_subplot(gs[1, 0], sharex=axr)

                    # space out the subplots
                    # gs.update(hspace=.05, wspace=.02, left=.07)


                # place y coordinate labels in the same location
                axr.yaxis.set_label_coords(labelcoords[0], labelcoords[1])
                axp.yaxis.set_label_coords(labelcoords[0], labelcoords[1])

                # --> plot tipper
                try:
                    axt = self.fig.add_subplot(gs[pdict['tip'], :], sharey=axt)
                    axt.yaxis.set_label_coords(labelcoords[0], labelcoords[1])
                except KeyError:
                    pass

                # --> plot phase tensors
                try:
                    # can't share axis because not on the same scale
                    axpt = self.fig.add_subplot(gs[pdict['pt'], :],
                                                aspect='equal', sharey=axpt)
                    axpt.yaxis.set_label_coords(labelcoords[0], labelcoords[1])
                except KeyError:
                    pass

                # --> plot strike
                try:
                    axst = self.fig.add_subplot(gs[pdict['strike'], :],
                                                sharex=axr, sharey=axst)
                    axst.yaxis.set_label_coords(labelcoords[0], labelcoords[1])
                except KeyError:
                    pass

                # --> plot skew
                try:
                    axsk = self.fig.add_subplot(gs[pdict['skew'], :],
                                                sharex=axr, sharey=axsk)
                    axsk.yaxis.set_label_coords(labelcoords[0], labelcoords[1])
                except KeyError:
                    pass

                # ---------plot the apparent resistivity----------------------
                # --> plot as error bars and just as points xy-blue, yx-red
                # res_xy
                ebxyr = axr.errorbar(mt.period,
                                     mt.Z.res_xy,
                                     marker=self.xy_marker,
                                     ms=self.marker_size,
                                     mfc=self.xy_mfc,
                                     mec=self.xy_color,
                                     mew=self.marker_lw,
                                     ls=self.xy_ls,
                                     yerr=mt.Z.res_err_xy,
                                     ecolor=self.xy_color,
                                     capsize=self.marker_size,
                                     elinewidth=self.marker_lw)

                # res_yx
                ebyxr = axr.errorbar(mt.period,
                                     mt.Z.res_yx,
                                     marker=self.yx_marker,
                                     ms=self.marker_size,
                                     mfc=self.yx_mfc,
                                     mec=self.yx_color,
                                     mew=self.marker_lw,
                                     ls=self.yx_ls,
                                     yerr=mt.Z.res_err_yx,
                                     ecolor=self.yx_color,
                                     capsize=self.marker_size,
                                     elinewidth=self.marker_lw)

                # --> set axes properties
                plt.setp(axr.get_xticklabels(), visible=False)
                axr.set_yscale('log', nonposy='clip')
                axr.set_xscale('log', nonposx='clip')
                axr.set_xlim(self.x_limits)
                axr.set_ylim(self.res_limits)
                axr.grid(True, alpha=.25, which='both',
                         color=(.25, .25, .25),
                         lw=.25)
                if ii == 0:
                    axr.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                   fontdict=fontdict)
                    axr.legend((ebxyr[0], ebyxr[0]),
                               ('$Z_{xy}$', '$Z_{yx}$'),
                               loc=3,
                               markerscale=1,
                               borderaxespad=.01,
                               labelspacing=.07,
                               handletextpad=.2,
                               borderpad=.02)
                else:
                    plt.setp(axr.get_yticklabels(), visible=False)

                # -----Plot the phase----------------------------------------
                # phase_xy
                ebxyp = axp.errorbar(mt.period,
                                     mt.Z.phase_xy,
                                     marker=self.xy_marker,
                                     ms=self.marker_size,
                                     mfc=self.xy_mfc,
                                     mec=self.xy_color,
                                     mew=self.marker_lw,
                                     ls=self.xy_ls,
                                     yerr=mt.Z.phase_err_xy,
                                     ecolor=self.xy_color,
                                     capsize=self.marker_size,
                                     elinewidth=self.marker_lw)

                # phase_yx:
                ebyxp = axp.errorbar(mt.period,
                                     mt.Z.phase_yx + 180,
                                     marker=self.yx_marker,
                                     ms=self.marker_size,
                                     mfc=self.yx_mfc,
                                     mec=self.yx_color,
                                     mew=self.marker_lw,
                                     ls=self.yx_ls,
                                     yerr=mt.Z.phase_err_yx,
                                     ecolor=self.yx_color,
                                     capsize=self.marker_size,
                                     elinewidth=self.marker_lw)

                # check the phase to see if any point are outside of [0:90]
                if self.phase_limits is None:
                    pymin = min(0, min([min(rp.phase_xy), min(rp.phase_yx)]))
                    pymax = max(89.9, max([max(rp.phase_xy), max(rp.phase_yx)]))
                    self.phase_limits = (pymin, pymax)
                #     self.phase_limits = (pymin, pymax)
                # else:
                #     self.phase_limits = (min(self.phase_limits[0], pymin),
                #                          max(self.phase_limits[1], pymax))
                # if self.phase_limits is None:
                #     if min(rp.phasexy) < 0 or min(rp.phase_yx) < 0:
                #         pymin = min([min(rp.phase_xy),
                #                      min(rp.phase_yx)])
                #         if pymin > 0:
                #             pymin = 0
                #     else:
                #         pymin = 0
                #
                #     if max(rp.phasexy) > 90 or max(rp.phase_yx) > 90:
                #         pymax = min([max(rp.phase_xy),  # YG: should use max instead ??
                #                      max(rp.phase_yx)])
                #         if pymax < 91:
                #             pymax = 89.9  # YG: why??
                #     else:
                #         pymax = 89.9
                #
                #     self.phase_limits = (pymin, pymax)

                # --> set axes properties
                if ii == 0:
                    axp.set_ylabel('Phase (deg)', fontdict)
                else:
                    plt.setp(axp.get_yticklabels(), visible=False)

                if self.plot_tipper == 'n' and self.plot_skew == 'n' and \
                        self.plot_strike == 'n':
                    axp.set_xlabel('Period (s)', fontdict)

                axp.set_xscale('log', nonposx='clip')
                if self.phase_limits is None:
                    self.phase_limits = (-179.9,179.9)
                axp.set_ylim(self.phase_limits)
                axp.set_xlim(self.x_limits)
                axp.yaxis.set_major_locator(MultipleLocator(15))
                axp.yaxis.set_minor_locator(MultipleLocator(5))
                axp.grid(True, alpha=.25, which='both',
                         color=(.25, .25, .25),
                         lw=.25)

                tklabels = [mtpl.labeldict[tt]
                            for tt in np.arange(np.log10(self.xlimits[0]),
                                                np.log10(self.xlimits[1]) + 1)]
                tklabels[0] = ''
                tklabels[-1] = ''

                axp.set_xticklabels(tklabels,
                                    fontdict={'size': self.font_size})

                if len(list(pdict.keys())) > 2:
                    plt.setp(axp.xaxis.get_ticklabels(), visible=False)
                    plt.setp(axp.xaxis.get_label(), visible=False)

                # -----plot tipper--------------------------------------------
                if self._plot_tipper.find('y') == 0:
                    plt.setp(axp.xaxis.get_ticklabels(), visible=False)

                    tp = mt.Tipper

                    txr = tp.mag_real * np.sin(tp.angle_real * np.pi / 180 + \
                                               np.pi * self.arrow_direction)
                    tyr = tp.mag_real * np.cos(tp.angle_real * np.pi / 180 + \
                                               np.pi * self.arrow_direction)

                    txi = tp.mag_imag * np.sin(tp.angle_imag * np.pi / 180 + \
                                               np.pi * self.arrow_direction)
                    tyi = tp.mag_imag * np.cos(tp.angle_imag * np.pi / 180 + \
                                               np.pi * self.arrow_direction)

                    nt = len(txr)

                    tiplist = []
                    tiplabel = []

                    for aa in range(nt):
                        xlenr = txr[aa] * mt.period[aa]
                        xleni = txi[aa] * mt.period[aa]

                        # --> plot real arrows
                        if self._plot_tipper.find('r') > 0:
                            axt.arrow(np.log10(mt.period[aa]),
                                      0,
                                      xlenr,
                                      tyr[aa],
                                      lw=self.arrow_lw,
                                      facecolor=self.arrow_color_real,
                                      edgecolor=self.arrow_color_real,
                                      head_width=self.arrow_head_width,
                                      head_length=self.arrow_head_length,
                                      length_includes_head=False)

                            if aa == 0:
                                line1 = axt.plot(0, 0, self.arrow_color_real)
                                tiplist.append(line1[0])
                                tiplabel.append('real')

                        # --> plot imaginary arrows
                        if self.plot_tipper.find('i') > 0:
                            axt.arrow(np.log10(mt.period[aa]),
                                      0,
                                      xleni,
                                      tyi[aa],
                                      facecolor=self.arrow_color_imag,
                                      edgecolor=self.arrow_color_imag,
                                      lw=self.arrow_lw,
                                      head_width=self.arrow_head_width,
                                      head_length=self.arrow_head_length,
                                      length_includes_head=False)
                            if aa == 0:
                                line2 = axt.plot(0, 0, self.arrow_color_imag)
                                tiplist.append(line2[0])
                                tiplabel.append('imag')

                    # make a line at 0 for reference
                    axt.plot(mt.period, [0] * nt, 'k', lw=.5)

                    if ii == 0:
                        axt.legend(tiplist, tiplabel,
                                   loc='upper left',
                                   markerscale=1,
                                   borderaxespad=.01,
                                   labelspacing=.07,
                                   handletextpad=.2,
                                   borderpad=.1,
                                   prop={'size': self.font_size})

                        axt.set_ylabel('Tipper', fontdict=fontdict)
                    else:
                        plt.setp(axt.get_yticklabels(), visible=False)

                    # set axis properties
                    axt.yaxis.set_major_locator(MultipleLocator(.2))
                    axt.yaxis.set_minor_locator(MultipleLocator(.1))
                    axt.set_xlabel('Period (s)', fontdict=fontdict)

                    axt.set_xscale('log', nonposx='clip')
                    if self.tipper_limits is None:
                        tmax = max([tyr.max(), tyi.max()])
                        if tmax > 1:
                            tmax = .899

                        tmin = min([tyr.min(), tyi.min()])
                        if tmin < -1:
                            tmin = -.899

                        self.tipper_limits = (tmin - .1, tmax + .1)

                    axt.set_ylim(self.tipper_limits)
                    axt.grid(True, alpha=.25, which='both',
                             color=(.25, .25, .25),
                             lw=.25)

                    tklabels = []
                    xticks = []
                    for tk in axt.get_xticks():
                        try:
                            tklabels.append(mtpl.labeldict[tk])
                            xticks.append(tk)
                        except KeyError:
                            pass
                    axt.set_xticks(xticks)
                    axt.set_xticklabels(tklabels,
                                        fontdict={'size': self.font_size})

                    if pdict['tip'] != nrows - 1:
                        plt.setp(axt.get_yticklabels(), visible=False)

                    # need to reset the xlimits caouse they get reset when calling
                    # set_ticks for some reason
                    axt.set_xlim(np.log10(self.xlimits[0]),
                                 np.log10(self.xlimits[1]))

                # ------plot strike angles----------------------------------------------
                if self._plot_strike.find('y') == 0:

                    if self._plot_strike.find('i') > 0:
                        # strike from invariants
                        zinv = Zinvariants(mt.Z)
                        s1 = zinv.strike

                        # fold angles so go from -90 to 90
                        s1[np.where(s1 > 90)] -= -180
                        s1[np.where(s1 < -90)] += 180

                        # plot strike with error bars
                        ps1 = axst.errorbar(mt.period,
                                            s1,
                                            marker=self.strike_inv_marker,
                                            ms=self.marker_size,
                                            mfc=self.strike_inv_color,
                                            mec=self.strike_inv_color,
                                            mew=self.marker_lw,
                                            ls='none',
                                            yerr=zinv.strike_err,
                                            ecolor=self.strike_inv_color,
                                            capsize=self.marker_size,
                                            elinewidth=self.marker_lw)

                        stlist.append(ps1[0])
                        stlabel.append('Z_inv')
                        st_maxlist.append(s1.max())
                        st_minlist.append(s1.min())

                    if self._plot_strike.find('p') > 0:
                        # strike from phase tensor
                        pt = mt.pt  # type: PhaseTensor
                        s2, s2_err = pt.azimuth, pt.azimuth_err

                        # fold angles to go from -90 to 90
                        s2[np.where(s2 > 90)] -= 180
                        s2[np.where(s2 < -90)] += 180

                        # plot strike with error bars
                        ps2 = axst.errorbar(mt.period,
                                            s2,
                                            marker=self.strike_pt_marker,
                                            ms=self.marker_size,
                                            mfc=self.strike_pt_color,
                                            mec=self.strike_pt_color,
                                            mew=self.marker_lw,
                                            ls='none',
                                            yerr=s2_err,
                                            ecolor=self.strike_pt_color,
                                            capsize=self.marker_size,
                                            elinewidth=self.marker_lw)

                        stlist.append(ps2[0])
                        stlabel.append('PT')
                        st_maxlist.append(s2.max())
                        st_minlist.append(s2.min())

                    if self._plot_strike.find('t') > 0:
                        # strike from tipper
                        tp = mt.Tipper
                        s3 = tp.angle_real + 90

                        # fold to go from -90 to 90
                        s3[np.where(s3 > 90)] -= 180
                        s3[np.where(s3 < -90)] += 180

                        # plot strike with error bars
                        ps3 = axst.errorbar(mt.period,
                                            s3,
                                            marker=self.strike_tip_marker,
                                            ms=self.marker_size,
                                            mfc=self.strike_tip_color,
                                            mec=self.strike_tip_color,
                                            mew=self.marker_lw,
                                            ls='none',
                                            yerr=np.zeros_like(s3),
                                            ecolor=self.strike_tip_color,
                                            capsize=self.marker_size,
                                            elinewidth=self.marker_lw)

                        stlist.append(ps3[0])
                        stlabel.append('Tip')
                        st_maxlist.append(s3.max())
                        st_minlist.append(s3.min())

                    # --> set axes properties
                    if self.strike_limits is None:
                        stmin = min(st_minlist)
                        if stmin - 3 < -90:
                            stmin -= 3
                        else:
                            stmin = -89.99

                        stmax = max(st_maxlist)
                        if stmin + 3 < 90:
                            stmin += 3
                        else:
                            stmin = 89.99
                        self.strike_limits = (-max([abs(stmin), abs(stmax)]),
                                              max([abs(stmin), abs(stmax)]))

                    axst.plot(axr.get_xlim(), [0, 0], color='k', lw=.5)
                    if ii == 0:
                        axst.set_ylabel('Strike',
                                        fontdict=fontdict)
                    else:
                        plt.setp(axst.get_yticklabels(), visible=False)
                    axst.set_xlabel('Period (s)',
                                    fontdict=fontdict)
                    axst.set_ylim(self.strike_limits)
                    axst.yaxis.set_major_locator(MultipleLocator(30))
                    axst.yaxis.set_minor_locator(MultipleLocator(5))
                    axst.set_xscale('log', nonposx='clip')
                    axst.grid(True, alpha=.25, which='both',
                              color=(.25, .25, .25),
                              lw=.25)
                    if ii == 0:
                        try:
                            axst.legend(stlist,
                                        stlabel,
                                        loc=3,
                                        markerscale=1,
                                        borderaxespad=.01,
                                        labelspacing=.07,
                                        handletextpad=.2,
                                        borderpad=.02,
                                        prop={'size': self.font_size - 1})
                        except:
                            pass

                    # set th xaxis tick labels to invisible
                    if pdict['strike'] != nrows - 1:
                        plt.setp(axst.xaxis.get_ticklabels(), visible=False)

                # ------plot skew angle---------------------------------------------
                if self._plot_skew == 'y':
                    # strike from phase tensor
                    pt = mt.pt
                    sk, sk_err = pt.beta, pt.beta_err

                    ps4 = axsk.errorbar(mt.period,
                                        sk,
                                        marker=self.skew_marker,
                                        ms=self.marker_size,
                                        mfc=self.skew_color,
                                        mec=self.skew_color,
                                        mew=self.marker_lw,
                                        ls='none',
                                        yerr=sk_err,
                                        ecolor=self.skew_color,
                                        capsize=self.marker_size,
                                        elinewidth=self.marker_lw)
                    stlist.append(ps4[0])
                    stlabel.append('Skew')
                    if self.skew_limits is None:
                        self.skew_limits = (-9, 9)

                    axsk.set_ylim(self.skew_limits)
                    axsk.yaxis.set_major_locator(MultipleLocator(3))
                    axsk.yaxis.set_minor_locator(MultipleLocator(1))
                    if ii ==0:
                        axsk.set_ylabel('Skew', fontdict)
                    else:
                        plt.setp(axsk.get_yticklabels(), visible=False)
                    axsk.set_xlabel('Period (s)', fontdict)
                    axsk.set_xscale('log', nonposx='clip')

                    # set th xaxis tick labels to invisible
                    if pdict['skew'] != nrows - 1:
                        plt.setp(axsk.xaxis.get_ticklabels(), visible=False)

                # ----plot phase tensor ellipse---------------------------------------
                if self._plot_pt == 'y':
                    # get phase tensor instance
                    pt = mt.pt

                    cmap = self.ellipse_cmap
                    ckmin = self.ellipse_range[0]
                    ckmax = self.ellipse_range[1]
                    try:
                        ckstep = float(self.ellipse_range[2])
                    except IndexError:
                        ckstep = 3

                    if cmap == 'mt_seg_bl2wh2rd':
                        bounds = np.arange(ckmin, ckmax + ckstep, ckstep)
                        nseg = float((ckmax - ckmin) / (2 * ckstep))

                    # get the properties to color the ellipses by
                    if self.ellipse_colorby == 'phiminang' or \
                            self.ellipse_colorby == 'phimin':
                        colorarray = pt.phimin

                    elif self.ellipse_colorby == 'phidet':
                        colorarray = np.sqrt(abs(pt.det)) * (180 / np.pi)

                    elif self.ellipse_colorby == 'skew' or \
                            self.ellipse_colorby == 'skew_seg':
                        colorarray = pt.beta

                    elif self.ellipse_colorby == 'ellipticity':
                        colorarray = pt.ellipticity

                    else:
                        raise NameError(self.ellipse_colorby + ' is not supported')

                    # -------------plot ellipses-----------------------------------
                    for kk, ff in enumerate(mt.period):
                        # make sure the ellipses will be visable
                        eheight = pt.phimin[kk] / pt.phimax[kk] * \
                                  self.ellipse_size
                        ewidth = pt.phimax[kk] / pt.phimax[kk] * \
                                 self.ellipse_size

                        # create an ellipse scaled by phimin and phimax and
                        # oriented along the azimuth which is calculated as
                        # clockwise but needs to be plotted counter-clockwise
                        # hence the negative sign.
                        ellipd = patches.Ellipse((np.log10(ff) * \
                                                  self.ellipse_spacing,
                                                  0),
                                                 width=ewidth,
                                                 height=eheight,
                                                 angle=90 - pt.azimuth[kk])

                        axpt.add_patch(ellipd)

                        # get ellipse color
                        if cmap.find('seg') > 0:
                            ellipd.set_facecolor(mtcl.get_plot_color(
                                colorarray[kk],
                                self.ellipse_colorby,
                                cmap,
                                ckmin,
                                ckmax,
                                bounds=bounds))
                        else:
                            ellipd.set_facecolor(mtcl.get_plot_color(
                                colorarray[kk],
                                self.ellipse_colorby,
                                cmap,
                                ckmin,
                                ckmax))

                    # ----set axes properties-----------------------------------------------
                    # --> set tick labels and limits
                    axpt.set_xlim(np.floor(np.log10(self.xlimits[0])),
                                  np.ceil(np.log10(self.xlimits[1])))

                    tklabels = []
                    xticks = []
                    for tk in axpt.get_xticks():
                        try:
                            tklabels.append(mtpl.labeldict[tk])
                            xticks.append(tk)
                        except KeyError:
                            pass
                    axpt.set_xticks(xticks)
                    axpt.set_xticklabels(tklabels,
                                         fontdict={'size': self.font_size})
                    axpt.set_xlabel('Period (s)', fontdict=fontdict)
                    axpt.set_ylim(ymin=-1.5 * self.ellipse_size,
                                  ymax=1.5 * self.ellipse_size)

                    axpt.grid(True,
                              alpha=.25,
                              which='major',
                              color=(.25, .25, .25),
                              lw=.25)

                    plt.setp(axpt.get_yticklabels(), visible=False)
                    if pdict['pt'] != nrows - 1:
                        plt.setp(axpt.get_xticklabels(), visible=False)

                    # add colorbar for PT only for first plot
                    if ii == 0:
                        axpos = axpt.get_position()
                        cb_position = (axpos.bounds[0] - .0575,
                                       axpos.bounds[1] + .02,
                                       .01,
                                       axpos.bounds[3] * .75)
                        cbax = self.fig.add_axes(cb_position)
                        if cmap == 'mt_seg_bl2wh2rd':
                            # make a color list
                            clist = [(cc, cc, 1)
                                     for cc in np.arange(0,
                                                         1 + 1. / (nseg),
                                                         1. / (nseg))] + \
                                    [(1, cc, cc)
                                     for cc in np.arange(1,
                                                         -1. / (nseg),
                                                         -1. / (nseg))]

                            # make segmented colormap
                            mt_seg_bl2wh2rd = colors.ListedColormap(clist)

                            # make bounds so that the middle is white
                            bounds = np.arange(ckmin - ckstep, ckmax + 2 * ckstep,
                                               ckstep)

                            # normalize the colors
                            norms = colors.BoundaryNorm(bounds,
                                                        mt_seg_bl2wh2rd.N)

                            # make the colorbar
                            cbpt = mcb.ColorbarBase(cbax,
                                                    cmap=mt_seg_bl2wh2rd,
                                                    norm=norms,
                                                    orientation='vertical',
                                                    ticks=bounds[1:-1])
                        else:
                            cbpt = mcb.ColorbarBase(cbax,
                                                    cmap=mtcl.cmapdict[cmap],
                                                    norm=colors.Normalize(vmin=ckmin,
                                                                          vmax=ckmax),
                                                    orientation='vertical')
                        cbpt.set_ticks([ckmin, (ckmax - ckmin) / 2, ckmax])
                        cbpt.set_ticklabels(['{0:.0f}'.format(ckmin),
                                             '{0:.0f}'.format((ckmax - ckmin) / 2),
                                             '{0:.0f}'.format(ckmax)])
                        cbpt.ax.yaxis.set_label_position('left')
                        cbpt.ax.yaxis.set_label_coords(-1.05, .5)
                        cbpt.ax.yaxis.tick_right()
                        cbpt.ax.tick_params(axis='y', direction='in')
                        cbpt.set_label(mtpl.ckdict[self.ellipse_colorby],
                                       fontdict={'size': self.font_size})

                # ==  == Plot the Z_xx, Z_yy components if desired ==
                if self.plot_num == 2:
                    # ---------plot the apparent resistivity----------------
                    axr2 = self.fig.add_subplot(gs[0, 1], sharex=axr, sharey=axr)
                    axr2.yaxis.set_label_coords(-.1, 0.5)

                    # res_xx
                    ebxxr = axr2.errorbar(mt.period,
                                          mt.Z.res_xx,
                                          marker=self.xy_marker,
                                          ms=self.marker_size,
                                          mfc=self.xy_mfc,
                                          mec=self.xy_color,
                                          mew=self.marker_lw,
                                          ls=self.xy_ls,
                                          yerr=mt.Z.res_err_xx,
                                          ecolor=self.xy_color,
                                          capsize=self.marker_size,
                                          elinewidth=self.marker_lw)

                    # res_yy
                    ebyyr = axr2.errorbar(mt.period,
                                          mt.Z.res_yy,
                                          marker=self.yx_marker,
                                          ms=self.marker_size,
                                          mfc=self.yx_mfc,
                                          mec=self.yx_color,
                                          mew=self.marker_lw,
                                          ls=self.yx_ls,
                                          yerr=mt.Z.res_err_yy,
                                          ecolor=self.yx_color,
                                          capsize=self.marker_size,
                                          elinewidth=self.marker_lw)

                    # --> set axes properties
                    plt.setp(axr2.get_xticklabels(), visible=False)
                    plt.setp(axr2.get_yticklabels(), visible=False)
                    axr2.set_yscale('log', nonposy='clip')
                    axr2.set_xscale('log', nonposx='clip')
                    axr2.set_xlim(self.x_limits)
                    axr2.set_ylim(self.res_limits)
                    axr2.grid(True,
                              alpha=.25,
                              which='both',
                              color=(.25, .25, .25),
                              lw=.25)
                    if ii == 0:
                        axr2.legend((ebxxr[0], ebyyr[0]),
                                    ('$Z_{xx}$', '$Z_{yy}$'),
                                    loc=3,
                                    markerscale=1,
                                    borderaxespad=.01,
                                    labelspacing=.07,
                                    handletextpad=.2,
                                    borderpad=.02)

                    # -----Plot the phase-----------------------------------
                    axp2 = self.fig.add_subplot(gs[1, 1], sharex=axr, sharey=axp)
                    axp2.yaxis.set_label_coords(-.1, 0.5)

                    # phase_xx
                    ebxxp = axp2.errorbar(mt.period,
                                          mt.Z.phase_xx,
                                          marker=self.xy_marker,
                                          ms=self.marker_size,
                                          mfc=self.xy_mfc,
                                          mec=self.xy_color,
                                          mew=self.marker_lw,
                                          ls=self.xy_ls,
                                          yerr=mt.Z.phase_err_xx,
                                          ecolor=self.xy_color,
                                          capsize=self.marker_size,
                                          elinewidth=self.marker_lw)

                    # phase_yy
                    ebyyp = axp2.errorbar(mt.period,
                                          mt.Z.phase_yy,
                                          marker=self.yx_marker,
                                          ms=self.marker_size,
                                          mfc=self.yx_mfc,
                                          mec=self.yx_color,
                                          mew=self.marker_lw,
                                          ls=self.yx_ls,
                                          yerr=mt.Z.phase_err_yy,
                                          ecolor=self.yx_color,
                                          capsize=self.marker_size,
                                          elinewidth=self.marker_lw)

                    # --> set axes properties
                    plt.setp(axp2.get_xticklabels(), visible=False)
                    plt.setp(axp2.get_yticklabels(), visible=False)
                    axp2.set_xlabel('Period (s)', fontdict)
                    axp2.set_xscale('log', nonposx='clip')
                    if self.phase_limits is None:
                        self.phase_limits=(-179.9,179.9)
                    axp2.set_ylim(self.phase_limits)
                    axp2.set_xlim(self.x_limits)
                    axp2.yaxis.set_major_locator(MultipleLocator(30))
                    axp2.yaxis.set_minor_locator(MultipleLocator(5))
                    # axp2.set_xticklabels(tklabels,
                    #                      fontdict={'size': self.font_size})
                    axp2.grid(True,
                              alpha=.25,
                              which='both',
                              color=(.25, .25, .25),
                              lw=.25)

                    if len(list(pdict.keys())) > 2:
                        plt.setp(axp2.xaxis.get_ticklabels(), visible=False)
                        plt.setp(axp2.xaxis.get_label(), visible=False)

                # == =Plot the Determinant if desired ==  ==  ==  ==
                if self.plot_num == 3:

                    # res_det
                    ebdetr = axr.errorbar(mt.period,
                                          rp.res_det,
                                          marker=self.det_marker,
                                          ms=self.marker_size,
                                          mfc=self.det_mfc,
                                          mec=self.det_color,
                                          mew=self.marker_lw,
                                          ls=self.det_ls,
                                          yerr=rp.res_det_err,
                                          ecolor=self.det_color,
                                          capsize=self.marker_size,
                                          elinewidth=self.marker_lw)

                    # phase_det
                    ebdetp = axp.errorbar(mt.period,
                                          rp.phase_det,
                                          marker=self.det_marker,
                                          ms=self.marker_size,
                                          mfc=self.det_mfc,
                                          mec=self.det_color,
                                          mew=self.marker_lw,
                                          ls=self.det_ls,
                                          yerr=rp.phase_det_err,
                                          ecolor=self.det_color,
                                          capsize=self.marker_size,
                                          elinewidth=self.marker_lw)

                    # --> set axes properties
                    plt.setp(axr.get_xticklabels(), visible=False)
                    if ii == 0:
                        axr.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                       fontdict=fontdict)
                    else:
                        plt.setp(axr.get_yticklabels(), visible=False)

                    axr.set_yscale('log', nonposy='clip')
                    axr.set_xscale('log', nonposx='clip')
                    axr.set_ylim(self.res_limits)
                    axr.set_xlim(self.xlimits)
                    axr.grid(True, alpha=.25,
                             which='both',
                             color=(.25, .25, .25),
                             lw=.25)

                    # --> set axes properties
                    axp.set_xlabel('Period (s)', fontdict)

                    if ii == 0:
                        axp.set_ylabel('Phase (deg)', fontdict)

                    else:
                        plt.setp(axp.get_yticklabels(), visible=False)

                    axp.set_xscale('log', nonposx='clip')
                    axp.set_ylim(self.phase_limits)
                    axp.yaxis.set_major_locator(MultipleLocator(15))
                    axp.yaxis.set_minor_locator(MultipleLocator(5))
                    tklabels = [mtpl.labeldict[tt]
                                for tt in np.arange(np.log10(self.xlimits[0]),
                                                    np.log10(self.xlimits[1]) + 1)]
                    tklabels[0] = ''
                    tklabels[-1] = ''

                    axp.set_xticklabels(tklabels,
                                        fontdict={'size': self.font_size})
                    axp.grid(True, alpha=.25,
                             which='both',
                             color=(.25, .25, .25),
                             lw=.25)

                # make title and show

                axr.set_title(mt.station, fontsize=self.font_size,
                              fontweight='bold')
            if show:
                plt.show()

        # ===Plot all responses into one plot to compare changes ==
        if self.plot_style == 'compare':
            ns = len(self.mt_list)

            # make color lists for the plots going light to dark
            cxy = [(0, 0 + float(cc) / ns, 1 - float(cc) / ns) for cc in range(ns)]
            cyx = [(1, float(cc) / ns, 0) for cc in range(ns)]
            cdet = [(0, 1 - float(cc) / ns, 0) for cc in range(ns)]
            ctipr = [(.75 * cc / ns, .75 * cc / ns, .75 * cc / ns) for cc in range(ns)]
            ctipi = [(float(cc) / ns, 1 - float(cc) / ns, .25) for cc in range(ns)]
            cst = [(.5 * cc / ns, 0, .5 * cc / ns) for cc in range(ns)]

            # make marker lists for the different components
            mxy = ['s', 'D', 'x', '+', '*', '1', '3', '4'] * 5
            myx = ['o', 'h', '8', 'p', 'H', 7, 4, 6] * 5

            legendlistxy = []
            legendlistyx = []
            stationlist = []
            tiplist = []
            stlist = []
            sklist = []

            # set some parameters of the figure and subplot spacing
            plt.rcParams['font.size'] = self.font_size
            plt.rcParams['figure.subplot.bottom'] = .1
            plt.rcParams['figure.subplot.top'] = .97
            plt.rcParams['figure.subplot.left'] = .08
            plt.rcParams['figure.subplot.right'] = .98

            # set the font properties for the axis labels
            fontdict = {'size': self.font_size + 1, 'weight': 'bold'}

            # set figure size according to what the plot will be.
            if self.fig_size is None:
                if self.plot_num == 1 or self.plot_num == 3:
                    self.fig_size = [8, 6]
                    pass

                elif self.plot_num == 2:
                    self.fig_size = [8, 6]
                    nrows += 1

            # make a figure instance
            self.fig = plt.figure(self.fig_num, self.fig_size,
                                  dpi=self.fig_dpi)

            # make a grid as usual, but put xy and yx in different plots
            # otherwise the plot is too busy to see what's going on.
            hr = [2, 1.5] + [1] * (nrows - 2)
            gs = gridspec.GridSpec(nrows, 2, height_ratios=hr, hspace=.05)

            # --> make figure for xy,yx components
            if self.plot_num == 1 or self.plot_num == 3:
                # set label coordinates
                labelcoords = (-0.125, 0.5)

                # space out the subplots
                gs.update(hspace=.05, wspace=.02, left=.1)

            # --> make figure for all 4 components
            elif self.plot_num == 2:
                # set label coordinates
                labelcoords = (-0.125, 0.5)

                # space out the subplots
                gs.update(hspace=.05, wspace=.02, left=.07)

                for key in pdict:
                    if key != 'res' and key != 'phase':
                        pdict[key] += 1

            # --> create the axes instances
            # apparent resistivity axis
            self.axrxy = self.fig.add_subplot(gs[0, 0])
            self.axryx = self.fig.add_subplot(gs[0, 1], sharex=self.axrxy,
                                              sharey=self.axrxy)

            # phase axis that shares period axis with resistivity
            self.axpxy = self.fig.add_subplot(gs[1, 0], sharex=self.axrxy)
            self.axpyx = self.fig.add_subplot(gs[1, 1], sharex=self.axrxy,
                                              sharey=self.axpxy)

            # place y coordinate labels in the same location
            # self.axrxy.yaxis.set_label_coords(labelcoords[0], labelcoords[1])
            # self.axpxy.yaxis.set_label_coords(labelcoords[0], labelcoords[1])

            # --> plot tipper
            try:
                self.axt = self.fig.add_subplot(gs[pdict['tip'], :])
                # self.axt.yaxis.set_label_coords(labelcoords[0] * .5,
                #                                 labelcoords[1])
            except KeyError:
                pass

            # --> plot phase tensors
            try:
                # can't share axis because not on the same scale
                self.axpt = self.fig.add_subplot(gs[pdict['pt'], :],
                                                 aspect='equal')
                # self.axpt.yaxis.set_label_coords(labelcoords[0] * .5,
                #                                  labelcoords[1])
            except KeyError:
                pass

            # --> plot strike
            try:
                self.axst = self.fig.add_subplot(gs[pdict['strike'], :],
                                                 sharex=self.axrxy)
                # self.axst.yaxis.set_label_coords(labelcoords[0] * .5,
                #                                  labelcoords[1])
            except KeyError:
                pass

            # --> plot skew
            try:
                self.axsk = self.fig.add_subplot(gs[pdict['skew'], :],
                                                 sharex=self.axrxy)
                # self.axsk.yaxis.set_label_coords(labelcoords[0] * .5,
                #                                  labelcoords[1])
            except KeyError:
                pass

            for ii, mt in enumerate(self.mt_list):
                # get the reistivity and phase object

                # set x-axis limits from short period to long period
                if self.xlimits is None:
                    self.xlimits = (10 ** (np.floor(np.log10(mt.period.min()))),
                                    10 ** (np.ceil(np.log10(mt.period.max()))))
                else:
                    self.xlimits = (10 ** min([np.floor(np.log10(self.xlimits[0])),
                                               np.floor(np.log10(mt.period.min()))]),
                                    10 ** max([np.ceil(np.log10(self.xlimits[1])),
                                               np.ceil(np.log10(mt.period.max()))]))
                if self.phase_limits is None:
                    self.phase_limits = (0, 89.9)

                stationlist.append(mt.station)

                # ==  ==  ==  == =Plot Z_xy and Z_yx ==
                # if self.plot_num == 1 or self.plot_num == 2:
                # ---------plot the apparent resistivity--------------------
                # --> plot as error bars and just as points xy-blue, yx-red
                # res_xy
                ebxyr = self.axrxy.errorbar(mt.period,
                                            mt.Z.res_xy,
                                            color=cxy[ii],
                                            marker=mxy[ii % len(mxy)],
                                            ms=self.marker_size,
                                            mfc='None',
                                            mec=cxy[ii],
                                            mew=self.marker_lw,
                                            ls=self.xy_ls,
                                            yerr=mt.Z.res_err_xy,
                                            ecolor=cxy[ii],
                                            capsize=self.marker_size,
                                            elinewidth=self.marker_lw)

                # res_yx
                ebyxr = self.axryx.errorbar(mt.period,
                                            mt.Z.res_yx,
                                            color=cyx[ii],
                                            marker=myx[ii % len(myx)],
                                            ms=self.marker_size,
                                            mfc='None',
                                            mec=cyx[ii],
                                            mew=self.marker_lw,
                                            ls=self.yx_ls,
                                            yerr=mt.Z.res_err_yx,
                                            ecolor=cyx[ii],
                                            capsize=self.marker_size,
                                            elinewidth=self.marker_lw)

                # -----Plot the phase---------------------------------------
                # phase_xy
                self.axpxy.errorbar(mt.period,
                                    mt.Z.phase_xy,
                                    color=cxy[ii],
                                    marker=mxy[ii % len(mxy)],
                                    ms=self.marker_size,
                                    mfc='None',
                                    mec=cxy[ii],
                                    mew=self.marker_lw,
                                    ls=self.xy_ls,
                                    yerr=mt.Z.phase_err_xy,
                                    ecolor=cxy[ii],
                                    capsize=self.marker_size,
                                    elinewidth=self.marker_lw)

                # phase_yx: Note add 180 to place it in same quadrant as
                # phase_xy
                self.axpyx.errorbar(mt.period,
                                    mt.Z.phase_yx + 180,
                                    color=cyx[ii],
                                    marker=myx[ii % len(myx)],
                                    ms=self.marker_size,
                                    mfc='None',
                                    mec=cyx[ii],
                                    mew=self.marker_lw,
                                    ls=self.yx_ls,
                                    yerr=mt.Z.phase_err_yx,
                                    ecolor=cyx[ii],
                                    capsize=self.marker_size,
                                    elinewidth=self.marker_lw)

                legendlistxy.append(ebxyr)
                legendlistyx.append(ebyxr)

                # ==== Plot the Z_xx, Z_yy components if desired ==
                if self.plot_num == 2:
                    # ---------plot the apparent resistivity----------------
                    self.axr2xx = self.fig.add_subplot(gs[2, 0],
                                                       sharex=self.axrxy)
                    self.axr2xx.yaxis.set_label_coords(-.095, 0.5)
                    self.axr2yy = self.fig.add_subplot(gs[2, 1],
                                                       sharex=self.axrxy)

                    # res_xx
                    ebxxr = self.axr2xx.errorbar(mt.period,
                                                 mt.Z.res_xx,
                                                 color=cxy[ii],
                                                 marker=mxy[ii % len(mxy)],
                                                 ms=self.marker_size,
                                                 mfc='None',
                                                 mec=cxy[ii],
                                                 mew=self.marker_lw,
                                                 ls=self.xy_ls,
                                                 yerr=mt.Z.res_err_xx,
                                                 ecolor=cxy[ii],
                                                 capsize=self.marker_size,
                                                 elinewidth=self.marker_lw)

                    # res_yy
                    ebyyr = self.axr2yy.errorbar(mt.period,
                                                 mt.Z.res_yy,
                                                 color=cyx[ii],
                                                 marker=myx[ii % len(myx)],
                                                 ms=self.marker_size,
                                                 mfc='None',
                                                 mec=cyx[ii],
                                                 mew=self.marker_lw,
                                                 ls=self.yx_ls,
                                                 yerr=mt.Z.res_err_yy,
                                                 ecolor=cyx[ii],
                                                 capsize=self.marker_size,
                                                 elinewidth=self.marker_lw)

                    # -----Plot the phase-----------------------------------
                    self.axp2xx = self.fig.add_subplot(gs[2, 0],
                                                       sharex=self.axrxy)
                    self.axp2xx.yaxis.set_label_coords(-.095, 0.5)
                    self.axp2yy = self.fig.add_subplot(gs[2, 1],
                                                       sharex=self.axrxy)

                    # phase_xx
                    ebxxp = self.axp2xx.errorbar(mt.period,
                                                 mt.Z.phase_xx,
                                                 color=cxy[ii],
                                                 marker=mxy[ii % len(mxy)],
                                                 ms=self.marker_size,
                                                 mfc='None',
                                                 mec=cxy[ii],
                                                 mew=self.marker_lw,
                                                 ls=self.xy_ls,
                                                 yerr=mt.Z.phase_err_xx,
                                                 ecolor=cxy[ii],
                                                 capsize=self.marker_size,
                                                 elinewidth=self.marker_lw)

                    # phase_yy
                    ebyyp = self.axp2yy.errorbar(mt.period,
                                                 mt.Z.phase_yy,
                                                 color=cyx[ii],
                                                 marker=myx[ii % len(mxy)],
                                                 ms=self.marker_size,
                                                 mfc='None',
                                                 mec=cyx[ii],
                                                 mew=self.marker_lw,
                                                 ls=self.yx_ls,
                                                 yerr=mt.Z.phase_err_yy,
                                                 ecolor=cyx[ii],
                                                 capsize=self.marker_size,
                                                 elinewidth=self.marker_lw)

                # ===Plot the Determinant if desired ==
                if self.plot_num == 3:
                    # res_det
                    ebdetr = self.axrxy.errorbar(mt.period,
                                                 mt.Z.res_det,
                                                 color=cxy[ii],
                                                 marker=mxy[ii % len(mxy)],
                                                 ms=self.marker_size,
                                                 mfc='None',
                                                 mec=cdet[ii],
                                                 mew=self.marker_lw,
                                                 ls=self.det_ls,
                                                 yerr=mt.Z.res_det_err,
                                                 ecolor=cdet[ii],
                                                 capsize=self.marker_size,
                                                 elinewidth=self.marker_lw)

                    # phase_det
                    ebdetp = self.axpxy.errorbar(mt.period,
                                                 mt.Z.phase_det,
                                                 color=cyx[ii],
                                                 marker=mxy[ii % len(mxy)],
                                                 ms=self.marker_size,
                                                 mfc='None',
                                                 mec=cdet[ii],
                                                 mew=self.marker_lw,
                                                 ls=self.det_ls,
                                                 yerr=mt.Z.phase_det_err,
                                                 ecolor=cdet[ii],
                                                 capsize=self.marker_size,
                                                 elinewidth=self.marker_lw)

                    legendlistxy.append(ebdetr)

                # -----plot tipper----------------------------------------------
                if self._plot_tipper.find('y') == 0:

                    txr = mt.Tipper.mag_real * np.sin(mt.Tipper.angle_real * np.pi / 180 + \
                                                      np.pi * self.arrow_direction)
                    tyr = mt.Tipper.mag_real * np.cos(mt.Tipper.angle_real * np.pi / 180 + \
                                                      np.pi * self.arrow_direction)

                    txi = mt.Tipper.mag_imag * np.sin(mt.Tipper.angle_imag * np.pi / 180 + \
                                                      np.pi * self.arrow_direction)
                    tyi = mt.Tipper.mag_imag * np.cos(mt.Tipper.angle_imag * np.pi / 180 + \
                                                      np.pi * self.arrow_direction)

                    nt = len(txr)

                    for aa in range(nt):
                        xlenr = txr[aa] * np.log10(mt.period[aa])
                        xleni = txi[aa] * np.log10(mt.period[aa])

                        if self.tipper_limits is None:
                            tmax = max([tyr.max(), tyi.max()])
                            tmin = min([tyr.min(), tyi.min()])
                            if np.isnan(tmax):
                                tmax = 1.0
                            if np.isnan(tmin):
                                tmin = -1.0
                            tmin = max([-1, tmin])
                            tmax = min([1, tmax])
                            self.tipper_limits = (tmin - .1, tmax + .1)
                        else:
                            tmax = max([tyr.max(), tyi.max(), self.tipper_limits[1] - .1]) + .1
                            tmin = min([tyr.min(), tyi.min(), self.tipper_limits[0] + .1]) - .1
                            if np.isnan(tmax):
                                tmax = 1.0
                            if np.isnan(tmin):
                                tmin = -1.0
                            tmin = max([-1, tmin])
                            tmax = min([1, tmax])
                            self.tipper_limits = (tmin, tmax)

                        # --> plot real arrows
                        if self._plot_tipper.find('r') > 0:
                            self.axt.arrow(np.log10(mt.period[aa]),
                                           0,
                                           xlenr,
                                           tyr[aa],
                                           lw=self.arrow_lw,
                                           facecolor=ctipr[ii],
                                           edgecolor=ctipr[ii],
                                           head_width=self.arrow_head_width,
                                           head_length=self.arrow_head_length,
                                           length_includes_head=False)

                        # --> plot imaginary arrows
                        if self._plot_tipper.find('i') > 0:
                            self.axt.arrow(np.log10(mt.period[aa]),
                                           0,
                                           xleni,
                                           tyi[aa],
                                           lw=self.arrow_lw,
                                           head_width=self.arrow_head_width,
                                           head_length=self.arrow_head_length,
                                           length_includes_head=False)

                    lt = self.axt.plot(0, 0, lw=1, color=ctipr[ii])
                    tiplist.append(lt[0])

                # ------plot strike angles----------------------------------------------
                if self._plot_strike.find('y') == 0:

                    #                    if self._plot_strike.find('i') > 0:
                    #                        #strike from invariants
                    #                        zinv = mt.Z.invariants
                    #                        s1 = zinv.strike
                    #
                    #                        #fold angles so go from -90 to 90
                    #                        s1[np.where(s1>90)] -= -180
                    #                        s1[np.where(s1<-90)] += 180
                    #
                    #                        #plot strike with error bars
                    #                        ps1 = self.axst.errorbar(mt.period,
                    #                                                s1,
                    #                                                marker=mxy[ii % len(mxy)],
                    #                                                ms=self.marker_size,
                    #                                                mfc=cst[ii],
                    #                                                mec=cst[ii],
                    #                                                mew=self.marker_lw,
                    #                                                ls='none',
                    #                                                yerr=zinv.strike_err,
                    #                                                ecolor=cst[ii],
                    #                                                capsize=self.marker_size,
                    #                                                elinewidth=self.marker_lw)
                    #
                    #                        stlist.append(ps1[0])

                    if self._plot_strike.find('p') > 0:
                        # strike from phase tensor
                        s2 = mt.pt.azimuth
                        s2_err = mt.pt.azimuth_err

                        # fold angles to go from -90 to 90
                        s2[np.where(s2 > 90)] -= 180
                        s2[np.where(s2 < -90)] += 180

                        # plot strike with error bars
                        ps2 = self.axst.errorbar(mt.period,
                                                 s2,
                                                 marker=myx[ii % len(myx)],
                                                 ms=self.marker_size,
                                                 mfc=cxy[ii],
                                                 mec=cxy[ii],
                                                 mew=self.marker_lw,
                                                 ls='none',
                                                 yerr=s2_err,
                                                 ecolor=cxy[ii],
                                                 capsize=self.marker_size,
                                                 elinewidth=self.marker_lw)

                        stlist.append(ps2[0])

                    if self._plot_strike.find('t') > 0:
                        # strike from tipper
                        s3 = mt.Tipper.angle_real + 90

                        # fold to go from -90 to 90
                        s3[np.where(s3 > 90)] -= 180
                        s3[np.where(s3 < -90)] += 180

                        # plot strike with error bars
                        ps3 = self.axst.errorbar(mt.period,
                                                 s3,
                                                 marker=mxy[ii % len(mxy)],
                                                 ms=self.marker_size,
                                                 mfc=ctipr[ii],
                                                 mec=ctipr[ii],
                                                 mew=self.marker_lw,
                                                 ls='none',
                                                 yerr=np.zeros_like(s3),
                                                 ecolor=ctipr[ii],
                                                 capsize=self.marker_size,
                                                 elinewidth=self.marker_lw)

                        stlist.append(ps3[0])

                # ------plot skew angle---------------------------------------------
                if self._plot_skew == 'y':
                    # strike from phase tensor
                    sk = mt.pt.beta
                    sk_err = mt.pt.beta_err

                    ps4 = self.axsk.errorbar(mt.period,
                                             sk,
                                             marker=mxy[ii % len(mxy)],
                                             ms=self.marker_size,
                                             mfc=cxy[ii],
                                             mec=cxy[ii],
                                             mew=self.marker_lw,
                                             ls='none',
                                             yerr=sk_err,
                                             ecolor=cxy[ii],
                                             capsize=self.marker_size,
                                             elinewidth=self.marker_lw)
                    stlist.append(ps4[0])

                # ----plot phase tensor ellipse---------------------------------------
                if self._plot_pt == 'y':
                    # get phase tensor instance
                    pt = mt.pt

                    cmap = self.ellipse_cmap
                    ckmin = self.ellipse_range[0]
                    ckmax = self.ellipse_range[1]
                    try:
                        ckstep = float(self.ellipse_range[2])
                    except IndexError:
                        ckstep = 3

                    if cmap == 'mt_seg_bl2wh2rd':
                        bounds = np.arange(ckmin, ckmax + ckstep, ckstep)
                        nseg = float((ckmax - ckmin) / (2 * ckstep))

                    # get the properties to color the ellipses by
                    if self.ellipse_colorby == 'phiminang' or \
                            self.ellipse_colorby == 'phimin':
                        colorarray = mt.pt.phimin

                    elif self.ellipse_colorby == 'phidet':
                        colorarray = np.sqrt(abs(mt.pt.det)) * (180 / np.pi)

                    elif self.ellipse_colorby == 'skew' or \
                            self.ellipse_colorby == 'skew_seg':
                        colorarray = mt.pt.beta

                    elif self.ellipse_colorby == 'ellipticity':
                        colorarray = mt.pt.ellipticity

                    else:
                        raise NameError(self.ellipse_colorby + ' is not supported')

                    # -------------plot ellipses-----------------------------------
                    for kk, ff in enumerate(mt.period):
                        # make sure the ellipses will be visable
                        eheight = mt.pt.phimin[kk] / mt.pt.phimax[kk] * \
                                  self.ellipse_size
                        ewidth = mt.pt.phimax[kk] / mt.pt.phimax[kk] * \
                                 self.ellipse_size

                        # create an ellipse scaled by phimin and phimax and oriented
                        # along the azimuth which is calculated as clockwise but needs
                        # to be plotted counter-clockwise hence the negative sign.
                        ellipd = patches.Ellipse((np.log10(ff) * self.ellipse_spacing,
                                                  ii * self.ellipse_size * 1.5),
                                                 width=ewidth,
                                                 height=eheight,
                                                 angle=90 - pt.azimuth[kk])

                        self.axpt.add_patch(ellipd)

                        # get ellipse color
                        if cmap.find('seg') > 0:
                            ellipd.set_facecolor(mtcl.get_plot_color(colorarray[kk],
                                                                     self.ellipse_colorby,
                                                                     cmap,
                                                                     ckmin,
                                                                     ckmax,
                                                                     bounds=bounds))
                        else:
                            ellipd.set_facecolor(mtcl.get_plot_color(colorarray[kk],
                                                                     self.ellipse_colorby,
                                                                     cmap,
                                                                     ckmin,
                                                                     ckmax))
                        ellipd.set_edgecolor(cxy[ii])

            # -------set axis properties---------------------------------------
            self.axrxy.set_yscale('log', nonposy='clip')
            self.axrxy.set_xscale('log', nonposx='clip')
            self.axrxy.set_ylim(self.res_limits)
            self.axrxy.set_xlim(self.x_limits)
            self.axrxy.grid(True, alpha=.25,
                            which='both',
                            color=(.25, .25, .25),
                            lw=.25)

            # make a text label in upper left hand corner
            # label the plot with a text box
            if self.text_location is None:
                txloc = self.xlimits[0] * self.text_xpad
                tyloc = self.axrxy.get_ylim()[1] * self.text_ypad
            else:
                txloc = self.text_location[0]
                tyloc = self.text_location[1]

            self.text = self.axrxy.text(txloc,
                                        tyloc,
                                        '$Z_{xy}$',
                                        fontdict={'size': self.text_size,
                                                  'weight': self.text_weight},
                                        verticalalignment='top',
                                        horizontalalignment='left',
                                        bbox={'facecolor': 'white', 'alpha': 1})

            plt.setp(self.axrxy.get_xticklabels(), visible=False)

            self.axrxy.set_ylabel('App. Resistivity($\Omega \cdot$m)',
                                  fontdict=fontdict)

            self.axryx.set_yscale('log', nonposy='clip')
            self.axryx.set_xscale('log', nonposx='clip')
            self.axryx.set_ylim(self.res_limits)
            self.axryx.set_xlim(self.x_limits)
            self.axryx.grid(True, alpha=.25,
                            which='both',
                            color=(.25, .25, .25),
                            lw=.25)

            self.text = self.axryx.text(txloc,
                                        tyloc,
                                        '$Z_{yx}$',
                                        fontdict={'size': self.text_size,
                                                  'weight': self.text_weight},
                                        verticalalignment='top',
                                        horizontalalignment='left',
                                        bbox={'facecolor': 'white', 'alpha': 1})

            plt.setp(self.axryx.get_xticklabels(), visible=False)
            plt.setp(self.axryx.get_yticklabels(), visible=False)

            # check the phase to see if any point are outside of [0:90]
            if self.phase_limits is None:
                self.phase_limits = (0, 89.99)

            # --> set axes properties
            self.axpxy.set_xlabel('Period(s)', fontdict=fontdict)
            self.axpxy.set_ylabel('Phase(deg)', fontdict=fontdict)
            self.axpxy.set_xscale('log', nonposx='clip')
            self.axpxy.set_ylim(self.phase_limits)
            self.axpxy.yaxis.set_major_locator(MultipleLocator(15))
            self.axpxy.yaxis.set_minor_locator(MultipleLocator(5))
            self.axpxy.grid(True, alpha=.25,
                            which='both',
                            color=(.25, .25, .25),
                            lw=.25)
            if len(list(pdict.keys())) > 2:
                plt.setp(self.axpxy.xaxis.get_ticklabels(), visible=False)
                self.axpxy.set_xlabel('')

            self.axpyx.set_xlabel('Period(s)', fontdict=fontdict)
            self.axpyx.set_xscale('log', nonposx='clip')
            self.axpyx.set_ylim(self.phase_limits)
            self.axpyx.yaxis.set_major_locator(MultipleLocator(15))
            self.axpyx.yaxis.set_minor_locator(MultipleLocator(5))
            self.axpyx.grid(True, alpha=.25,
                            which='both',
                            color=(.25, .25, .25),
                            lw=.25)
            plt.setp(self.axpyx.yaxis.get_ticklabels(), visible=False)

            if len(list(pdict.keys())) > 2:
                plt.setp(self.axpyx.xaxis.get_ticklabels(), visible=False)
                self.axpyx.set_xlabel('')

            # make legend
            if self.plot_num == 1:
                self.axrxy.legend(legendlistxy,
                                  stationlist,
                                  loc=3,
                                  ncol=2,
                                  markerscale=.75,
                                  borderaxespad=.01,
                                  labelspacing=.07,
                                  handletextpad=.2,
                                  borderpad=.25)

                self.axryx.legend(legendlistyx,
                                  stationlist,
                                  loc=3,
                                  ncol=2,
                                  markerscale=.75,
                                  borderaxespad=.01,
                                  labelspacing=.07,
                                  handletextpad=.2,
                                  borderpad=.25)

            elif self.plot_num == 3:
                llist = [ll[0] for ll in legendlistxy]
                slist = [ss + '_det' for ss in stationlist]

                self.axrxy.legend(llist,
                                  slist,
                                  loc=3,
                                  markerscale=.75,
                                  borderaxespad=.01,
                                  labelspacing=.07,
                                  handletextpad=.2,
                                  borderpad=.25)
                self.axryx.legend(llist,
                                  slist,
                                  loc=3,
                                  markerscale=.75,
                                  borderaxespad=.01,
                                  labelspacing=.07,
                                  handletextpad=.2,
                                  borderpad=.25)

            if self.plot_num == 2:
                # --> set axes properties for resxx
                self.axrxy.set_yscale('log', nonposy='clip')
                self.axrxy.set_xscale('log', nonposx='clip')
                self.axrxy.set_xlim(self.x_limits)
                self.axrxy.grid(True,
                                alpha=.25,
                                which='both',
                                color=(.25, .25, .25),
                                lw=.25)
                plt.setp(self.axrxy.get_xticklabels(), visible=False)

                # --> set axes properties for resyy
                self.axryx.set_yscale('log', nonposy='clip')
                self.axryx.set_xscale('log', nonposx='clip')
                self.axryx.set_xlim(self.x_limits)
                self.axryx.grid(True,
                                alpha=.25,
                                which='both',
                                color=(.25, .25, .25),
                                lw=.25)
                plt.setp(self.axryx.get_xticklabels(), visible=False)

                # --> set axes properties Phasexx
                self.axpxy.set_xlabel('Period(s)', fontdict)
                self.axpxy.set_xscale('log', nonposx='clip')
                self.axpxy.set_ylim(ymin=-179.9, ymax=179.9)
                self.axpxy.yaxis.set_major_locator(MultipleLocator(30))
                self.axpxy.yaxis.set_minor_locator(MultipleLocator(5))
                self.axpxy.grid(True,
                                alpha=.25,
                                which='both',
                                color=(.25, .25, .25),
                                lw=.25)

                # --> set axes properties Phaseyy
                self.axpyx.set_xlabel('Period(s)', fontdict)
                self.axpyx.set_xscale('log', nonposx='clip')
                self.axpyx.set_ylim(ymin=-179.9, ymax=179.9)
                self.axpyx.yaxis.set_major_locator(MultipleLocator(30))
                self.axpyx.yaxis.set_minor_locator(MultipleLocator(5))
                self.axpyx.grid(True,
                                alpha=.25,
                                which='both',
                                color=(.25, .25, .25),
                                lw=.25)
                if len(list(pdict.keys())) > 3:
                    plt.setp(self.axpxy.xaxis.get_ticklabels(), visible=False)
                    self.axpxy.set_xlabel('')
                    plt.setp(self.axpyx.xaxis.get_ticklabels(), visible=False)
                    self.axpyx.set_xlabel('')

            if self._plot_tipper.find('y') == 0:
                self.axt.plot(self.axt.get_xlim(), [0, 0], color='k', lw=.5)
                # --> set axis properties Tipper
                if self.plot_num == 2:
                    plt.setp(self.axpxy.get_xticklabels(), visible=False)
                    self.axpxy.set_xlabel('')
                    plt.setp(self.axpyx.get_xticklabels(), visible=False)
                    self.axpyx.set_xlabel('')

                self.axt.yaxis.set_major_locator(MultipleLocator(.2))
                self.axt.yaxis.set_minor_locator(MultipleLocator(.1))
                self.axt.set_xlabel('Period(s)', fontdict=fontdict)
                self.axt.set_ylabel('Tipper', fontdict=fontdict)
                self.axt.set_xlim(np.log10(self.xlimits[0]),
                                  np.log10(self.xlimits[1]))
                tklabels = []
                xticks = []
                for tk in self.axt.get_xticks():
                    try:
                        tklabels.append(mtpl.labeldict[tk])
                        xticks.append(tk)
                    except KeyError:
                        pass
                self.axt.set_xticks(xticks)
                self.axt.set_xticklabels(tklabels,
                                         fontdict={'size': self.font_size})

                self.axt.set_ylim(self.tipper_limits)
                self.axt.grid(True, alpha=.25,
                              which='both',
                              color=(.25, .25, .25),
                              lw=.25)

                self.axt.legend(tiplist,
                                stationlist,
                                loc=3,
                                ncol=2,
                                markerscale=1,
                                borderaxespad=.01,
                                labelspacing=.07,
                                handletextpad=.2,
                                borderpad=.02)

                # need to reset the xlimits caouse they get reset when calling
                # set_ticks for some reason
                self.axt.set_xlim(np.log10(self.xlimits[0]),
                                  np.log10(self.xlimits[1]))

                if pdict['tip'] != nrows - 1:
                    plt.setp(self.axt.xaxis.get_ticklabels(), visible=False)
                    self.axt.set_xlabel(' ')

            # --> set axes properties for strike and skew
            if self._plot_strike[0] == 'y':

                if self.strike_limits is None:
                    self.strike_limits = (-89.99, 89.99)
                self.axst.plot(self.axrxy.get_xlim(), [0, 0], color='k', lw=.5)

                self.axst.set_ylabel('Strike(deg)',
                                     fontdict=fontdict)
                self.axst.set_xlabel('Period(s)',
                                     fontdict=fontdict)
                self.axst.set_ylim(self.strike_limits)
                self.axst.yaxis.set_major_locator(MultipleLocator(30))
                self.axst.yaxis.set_minor_locator(MultipleLocator(5))
                self.axst.set_xscale('log', nonposx='clip')
                self.axst.grid(True, alpha=.25,
                               which='both',
                               color=(.25, .25, .25),
                               lw=.25)
                # self.axst.legend(stlist,
                #                stationlist,
                #                loc=3,
                #                ncol=2,
                #                markerscale=1,
                #                borderaxespad=.01,
                #                labelspacing=.07,
                #                handletextpad=.2,
                #                borderpad=.02)
                if pdict['strike'] != nrows - 1:
                    plt.setp(self.axst.xaxis.get_ticklabels(), visible=False)
                    self.axst.set_xlabel(' ')

            # --> set axes properties for skew
            if self._plot_skew == 'y':
                self.axsk.set_ylim(self.skew_limits)
                self.axsk.yaxis.set_major_locator(MultipleLocator(3))
                self.axsk.yaxis.set_minor_locator(MultipleLocator(1))
                self.axsk.set_ylabel('Skew(deg)', fontdict)
                self.axsk.set_xlabel('Period(s)',
                                     fontdict=fontdict)
                self.axsk.set_xscale('log', nonposx='clip')
                self.axsk.grid(True, alpha=.25,
                               which='both',
                               color=(.25, .25, .25),
                               lw=.25)

                # self.axsk.legend(sklist,
                #                 stationlist,
                #                 loc=4,
                #                 ncol=2,
                #                 markerscale=1,
                #                 borderaxespad=.01,
                #                 labelspacing=.07,
                #                 handletextpad=.2,
                #                 borderpad=.02)
                if pdict['skew'] != nrows - 1:
                    plt.setp(self.axsk.xaxis.get_ticklabels(), visible=False)
                    self.axsk.set_xlabel(' ')
            # ----set axes properties for pt-----------------------------------
            if self._plot_pt == 'y':
                self.axpt.set_xlim(np.floor(np.log10(self.xlimits[0])) * \
                                   self.ellipse_spacing,
                                   np.ceil(np.log10(self.xlimits[1])) * \
                                   self.ellipse_spacing)

                tklabels = []
                xticks = []
                for tk in self.axpt.get_xticks():
                    try:
                        tklabels.append(mtpl.labeldict[tk / self.ellipse_spacing])
                        xticks.append(tk)
                    except KeyError:
                        pass
                self.axpt.set_xticks(xticks)
                self.axpt.set_xticklabels(tklabels,
                                          fontdict={'size': self.font_size})
                self.axpt.set_xlabel('Period (s)', fontdict=fontdict)
                self.axpt.set_ylim(ymin=-1.5 * self.ellipse_size,
                                   ymax=1.5 * self.ellipse_size * (ii + 1))

                self.axpt.grid(True,
                               alpha=.25,
                               which='major',
                               color=(.25, .25, .25),
                               lw=.25)

                plt.setp(self.axpt.get_yticklabels(), visible=False)
                if pdict['pt'] != nrows - 1:
                    plt.setp(self.axpt.get_xticklabels(), visible=False)

                # add colorbar for PT
                axpos = self.axpt.get_position()
                cb_position = (axpos.bounds[0] - .0575,
                               axpos.bounds[1] + .02,
                               .01,
                               axpos.bounds[3] * .75)
                self.cbax = self.fig.add_axes(cb_position)
                if self.ellipse_cmap == 'mt_seg_bl2wh2rd':
                    # make a color list
                    clist = [(cc, cc, 1)
                             for cc in np.arange(0, 1 + 1. / (nseg), 1. / (nseg))] + \
                            [(1, cc, cc)
                             for cc in np.arange(1, -1. / (nseg), -1. / (nseg))]

                    # make segmented colormap
                    mt_seg_bl2wh2rd = colors.ListedColormap(clist)

                    # make bounds so that the middle is white
                    bounds = np.arange(ckmin - ckstep, ckmax + 2 * ckstep, ckstep)

                    # normalize the colors
                    norms = colors.BoundaryNorm(bounds, mt_seg_bl2wh2rd.N)

                    # make the colorbar
                    self.cbpt = mcb.ColorbarBase(self.cbax,
                                                 cmap=mt_seg_bl2wh2rd,
                                                 norm=norms,
                                                 orientation='vertical',
                                                 ticks=bounds[1:-1])
                else:
                    self.cbpt = mcb.ColorbarBase(self.cbax,
                                                 cmap=mtcl.cmapdict[cmap],
                                                 norm=colors.Normalize(vmin=ckmin,
                                                                       vmax=ckmax),
                                                 orientation='vertical')
                self.cbpt.set_ticks([ckmin, (ckmax - ckmin) / 2, ckmax])
                self.cbpt.set_ticklabels(['{0:.0f}'.format(ckmin),
                                          '{0:.0f}'.format((ckmax - ckmin) / 2),
                                          '{0:.0f}'.format(ckmax)])
                self.cbpt.ax.yaxis.set_label_position('left')
                self.cbpt.ax.yaxis.set_label_coords(-1.05, .5)
                self.cbpt.ax.yaxis.tick_right()
                self.cbpt.ax.tick_params(axis='y', direction='in')
                self.cbpt.set_label(mtpl.ckdict[self.ellipse_colorby],
                                    fontdict={'size': self.font_size})

                if pdict['pt'] != nrows - 1:
                    plt.setp(self.axpt.xaxis.get_ticklabels(), visible=False)
                    self.axpt.set_xlabel(' ')
            if show:
                plt.show()

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

        plt.close('all')
        self.plot()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return "Plots resistivity and phase for the different modes of the MT \n" + \
               "response for multiple sites. At the moment it supports the \n" + \
               "input of an .edi file. Other formats that will be supported\n" + \
               "are the impedance tensor and errors with an array of periods\n" + \
               "and .j format.\n"
