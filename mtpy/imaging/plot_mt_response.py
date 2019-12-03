# -*- coding: utf-8 -*-
"""
=================
plot_mt_response
=================

Plots the resistivity and phase for different modes and components

Created 2017

@author: jpeacock
"""
# ==============================================================================
import os

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.colorbar as mcb
import matplotlib.gridspec as gridspec

import mtpy.imaging.mtcolors as mtcl

# ==============================================================================
#  Plot apparent resistivity and phase
# ==============================================================================
from mtpy import MtPyLog
from mtpy.imaging.mtplottools import PlotSettings


class PlotMTResponse(PlotSettings):
    """
    Plots Resistivity and phase for the different modes of the MT response.  At
    the moment it supports the input of an .edi file. Other formats that will
    be supported are the impedance tensor and errors with an array of periods
    and .j format.

    The normal use is to input an .edi file, however it would seem that not
    everyone uses this format, so you can input the data and put it into 
    arrays or objects like class mtpy.core.z.Z.  Or if the data is in 
    resistivity and phase format they can be input as arrays or a class
    mtpy.imaging.mtplot.ResPhase.  Or you can put it into a class
    mtpy.imaging.mtplot.MTplot.

    The plot places the apparent resistivity in log scale in the top panel(s), 
    depending on the plot_num.  The phase is below this, note that 180 degrees
    has been added to the yx phase so the xy and yx phases plot in the same
    quadrant.  Both the resistivity and phase share the same x-axis which is in
    log period, short periods on the left to long periods on the right.  So
    if you zoom in on the plot both plots will zoom in to the same 
    x-coordinates.  If there is tipper information, you can plot the tipper
    as a third panel at the bottom, and also shares the x-axis.  The arrows are
    in the convention of pointing towards a conductor.  The xx and yy 
    components can be plotted as well, this adds two panels on the right.  
    Here the phase is left unwrapped.  Other parameters can be added as 
    subplots such as strike, skew and phase tensor ellipses.

    To manipulate the plot you can change any of the attributes listed below
    and call redraw_plot().  If you know more aout matplotlib and want to 
    change axes parameters, that can be done by changing the parameters in the
    axes attributes and then call update_plot(), note the plot must be open.


    Arguments:
    ----------
        **fn**: string
                       filename containing impedance (.edi) is the only 
                       format supported at the moment

        **z_array**: np.ndarray((nf, 2, 2), dtype='complex')
                impedance tensor with length of nf -> the number of freq
                *default* is None

        **z_err_array**: np.ndarray((nf, 2, 2), dtype='real')
                    impedance tensor error estimates, same shape as z.
                    *default* is None

        **res_array**: np.ndarray((nf, 2, 2))
                        array of resistivity values in linear scale.
                        *default* is None

        **res_err_array**: np.ndarray((nf, 2, 2))
                            array of resistivity error estimates, same shape 
                            as res_array. *default* is None

        **phase_array**: np.ndarray((nf, 2, 2))
                          array of phase values in degrees, same shape as 
                          res_array. *default* is None

        **phase_err_array**: np.ndarray((nf, 2, 2))
                              array of phase error estimates, same shape as 
                              phase_array. *default* is None

        **tipper_array**: np.ndarray((nf, 1, 2), dtype='complex')
                           array of tipper values for tx, ty. *default* is None

        **tipper_err_array**: np.ndarray((nf, 1, 2))
                               array of tipper error estimates, same shape as
                               tipper_array. *default* is None

        **z_object**: class mtpy.core.z.Z
                      object of mtpy.core.z.  If this is input be sure the
                      attribute z.freq is filled.  *default* is None

        **tipper_object**: class mtpy.core.z.Tipper
                            object of mtpy.core.z. If this is input be sure the
                            attribute z.freq is filled.  
                            *default* is None 

        **mt_object** : class mtpy.imaging.mtplottools.MTplot
                        object of mtpy.imaging.mtplottools.MTplot
                        *default* is None

    Optional Key Words:
    --------------------

        *fig_num*: int
                     figure number
                     *default* is 1

        *fig_size*: [width, height] in inches of actual figure size

        *ffactor*: float
                      scaling factor for computing resistivity from 
                      impedances.
                      *Default* is 1

        *rotation_angle*: float
                   rotation angle of impedance tensor (deg or radians), 
                   *Note* : rotaion is clockwise positive
                   *default* is 0

        *plot_num*: [ 1 | 2 | 3 ]
                        * 1 for just Ex/By and Ey/Bx *default*
                        * 2 for all 4 components
                        * 3 for off diagonal plus the determinant

        *plot_title*: string
                         plot_title of plot
                         *default* is station name

        *plot_tipper*: [ 'yri' | 'yr' | 'yi' | 'n' ]
                          Plots the tipper in a bottom pannel
                          * 'yri'  --> plots the real and imaginar parts
                          * 'yr'   --> plots just the real part
                          * 'yi'   --> plots just the imaginary part

                          *Note:* the convention is to point towards a 
                          conductor.  Can change this by setting the
                          parameter arrow_direction = 1.

        *plot_strike*: [ 'y' | 1 | 2 | 3 | 'n' ]
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

        *plot_skew*: [ 'y' | 'n' ]
                        plots the skew angle calculated from the phase tensor
                            * 'y' --> plots skew angle
                            * 'n' --> does not plot skew angle *default*

        *plot_pt*: [ 'y' | 'n' ]
                    plots the phase tensor ellipses which have the properties
                    of ellipse_
                        * 'y' --> plots phase tensor ellipses
                        * 'n' --> does not plot ellipses *default*

        *fig_dpi*: int
                 dots-per-inch resolution, *default* is 300


        :Example: ::

            >>> import mtpy.core.mt.MT
            >>> import mtpy.imaging.plot_mt_response.PlotMTResponse
            >>> edifile = r"/home/MT01/MT01.edi"
            >>> mt_obj = MT(edifile)
            >>> rp1 = PlotMTResponse(mt_obj.Z, plot_num=2)
            >>> # plots all 4 components
            >>> rp1 = PlotMTResponse(mt_obj.Z,  t_object=mt_obj.Tipper,  plot_tipper='yr')
            >>> # plots the real part of the tipper

    Attributes:
    -----------
        -fn              filename to be plotted (only supports .edi so far) 
        -fig_num         figure number for plotting
        -fig_size        size of figure in inches [width, height]
        -plot_num        plot type, see arguments for details 
        -plot_title      title of the plot, *default* is station name
        -fig_dpi         Dots-per-inch resolution of plot, *default* is 300
        -rotation_angle           Rotate impedance tensor by this angle (deg) assuming
                         that North is 0 and angle is positive clockwise

        -plot_tipper    string to tell the program to plot tipper arrows or 
                        not, see accepted values above in arguments

        -plot_strike    string or integer telling the program to plot the 
                        strike angle, see values above in arguments  (YG: not implemented)

        -plot_skew      string to tell the program to plot skew angle.
                        The skew is plotted in the same subplot as the strike
                        angle at the moment  (YG: not implemented)


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
        -marker_lw       line width of marker, *default* is 100./fig_dpi
        -lw              line width of line and errorbar lines
        ..

         *For more on line and marker styles see matplotlib.lines.Line2D*

        -arrow_lw          line width of the arrow, *default* is 0.75
        -arrow_head_width  head width of the arrow, *default* is 0 for no arrow
                           head.  Haven't found a good way to scale the arrow
                           heads in a log scale.

        -arrow_head_length  head width of the arrow, *default* is 0 for no arrow
                            head.  Haven't found a good way to scale the arrow
                            heads in a log scale.

        -arrow_color_real  color of the real arrows, *default* is black
        -arrow_color_imag  color of the imaginary arrows, *default* is blue

        -arrow_direction   0 for pointing towards a conductor and -1 for 
                           pointing away from a conductor.


        -x_limits        limits on the x-limits (period), *default* is None
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

        -phase_quadrant [ 1 | 3 ] 
                        * 1 for both phases to be in 0 to 90, 
                        * 3 for xy to be in 0-90 and yx to be in -180 to 270  

        -tipper_limits  limits of the y-axis, *default* is (-1,1)

        -skew_limits    limits for skew angle, *default* is (-9,9)

        -strike_limits  limits for strike angle, *default* is (-90,90) 


    Methods:
    --------
        * *plot*: plots the pseudosection according to keywords
        * *redraw_plot*: redraws the plot, use if you change some of the 
                         attributes.
        * *update_plot*: updates the plot, use if you change some of the 
                         axes attributes, figure needs to be open to update.
        * *save_plot*: saves the plot to given filepath.


    """

    def __init__(self, z_object=None, t_object=None, pt_obj=None,
                 station='MT Response', **kwargs):
        super(PlotMTResponse, self).__init__()
        self._logger = MtPyLog.get_mtpy_logger(self.__class__.__module__ + "." + self.__class__.__name__)
        self.Z = z_object
        self.Tipper = t_object
        self.pt = pt_obj
        self.station = station

        self.phase_quadrant = 1

        self.plot_num = kwargs.pop('plot_num', 1)
        self.rotation_angle = kwargs.pop('rotation_angle', 0)

        if self.Tipper is not None:
            self.plot_tipper = 'yri'
        else:
            self.plot_tipper = 'n'

        if self.pt is not None:
            self.plot_pt = 'y'
        else:
            self.plot_pt = 'n'

        # set arrow properties
        self.arrow_size = 1
        self.arrow_head_length = 0.03
        self.arrow_head_width = 0.03
        self.arrow_lw = .5
        self.arrow_threshold = 2
        self.arrow_color_imag = 'b'
        self.arrow_color_real = 'k'
        self.arrow_direction = 0

        # ellipse_properties
        self.ellipse_size = 0.25
        self.ellipse_range = (0, 90, 10)
        self.ellipse_colorby = 'phimin'
        self.ellipse_cmap = 'mt_bl2gr2rd'
        self.ellipse_spacing = kwargs.pop('ellipse_spacing', 1)
        if self.ellipse_size == 2 and self.ellipse_spacing == 1:
            self.ellipse_size = 0.25

            # figure properties:
        self.fig_num = 1
        self.fig_dpi = 150
        self.fig_size = None

        self.font_size = 7
        self.marker_size = 3
        self.marker_lw = .5
        self.lw = .5
        self.plot_title = None

        # line styles:
        self.xy_ls = ':'
        self.yx_ls = ':'
        self.det_ls = ':'

        # marker styles:
        self.xy_marker = 's'
        self.yx_marker = 'o'
        self.det_marker = 'v'

        # marker color styles:
        self.xy_color = (0, 0, .75)
        self.yx_color = (.75, 0, 0)
        self.det_color = (0, .75, 0)

        # marker face color styles:
        self.xy_mfc = (0, 0, .75)
        self.yx_mfc = (.75, 0, 0)
        self.det_mfc = (0, .75, 0)

        # plot limits
        self.x_limits = None
        self.res_limits = None
        self.phase_limits = None
        self.tipper_limits = None
        self.pt_limits = None

        # layout params
        self.show_resphase_xticklabels = False

        self.plot_yn = 'y'

        for key in list(kwargs.keys()):
            if hasattr(self, key):
                setattr(self, key, kwargs[key])
            else:
                self._logger.warn("Argument {}={} is not supported thus not been set.".format(key, kwargs[key]))

        # plot on initializing
        if self.plot_yn == 'y':
            self.plot()

    @property
    def period(self):
        """
        plot period
        """
        if self.Z is not None:
            return 1. / self.Z.freq
        elif self.Tipper is not None:
            return 1. / self.Tipper.freq
        else:
            return None

    def plot(self, show=True, overlay_mt_obj=None):
        """
        plotResPhase(filename,fig_num) will plot the apparent resistivity and 
        phase for a single station. 

        """
        if overlay_mt_obj is None: # original case, no overlay edis
            Z2 = None
        else:
            Z2 = overlay_mt_obj.Z

        label_dict = dict([(ii, '$10^{' + str(ii) + '}$') for ii in range(-20, 21)])
        ckdict = {'phiminang': r'$\Phi_{min}$ (deg)',
                  'phimin': r'$\Phi_{min}$ (deg)',
                  'phimaxang': r'$\Phi_{max}$ (deg)',
                  'phimax': r'$\Phi_{max}$ (deg)',
                  'phidet': r'Det{$\Phi$} (deg)',
                  'skew': r'Skew (deg)',
                  'normalized_skew': r'Normalized Skew (deg)',
                  'ellipticity': r'Ellipticity',
                  'skew_seg': r'Skew (deg)',
                  'normalized_skew_seg': r'Normalized Skew (deg)',
                  'geometric_mean': r'$\sqrt{\Phi_{min} \cdot \Phi_{max}}$'}

        if self.plot_tipper.find('y') == 0:
            if self.Tipper is None or np.all(self.Tipper.tipper == 0 + 0j):
                print('No Tipper data for station {0}'.format(self.station))
                self.plot_tipper = 'n'

        if self.plot_pt == 'y':
            #if np.all(self.Z.z == 0 + 0j) or self.Z is None:
            if self.pt is None: # no phase tensor object provided
                print('No Tipper data for station {0}'.format(self.station))
                self.plot_pt = 'n'

        # set x-axis limits from short period to long period
        if self.x_limits is None:
            self.x_limits = (10 ** (np.floor(np.log10(self.period.min()))),
                             10 ** (np.ceil(np.log10((self.period.max())))))
        if self.phase_limits is None:
            pass

        if self.res_limits is None:
            self.res_limits = (10 ** (np.floor(
                np.log10(min([np.nanmin(self.Z.res_xy),
                              np.nanmin(self.Z.res_yx)])))),
                               10 ** (np.ceil(
                                   np.log10(max([np.nanmax(self.Z.res_xy),
                                                 np.nanmax(self.Z.res_yx)])))))

        # set some parameters of the figure and subplot spacing
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.bottom'] = .1
        plt.rcParams['figure.subplot.top'] = .93
        plt.rcParams['figure.subplot.left'] = .80
        plt.rcParams['figure.subplot.right'] = .98

        # set the font properties for the axis labels
        fontdict = {'size': self.font_size + 2, 'weight': 'bold'}

        # create a dictionary for the number of subplots needed
        pdict = {'res': 0,
                 'phase': 1}
        # start the index at 2 because resistivity and phase is permanent for
        # now
        index = 2
        if self.plot_tipper.find('y') >= 0:
            pdict['tip'] = index
            index += 1
        if self.plot_pt.find('y') >= 0:
            pdict['pt'] = index
            index += 1

        # get number of rows needed
        nrows = index

        # set height ratios of the subplots
        hr = [2, 1.5] + [1] * (len(list(pdict.keys())) - 2)

        # create a grid to place the figures into, set to have 2 rows and 2
        # columns to put any of the 4 components.  Make the phase plot
        # slightly shorter than the apparent resistivity plot and have the two
        # close to eachother vertically.  If there is tipper add a 3rd row and
        # if there is strike add another row
        gs = gridspec.GridSpec(nrows, 2, height_ratios=hr, hspace=.05)

        # make figure instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)

        # --> make figure for xy,yx components
        if self.plot_num == 1 or self.plot_num == 3:
            # set label coordinates
            labelcoords = (-0.075, 0.5)

            # space out the subplots
            gs.update(hspace=.05, wspace=.15, left=.1)

            # --> create the axes instances
            # apparent resistivity axis
            self.axr = self.fig.add_subplot(gs[0, :])

            # phase axis that shares period axis with resistivity
            self.axp = self.fig.add_subplot(gs[1, :], sharex=self.axr)



        # --> make figure for all 4 components
        elif self.plot_num == 2:
            # set label coordinates
            labelcoords = (-0.095, 0.5)

            # space out the subplots
            gs.update(hspace=.05, wspace=.15, left=.07)

            # --> create the axes instances
            # apparent resistivity axis
            self.axr = self.fig.add_subplot(gs[0, 0])

            # phase axis that shares period axis with resistivity
            self.axp = self.fig.add_subplot(gs[1, 0], sharex=self.axr)

        # place y coordinate labels in the same location
        self.axr.yaxis.set_label_coords(labelcoords[0], labelcoords[1])
        self.axp.yaxis.set_label_coords(labelcoords[0], labelcoords[1])

        # --> plot tipper
        try:
            self.axt = self.fig.add_subplot(gs[pdict['tip'], :], )
            self.axt.yaxis.set_label_coords(labelcoords[0], labelcoords[1])
        except KeyError:
            pass

        # --> plot phase tensors
        try:
            # can't share axis because not on the same scale
            self.axpt = self.fig.add_subplot(gs[pdict['pt'], :],
                                             aspect='equal')
            self.axpt.yaxis.set_label_coords(labelcoords[0], labelcoords[1])
        except KeyError:
            pass

        nz_xx = np.nonzero(self.Z.z[:, 0, 0])
        nz_xy = np.nonzero(self.Z.z[:, 0, 1])
        nz_yx = np.nonzero(self.Z.z[:, 1, 0])
        nz_yy = np.nonzero(self.Z.z[:, 1, 1])

        if self.Tipper is not None:  # fix github issue #24.
            # NOTE the following lines seems not have any effect anyway
            nz_tx = np.nonzero(self.Tipper.tipper[:, 0, 0])
            nz_ty = np.nonzero(self.Tipper.tipper[:, 0, 1])

        # ---------plot the apparent resistivity--------------------------------
        # --> plot as error bars and just as points xy, yx
        # res_xy
        self.ebxyr = self.axr.errorbar(self.period[nz_xy],
                                       self.Z.res_xy[nz_xy],
                                       marker=self.xy_marker,
                                       ms=self.marker_size,
                                       mew=self.lw,
                                       mec=self.xy_color,
                                       color=self.xy_color,
                                       ecolor=self.xy_color,
                                       ls=self.xy_ls,
                                       lw=self.lw,
                                       yerr=self.Z.res_err_xy[nz_xy],
                                       capsize=self.marker_size,
                                       capthick=self.lw)
        # FZ: overlay edi logic
        if(Z2 is not None):
            self.ebxyr2=self.axr.errorbar(self.period[nz_xy],
                                       Z2.res_xy[nz_xy],
                                       marker=self.xy_marker,
                                       ms=0.5*self.marker_size,
                                       mew=self.lw,
                                       mec=self.xy_color, # (0.5,0.5,0.9),
                                       color=self.xy_color,
                                       ecolor=self.xy_color,
                                       ls=self.xy_ls,
                                       lw=0.5*self.lw,
                                       yerr=Z2.res_err_xy[nz_xy],
                                       capsize=self.marker_size,
                                       capthick=self.lw)

        # res_yx
        self.ebyxr = self.axr.errorbar(self.period[nz_yx],
                                       self.Z.res_yx[nz_yx],
                                       marker=self.yx_marker,
                                       ms=self.marker_size,
                                       mew=self.lw,
                                       mec=self.yx_color,
                                       color=self.yx_color,
                                       ecolor=self.yx_color,
                                       ls=self.yx_ls,
                                       lw=self.lw,
                                       yerr=self.Z.res_err_yx[nz_yx],
                                       capsize=self.marker_size,
                                       capthick=self.lw)

        if (Z2 is not None):
            self.ebyxr2 = self.axr.errorbar(self.period[nz_yx],
                                           Z2.res_yx[nz_yx],
                                           marker=self.yx_marker,
                                           ms=0.5*self.marker_size,
                                           mew=self.lw,
                                           mec=self.yx_color,
                                           color=self.yx_color,
                                           ecolor=self.yx_color,
                                           ls=self.yx_ls,
                                           lw=self.lw,
                                           yerr=Z2.res_err_yx[nz_yx],
                                           capsize=0.5*self.marker_size,
                                           capthick=self.lw)

        # --> set axes properties
        plt.setp(self.axr.get_xticklabels(), visible=False)
        self.axr.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                            fontdict=fontdict)
        self.axr.set_yscale('log', nonposy='clip')
        self.axr.set_xscale('log', nonposx='clip')
        self.axr.set_xlim(self.x_limits)
        self.axr.set_ylim(self.res_limits)
        self.axr.grid(True, alpha=.25, which='both', color=(.25, .25, .25),
                      lw=.25)

        self.axr.legend((self.ebxyr[0], self.ebyxr[0]),
                        ('$Z_{xy}$', '$Z_{yx}$'),
                        loc=3,
                        markerscale=1,
                        borderaxespad=.01,
                        labelspacing=.07,
                        handletextpad=.2,
                        borderpad=.02)

        if Z2 is not None:
            self.axr.legend((self.ebxyr[0], self.ebyxr[0], self.ebxyr2, self.ebyxr2),
                            ('$Z_{xy}$', '$Z_{yx}$', '$Z2_{xy}$','$Z2_{yx}$'),
                            loc=3,
                            markerscale=1,
                            borderaxespad=.01,
                            labelspacing=.07,
                            handletextpad=.2,
                            borderpad=.02)
        # -----Plot the phase---------------------------------------------------
        # phase_xy
        self.ebxyp = self.axp.errorbar(self.period[nz_xy],
                                       self.Z.phase_xy[nz_xy],
                                       marker=self.xy_marker,
                                       ms=self.marker_size,
                                       mew=self.lw,
                                       mec=self.xy_color,
                                       color=self.xy_color,
                                       ecolor=self.xy_color,
                                       ls=self.xy_ls,
                                       lw=self.lw,
                                       yerr=self.Z.phase_err_xy[nz_xy],
                                       capsize=self.marker_size,
                                       capthick=self.lw)

        # phase_yx: Note add 180 to place it in same quadrant as phase_xy
        self.ebyxp = self.axp.errorbar(self.period[nz_yx],
                                       self.Z.phase_yx[nz_yx] + 180,
                                       marker=self.yx_marker,
                                       ms=self.marker_size,
                                       mew=self.lw,
                                       mec=self.yx_color,
                                       color=self.yx_color,
                                       ecolor=self.yx_color,
                                       ls=self.yx_ls,
                                       lw=self.lw,
                                       yerr=self.Z.phase_err_yx[nz_yx],
                                       capsize=self.marker_size,
                                       capthick=self.lw)

        # check the phase to see if any point are outside of [0:90]
        if self.phase_limits is None:
            if min(self.Z.phase_xy) < 0 or min(self.Z.phase_yx + 180) < 0:
                pymin = min([min(self.Z.phase_xy), min(self.Z.phase_yx)])
                if pymin > 0:
                    pymin = 0
            else:
                pymin = 0

            if max(self.Z.phase_xy) > 90 or max(self.Z.phase_yx + 180) > 90:
                pymax = min([max(self.Z.phase_xy), max(self.Z.phase_yx + 180)])
                if pymax < 91:
                    pymax = 89.9
            else:
                pymax = 89.9

            self.phase_limits = (pymin, pymax)

        # --> set axes properties
        self.axp.set_xlabel('Period (s)', fontdict)
        self.axp.set_ylabel('Phase (deg)', fontdict)
        self.axp.set_xscale('log', nonposx='clip')
        self.axp.set_ylim(self.phase_limits)
        self.axp.yaxis.set_major_locator(MultipleLocator(15))
        self.axp.yaxis.set_minor_locator(MultipleLocator(5))
        self.axp.grid(True, alpha=.25, which='both', color=(.25, .25, .25),
                      lw=.25)
        # set th xaxis tick labels to invisible
        if self.plot_tipper.find('y') >= 0 or self.plot_pt == 'y':
            plt.setp(self.axp.xaxis.get_ticklabels(), visible=False)
            self.axp.set_xlabel('')

        # -----plot tipper----------------------------------------------------
        if self.plot_tipper.find('y') == 0:

            txr = self.Tipper.mag_real * np.sin(self.Tipper.angle_real * np.pi / 180 + \
                                                np.pi * self.arrow_direction)
            tyr = self.Tipper.mag_real * np.cos(self.Tipper.angle_real * np.pi / 180 + \
                                                np.pi * self.arrow_direction)

            txi = self.Tipper.mag_imag * np.sin(self.Tipper.angle_imag * np.pi / 180 + \
                                                np.pi * self.arrow_direction)
            tyi = self.Tipper.mag_imag * np.cos(self.Tipper.angle_imag * np.pi / 180 + \
                                                np.pi * self.arrow_direction)

            nt = len(txr)

            tiplist = []
            tiplabel = []

            for aa in range(nt):
                xlenr = txr[aa] * np.log10(self.period[aa])
                xleni = txi[aa] * np.log10(self.period[aa])

                # --> plot real arrows
                if self.plot_tipper.find('r') > 0:
                    self.axt.arrow(np.log10(self.period[aa]),
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
                        line1 = self.axt.plot(0, 0, self.arrow_color_real)
                        tiplist.append(line1[0])
                        tiplabel.append('real')

                # --> plot imaginary arrows
                if self.plot_tipper.find('i') > 0:
                    self.axt.arrow(np.log10(self.period[aa]),
                                   0,
                                   xleni,
                                   tyi[aa],
                                   lw=self.arrow_lw,
                                   facecolor=self.arrow_color_imag,
                                   edgecolor=self.arrow_color_imag,
                                   head_width=self.arrow_head_width,
                                   head_length=self.arrow_head_length,
                                   length_includes_head=False)
                    if aa == 0:
                        line2 = self.axt.plot(0, 0, self.arrow_color_imag)
                        tiplist.append(line2[0])
                        tiplabel.append('imag')

            # make a line at 0 for reference
            self.axt.plot(np.log10(self.period), [0] * nt, 'k', lw=.5)

            self.axt.legend(tiplist, tiplabel,
                            loc='upper left',
                            markerscale=1,
                            borderaxespad=.01,
                            labelspacing=.07,
                            handletextpad=.2,
                            borderpad=.1,
                            prop={'size': self.font_size})

            # set axis properties

            self.axt.set_xlim(np.log10(self.x_limits[0]),
                              np.log10(self.x_limits[1]))

            tklabels = []
            xticks = []

            for tk in self.axt.get_xticks():
                try:
                    tklabels.append(label_dict[tk])
                    xticks.append(tk)
                except KeyError:
                    pass
            self.axt.set_xticks(xticks)
            self.axt.set_xticklabels(tklabels,
                                     fontdict={'size': self.font_size})
            self.axt.set_xlabel('Period (s)', fontdict=fontdict)
            # need to reset the x_limits caouse they get reset when calling
            # set_ticks for some reason
            self.axt.set_xlim(np.log10(self.x_limits[0]),
                              np.log10(self.x_limits[1]))

            self.axt.yaxis.set_major_locator(MultipleLocator(.2))
            self.axt.yaxis.set_minor_locator(MultipleLocator(.1))
            self.axt.set_xlabel('Period (s)', fontdict=fontdict)
            self.axt.set_ylabel('Tipper', fontdict=fontdict)

            # self.axt.set_xscale('log', nonposx='clip')
            if self.tipper_limits is None:
                tmax = max([np.nanmax(tyr), np.nanmax(tyi)])
                if tmax > 1:
                    tmax = .899

                tmin = min([np.nanmin(tyr), np.nanmin(tyi)])
                if tmin < -1:
                    tmin = -.899

                self.tipper_limits = (tmin - .1, tmax + .1)

            self.axt.set_ylim(self.tipper_limits)
            self.axt.grid(True, alpha=.25, which='both', color=(.25, .25, .25),
                          lw=.25)

            # set th xaxis tick labels to invisible
            if self.plot_pt == 'y':
                plt.setp(self.axt.xaxis.get_ticklabels(), visible=False)
                self.axt.set_xlabel('')

        # ----plot phase tensor ellipse---------------------------------------
        if self.plot_pt == 'y':

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
                colorarray = self.pt.phimin

            elif self.ellipse_colorby == 'phimaxang' or \
                    self.ellipse_colorby == 'phimax':
                colorarray = self.pt.phimax


            elif self.ellipse_colorby == 'phidet':
                colorarray = np.sqrt(abs(self.pt.det)) * (180 / np.pi)


            elif self.ellipse_colorby == 'skew' or \
                    self.ellipse_colorby == 'skew_seg':
                colorarray = self.pt.beta

            elif self.ellipse_colorby == 'ellipticity':
                colorarray = self.pt.ellipticity

            else:
                raise NameError(self.ellipse_colorby + ' is not supported')

            # -------------plot ellipses-----------------------------------
            for ii, ff in enumerate(self.period):
                # make sure the ellipses will be visable
                eheight = self.pt.phimin[ii] / self.pt.phimax[ii] * \
                          self.ellipse_size
                ewidth = self.pt.phimax[ii] / self.pt.phimax[ii] * \
                         self.ellipse_size

                # create an ellipse scaled by phimin and phimax and oriented
                # along the azimuth which is calculated as clockwise but needs
                # to be plotted counter-clockwise hence the negative sign.
                ellipd = patches.Ellipse((np.log10(ff) * self.ellipse_spacing,
                                          0),
                                         width=ewidth,
                                         height=eheight,
                                         angle=90 - self.pt.azimuth[ii])

                self.axpt.add_patch(ellipd)

                # get ellipse color
                if cmap.find('seg') > 0:
                    ellipd.set_facecolor(mtcl.get_plot_color(colorarray[ii],
                                                             self.ellipse_colorby,
                                                             cmap,
                                                             ckmin,
                                                             ckmax,
                                                             bounds=bounds))
                else:
                    ellipd.set_facecolor(mtcl.get_plot_color(colorarray[ii],
                                                             self.ellipse_colorby,
                                                             cmap,
                                                             ckmin,
                                                             ckmax))

            # ----set axes properties-----------------------------------------------
            # --> set tick labels and limits
            self.axpt.set_xlim(np.log10(self.x_limits[0]),
                               np.log10(self.x_limits[1]))

            tklabels = []
            xticks = []
            for tk in self.axpt.get_xticks():
                try:
                    tklabels.append(label_dict[tk])
                    xticks.append(tk)
                except KeyError:
                    pass
            self.axpt.set_xticks(xticks)
            self.axpt.set_xticklabels(tklabels,
                                      fontdict={'size': self.font_size})
            self.axpt.set_xlabel('Period (s)', fontdict=fontdict)
            self.axpt.set_ylim(ymin=-1.5 * self.ellipse_size,
                               ymax=1.5 * self.ellipse_size)
            # need to reset the x_limits caouse they get reset when calling
            # set_ticks for some reason
            self.axpt.set_xlim(np.log10(self.x_limits[0]),
                               np.log10(self.x_limits[1]))
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
            if cmap == 'mt_seg_bl2wh2rd':
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
            self.cbpt.set_label(ckdict[self.ellipse_colorby],
                                fontdict={'size': self.font_size})

        # ===Plot the xx, yy components if desired==============================
        if self.plot_num == 2:
            # ---------plot the apparent resistivity----------------------------
            self.axr2 = self.fig.add_subplot(gs[0, 1], sharex=self.axr)
            self.axr2.yaxis.set_label_coords(-.1, 0.5)

            # res_xx
            self.ebxxr = self.axr2.errorbar(self.period[nz_xx],
                                            self.Z.res_xx[nz_xx],
                                            marker=self.xy_marker,
                                            ms=self.marker_size,
                                            mew=self.lw,
                                            mec=self.xy_color,
                                            color=self.xy_color,
                                            ecolor=self.xy_color,
                                            ls=self.xy_ls,
                                            lw=self.lw,
                                            yerr=self.Z.res_err_xx[nz_xx],
                                            capsize=self.marker_size,
                                            capthick=self.lw)

            # res_yy
            self.ebyyr = self.axr2.errorbar(self.period[nz_yy],
                                            self.Z.res_yy[nz_yy],
                                            marker=self.yx_marker,
                                            ms=self.marker_size,
                                            mew=self.lw,
                                            mec=self.yx_color,
                                            color=self.yx_color,
                                            ecolor=self.yx_color,
                                            ls=self.yx_ls,
                                            lw=self.lw,
                                            yerr=self.Z.res_err_yy[nz_yy],
                                            capsize=self.marker_size,
                                            capthick=self.lw)
            if Z2 is not None:
                # res_xx of Z2, with smaller marker size
                self.ebxxr2 = self.axr2.errorbar(self.period[nz_xx],
                                                Z2.res_xx[nz_xx],
                                                marker=self.xy_marker,
                                                ms=0.5*self.marker_size,
                                                mew=self.lw,
                                                mec=self.xy_color,
                                                color=self.xy_color,
                                                ecolor=self.xy_color,
                                                ls=self.xy_ls,
                                                lw=self.lw,
                                                yerr=Z2.res_err_xx[nz_xx],
                                                capsize=0.5*self.marker_size,
                                                capthick=self.lw)

                # res_yy
                self.ebyyr2 = self.axr2.errorbar(self.period[nz_yy],
                                                Z2.res_yy[nz_yy],
                                                marker=self.yx_marker,
                                                ms=0.5*self.marker_size,
                                                mew=self.lw,
                                                mec=self.yx_color,
                                                color=self.yx_color,
                                                ecolor=self.yx_color,
                                                ls=self.yx_ls,
                                                lw=self.lw,
                                                yerr=Z2.res_err_yy[nz_yy],
                                                capsize=0.5*self.marker_size,
                                                capthick=self.lw)

            # --> set axes properties
            plt.setp(self.axr2.get_xticklabels(), visible=False)
            self.axr2.set_yscale('log', nonposy='clip')
            self.axr2.set_xscale('log', nonposx='clip')
            self.axr2.set_xlim(self.x_limits)
            self.axr2.grid(True, alpha=.25,
                           which='both',
                           color=(.25, .25, .25),
                           lw=.25)

            if Z2 is None:
                self.axr2.legend((self.ebxxr[0], self.ebyyr[0]),
                                 ('$Z_{xx}$', '$Z_{yy}$'),
                                 loc=3, markerscale=1,
                                 borderaxespad=.01,
                                 labelspacing=.07,
                                 handletextpad=.2,
                                 borderpad=.02)

            else:
                self.axr2.legend((self.ebxxr[0], self.ebyyr[0], self.ebxxr2[0], self.ebyyr2[0]),
                                 ('$Z_{xx}$', '$Z_{yy}$', '$Z2_{xx}$', '$Z2_{yy}$'),
                                 loc=3, markerscale=1,
                                 borderaxespad=.01,
                                 labelspacing=.07,
                                 handletextpad=.2,
                                 borderpad=.02)

            # -----Plot the phase-----------------------------------------------
            self.axp2 = self.fig.add_subplot(gs[1, 1], sharex=self.axr)

            self.axp2.yaxis.set_label_coords(-.1, 0.5)

            # phase_xx
            self.ebxxp = self.axp2.errorbar(self.period[nz_xx],
                                            self.Z.phase_xx[nz_xx],
                                            marker=self.xy_marker,
                                            ms=self.marker_size,
                                            mew=self.lw,
                                            mec=self.xy_color,
                                            color=self.xy_color,
                                            ecolor=self.xy_color,
                                            ls=self.xy_ls,
                                            lw=self.lw,
                                            yerr=self.Z.phase_err_xx[nz_xx],
                                            capsize=self.marker_size,
                                            capthick=self.lw)

            # phase_yy
            self.ebyyp = self.axp2.errorbar(self.period[nz_yy],
                                            self.Z.phase_yy[nz_yy],
                                            marker=self.yx_marker,
                                            ms=self.marker_size,
                                            mew=self.lw,
                                            mec=self.yx_color,
                                            color=self.yx_color,
                                            ecolor=self.yx_color,
                                            ls=self.yx_ls,
                                            lw=self.lw,
                                            yerr=self.Z.phase_err_yy[nz_yy],
                                            capsize=self.marker_size,
                                            capthick=self.lw)

            # --> set axes properties
            self.axp2.set_xlabel('Period (s)', fontdict)
            self.axp2.set_xscale('log', nonposx='clip')
            self.axp2.set_ylim(ymin=-179.9, ymax=179.9)
            self.axp2.yaxis.set_major_locator(MultipleLocator(30))
            self.axp2.yaxis.set_minor_locator(MultipleLocator(5))
            self.axp2.grid(True, alpha=.25,
                           which='both',
                           color=(.25, .25, .25),
                           lw=.25)

            if len(list(pdict.keys())) > 2:
                plt.setp(self.axp2.xaxis.get_ticklabels(), visible=False)
                plt.setp(self.axp2.xaxis.get_label(), visible=False)

        # ===Plot the Determinant if desired==================================
        if self.plot_num == 3:
            # res_det
            self.ebdetr = self.axr.errorbar(self.period,
                                            self.Z.res_det,
                                            marker=self.det_marker,
                                            ms=self.marker_size,
                                            mew=self.lw,
                                            mec=self.det_color,
                                            color=self.det_color,
                                            ecolor=self.det_color,
                                            ls=self.det_ls,
                                            lw=self.lw,
                                            yerr=self.Z.res_det_err,
                                            capsize=self.marker_size,
                                            capthick=self.lw)
            if Z2 is not None:
                self.ebdetr2 = self.axr.errorbar(self.period,
                                                Z2.res_det,
                                                marker=self.det_marker,
                                                ms=0.5*self.marker_size,
                                                mew=self.lw,
                                                mec=self.det_color,
                                                color=self.det_color,
                                                ecolor=self.det_color,
                                                ls=self.det_ls,
                                                lw=self.lw,
                                                yerr=Z2.res_det_err,
                                                capsize=0.5*self.marker_size,
                                                capthick=self.lw)

            # phase_det
            self.ebdetp = self.axp.errorbar(self.period,
                                            self.Z.phase_det,
                                            marker=self.det_marker,
                                            ms=self.marker_size,
                                            mew=self.lw,
                                            mec=self.det_color,
                                            color=self.det_color,
                                            ecolor=self.det_color,
                                            ls=self.det_ls,
                                            lw=self.lw,
                                            yerr=self.Z.phase_det_err,
                                            capsize=self.marker_size,
                                            capthick=self.lw)

            self.axr.legend((self.ebxyr[0], self.ebyxr[0], self.ebdetr[0]),
                            ('$Z_{xy}$', '$Z_{yx}$', '$\det(\mathbf{\hat{Z}})$'),
                            loc=3,
                            markerscale=1,
                            borderaxespad=.01,
                            labelspacing=.07,
                            handletextpad=.2,
                            borderpad=.02)

            if Z2 is not None:
                self.axr.legend((self.ebxyr[0], self.ebyxr[0], self.ebdetr[0], self.ebxyr2[0], self.ebyxr2[0], self.ebdetr2[0],),
                                ('$Z_{xy}$', '$Z_{yx}$','$\det(\mathbf{\hat{Z}})$', '$Z2_{xy}$','$Z2_{yx}$','$\det(\mathbf{\hat{Z2}})$'),
                                loc=3,
                                markerscale=1,
                                borderaxespad=.01,
                                labelspacing=.07,
                                handletextpad=.2,
                                borderpad=.02)

        if self.show_resphase_xticklabels:
            if self.plot_num in [1, 3]:
                gs.update(hspace=0.2, wspace=.15, left=.1)
            else:
                gs.update(hspace=0.2, wspace=.15, left=.07)
                plt.setp(self.axp2.xaxis.get_ticklabels(), visible=True)
                plt.setp(self.axr2.xaxis.get_ticklabels(), visible=True)
                self.axr2.tick_params(axis='x', pad=2, direction='in', which='both', labelsize=self.font_size - 1)
                self.axp2.tick_params(axis='x', pad=2, direction='in', which='both', labelsize=self.font_size - 1)
                self.axp2.set_xlabel('Period (s)', fontsize=self.font_size - 1, labelpad=0)  #
            plt.setp(self.axr.xaxis.get_ticklabels(), visible=True)
            plt.setp(self.axp.xaxis.get_ticklabels(), visible=True)
            self.axr.tick_params(axis='x', pad=2, direction='in', which='both', labelsize=self.font_size - 1)
            self.axp.tick_params(axis='x', pad=2, direction='in', which='both', labelsize=self.font_size - 1)
        #            self.axp.set_xlabel('Period (s)',fontsize=self.font_size-2,labelpad=0)

        # make plot_title and show
        if self.plot_title is None:
            self.plot_title = self.station

        self.fig.suptitle(self.plot_title, fontdict=fontdict)

        # be sure to show
        if show:
            plt.show()

    def save_plot(self, save_fn, file_format='pdf', orientation='portrait',
                  fig_dpi=None, close_plot='y'):
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
                          saved.  If None then the fig_dpi will be that at 
                          which the figure was made.  I don't think that 
                          it can be larger than fig_dpi of the figure.

            **close_plot** : [ y | n ]
                             * 'y' will close the plot after saving.
                             * 'n' will leave plot open

        :Example: ::

            >>> # to save plot as jpg
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotResPhase(r'/home/MT/mt01.edi')
            >>> p1.save_plot(r'/home/MT/figures', file_format='jpg')

        """

        if fig_dpi is None:
            fig_dpi = self.fig_dpi

        if not os.path.isdir(save_fn):
            file_format = save_fn[-3:]

        else:
            save_fn = os.path.join(save_fn, self.station + '_ResPhase.' +
                                   file_format)

        self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                         orientation=orientation)

        if close_plot == 'y':
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print('Saved figure to: ' + self.fig_fn)

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

        plt.close(self.fig)
        self.plot()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return "Plots Resistivity and phase for the different modes of the" + \
               "MT response."

#########################################################
# Plot the data from one or two EDI files.
# python /c/Githubz/mtpy/mtpy/imaging/plot_mt_response.py /c/Githubz/mtpy/data/edifiles/15125A.edi
# OR
# python /c/Githubz/mtpy/mtpy/imaging/plot_mt_response.py /c/Githubz/mtpy/data/edifiles/15125A.edi /c/Githubz/mtpy/data/edifiles/15129A.edi
####################################################################
if __name__ == "__main__":
    import sys
    from mtpy.core.mt import MT

    if len(sys.argv)<2:
        print("USAGE: python %s edifile1 [edifile2]"%sys.argv[0])
        sys.exit(1)
    elif len(sys.argv) == 2: # one edi file provided
        edifile = sys.argv[1]  # r"C:/Githubz/mtpy/data/edifiles/15125A.edi"
        mt_obj = MT(edifile)

        rp1 = PlotMTResponse(z_object=mt_obj.Z,  # this is mandatory
                             # t_object=mt_obj.Tipper,
                             # pt_obj=mt_obj.pt,
                             station=mt_obj.station,
                             plot_tipper='yr',  # plots the real part of the tipper
                             plot_num=3)  # plot_num =1 xy + yx; 2 all; 3 xy yx det

        # rp1.xy_color = (.5,.5,.9)
        # rp1.xy_marker = '*'
        # rp1.redraw_plot()
    elif(len(sys.argv)==3):    # overlay 2 edi files provided
        edifile = sys.argv[1]  # r"C:/Githubz/mtpy/data/edifiles/15125A.edi"
        mt_obj = MT(edifile)
        edifile2 = sys.argv[2]  # r"C:/Githubz/mtpy/data/edifiles/15126A.edi"
        mt_obj2 = MT(edifile2)

        rp1 = PlotMTResponse(z_object=mt_obj.Z,  # this is mandatory
                             # t_object=mt_obj.Tipper,
                             # pt_obj=mt_obj.pt,
                             station=mt_obj.station ,
                             plot_yn='n',
                             plot_tipper='yr',  # plots the real part of the tipper
                             plot_num=3)  # plot_num =1 xy + yx; 2 all; 3 xy yx det


        rp1.station = rp1.station + " and " + mt_obj2.station
        rp1.plot(overlay_mt_obj=mt_obj2)


