# -*- coding: utf-8 -*-
"""
Created on Thu May 30 17:07:50 2013

@author: jpeacock-pr
"""

#==============================================================================

import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import MultipleLocator
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.colorbar as mcb
import mtpy.imaging.mtcolors as mtcl
import mtpy.imaging.mtplottools as mtpl
#reload(mtpl)
#==============================================================================


class PlotPhaseTensor(mtpl.MTEllipse):
    """
    Will plot phase tensor, strike angle, min and max phase angle,
    azimuth, skew, and ellipticity as subplots on one plot.  It can plot
    the resistivity tensor along side the phase tensor for comparison.

    Arguments:
    ----------

        **fn** : string
               filename containing impedance (.edi) is the only
               format supported at the moment

        **z_object** : class mtpy.core.z.Z
                      object of mtpy.core.z.  If this is input be sure the
                      attribute z.freq is filled.  *default* is None

        **mt_object** : class mtpy.imaging.mtpl.MTplot.mtpl.MTplot
                        object of mtpy.imaging.mtpl.MTplot.mtpl.MTplot
                        *default* is None

        **pt_object** : class mtpy.analysis.pt
                        phase tensor object of mtpy.analysis.pt.  If this is
                        input then the ._mt attribute is set to None cause
                        at the moment cannot tranform the phase tensor to z
                        *default* is None

        **fignum** : int (figure number)

        **rot_z** : float (angle in degrees)
                     rotation angle clockwise positive assuming 0 is North.
                     *Default* is 0

        **plot_yn** : [ 'y' | 'n' ]


        **dpi** : int
                  Dots-per-inch resolution of figure.
                  *Default* is 300

        **ellipse_dict** : dictionary
                          dictionary of parameters for the phase tensor
                          ellipses with keys:
                          * 'size' -> size of ellipse in points
                                     *default* is .25

                          * 'colorby' : [ 'phimin' | 'phimax' | 'beta' |
                                    'skew_seg' | 'phidet' | 'ellipticity' ]

                                    - 'phimin' -> colors by minimum phase
                                    - 'phimax' -> colors by maximum phase
                                    - 'skew' -> colors by skew
                                    - 'skew_seg' -> colors by skew in
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



    :Example: ::

        #To plot just the phase tensor components
        >>> import mtpy.imaging.mtplot as mtplot
        >>> pt1 = mtplot.plot_pt(fn=r"/home/MT/edifiles/MT01.edi")

    Attributes:
    -----------
        -ax1      matplotlib.axes object for the phase tensor ellipses
        -cbpt     matplotlib.colors.ColorBarBase object for coloring ellipses
        -ax2      matplotlib.axes object for the strike angle
        -ax3      matplotlib.axes object for minimum and maximum phase
        -ax4      matplotlib.axes object for skew angle
        -ax5      matplotlib.axes object for ellipticity

        -font_size  size of font for the axes labels, titles will be +2
        -fignum     number of the figure instance
        -fig_size   size of figure in inches
        -plot_yn    boolean to tell the class to plot on instance creation
        -dpi        dots-per-inch resolution of figure

        -strike_inv_marker   marker for strike determined from invariants
                             in ax2
        -strike_inv_color    color for invariant marker in ax2
        -strike_pt_marker    marker for strike determined from pt in ax2
        -strike_pt_color     color for pt strike in ax2
        -strike_tp_marker    marker for strike determined from Tipper in ax2
        -strike_tp_color     color for tipper strike in ax2

        -ptmin_marker  marker for minimum phase in ax3
        -ptmin_color   color for minimum phase markers in ax3
        -ptmax_marker  marker for maximum phase in ax3
        -ptmax_color   color for maximum phase markers in ax3

        -skew_marker   marker for skew angle determined from pt in ax4
        -skew_color    color for skew angle in ax4

        -ellip_marker  marker for ellipticity determined from pt in ax5
        -ellip_color   color for ellipticity in ax5

        -marker_size   size of the marker in all plots
        -marker_lw     width of face lines for markers in all plots

        -pt_limits      limits on the minimu phase and maximum phase (deg)
        -strike_limits  limits on the strike angle in degrees, note the strike
                        is calculated to go from -90 to 90.
        -skew_limits    limits on skew angles (deg)
        -ellip_limits   limits on ellipticity ratio from [0,1]

        -skew_cutoff    plots a line in ax4 at positive and negative of this
                        value to visually recognize 3D effects
        -ellip_cutoff   plots a line in ax5 to represent the cutoff of 2D

        -ellipse_cmap     color map for coloring ellipses of ax1
        -ellipse_colorby  parameter to color the ellipses by
        -ellipse_range    min and max values for coloring ellipses
        -ellipse_size     scaling factor of ellipses
        -ellipse_spacing  spacing between ellipses

        -mt            mtpy.imaging.mtpl.MTplot.mtpl.MTplot object (PlotPhaseTensor._mt)

    """

    def __init__(self, **kwargs):

        fn = kwargs.pop('fn', None)
        z_object = kwargs.pop('z_object', None)
        mt_object = kwargs.pop('mt_object', None)
        pt_object = kwargs.pop('pt_object', None)

        #--> get mt object
        if fn is not None:
            self._mt = mtpl.MTplot(fn=fn)
        elif z_object is not None:
            self._mt = mtpl.MTplot(z_object=z_object)
        elif mt_object is not None:
            self._mt = mt_object
        elif pt_object is not None:
            self.pt = pt_object
            self._mt = mtpl.MTplot()
            self._mt.freq = self.pt.freq

        self.font_size = kwargs.pop('font_size', 7)
        self.fig_dpi = kwargs.pop('dpi', 300)
        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [8, 8])
        self.rot_z = kwargs.pop('rot_z', 0)
        self.plot_yn = kwargs.pop('plot_yn', 'y')

        self.ptmin_marker = kwargs.pop('ptmin_marker', 'o')
        self.ptmax_marker = kwargs.pop('ptmax_marker', 's')
        self.strike_inv_marker = kwargs.pop('strike_inv_marker', 's')
        self.strike_pt_marker = kwargs.pop('strike_pt_marker', 'o')
        self.strike_tp_marker = kwargs.pop('strike_tp_marker', 'v')
        self.skew_marker = kwargs.pop('skew_marker', 's')
        self.ellip_marker = kwargs.pop('ellip_marker', 's')

        self.ptmin_color = kwargs.pop('ptmin_color', 'r')
        self.ptmax_color = kwargs.pop('ptmax_color', 'b')
        self.strike_inv_color = kwargs.pop('strike_inv_color', 'c')
        self.strike_pt_color = kwargs.pop('strike_pt_color', 'purple')
        self.strike_tp_color = kwargs.pop('strike_tp_color', (.5, .5, 0))
        self.skew_color = kwargs.pop('skew_color', 'g')
        self.ellip_color = kwargs.pop('ellip_color', 'orange')

        self.marker_size = kwargs.pop('marker_size', 2)
        self.marker_lw = kwargs.pop('marker_lw', .5)

        self.pt_limits = kwargs.pop('pt_limits', None)
        self.strike_limits = kwargs.pop('strike_limits', None)
        self.ellip_limits = kwargs.pop('ellip_limits', None)
        self.skew_limits = kwargs.pop('skew_limits', None)

        self.skew_cutoff = kwargs.pop('skew_cutoff', 3)
        self.ellip_cutoff = kwargs.pop('ellip_cutoff', 0.2)

        # read ellipse dict
        ellipse_dict = kwargs.pop('ellipse_dict', None)
        if ellipse_dict is None:
            self._ellipse_dict = {'size': .25}
        else:
            self._ellipse_dict = ellipse_dict

        self._read_ellipse_dict()

        self.ellipse_spacing = kwargs.pop('ellipse_spacing', 1)

        self.cb_position = kwargs.pop('cb_position', (.045, .78, .015, .12))

        if self.plot_yn == 'y':
            self.plot()

    def plot(self):
        """
        plots the phase tensor elements
        """

        # Set plot parameters
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = .1
        plt.rcParams['figure.subplot.right'] = .98
        plt.rcParams['figure.subplot.bottom'] = .1
        plt.rcParams['figure.subplot.top'] = .95
        plt.rcParams['figure.subplot.wspace'] = .21
        plt.rcParams['figure.subplot.hspace'] = .5

        font_dict = {'size': self.font_size, 'weight': 'bold'}
        font_dictt = {'size': self.font_size + 2, 'weight': 'bold'}

        #--> create plot instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()

        # get phase tensor instance
        try:
            self.pt
            self.pt.rotate(self.rot_z)
        except AttributeError:
            self.pt = self._mt.pt
            self.pt.rotate(self.rot_z)
            self.zinv = self._mt.Z.invariants
            self.zinv.rotate(self.rot_z)

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
            colorarray = self.pt.phimin[0]

        elif self.ellipse_colorby == 'phidet':
            colorarray = np.sqrt(abs(self.pt.det[0])) * (180 / np.pi)

        elif self.ellipse_colorby == 'skew' or\
                self.ellipse_colorby == 'skew_seg':
            colorarray = self.pt.beta[0]

        elif self.ellipse_colorby == 'ellipticity':
            colorarray = self.pt.ellipticity[0]

        else:
            raise NameError(self.ellipse_colorby + ' is not supported')

        #-------------plotPhaseTensor-----------------------------------
        self.ax1 = self.fig.add_subplot(3, 1, 1, aspect='equal')
        self._mt._period = 1. / self._mt.freq

        for ii, ff in enumerate(self._mt.period):
            # make sure the ellipses will be visable
            if self.pt.phimax[0][ii] != 0:
                eheight = self.pt.phimin[0][ii] / self.pt.phimax[0][ii] *\
                    self.ellipse_size

                ewidth = self.ellipse_size

            else:
                # tiny instead of nothing
                eheight = 0.01 * self.ellipse_size
                ewidth = 0.01 * self.ellipse_size

            # ewidth = self.pt.phimax[0][ii]/self.pt.phimax[0][ii]*\
            #                                                  self.ellipse_size

            # alternative scaling
            # eheight = self.pt.phimin[0][ii]/max(np.abs(self.pt.phimax[0]))*\
            #                                                   self.ellipse_size
            # ewidth = self.pt.phimax[0][ii]/max(np.abs(self.pt.phimax[0]))*\
            #                                                   self.ellipse_size

            # create an ellipse scaled by phimin and phimax and oriented along
            # the azimuth which is calculated as clockwise but needs to
            # be plotted counter-clockwise hence the negative sign.
            ellipd = patches.Ellipse((np.log10(ff) * self.ellipse_spacing, 0),
                                     width=ewidth,
                                     height=eheight,
                                     angle=90 - self.pt.azimuth[0][ii])

            self.ax1.add_patch(ellipd)

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

        #----set axes properties-----------------------------------------------
        #--> set tick labels and limits
        xlimits = (np.floor(np.log10(self._mt.period[0])),
                   np.ceil(np.log10(self._mt.period[-1])))

        self.ax1.set_xlim(xlimits)
        tklabels = []
        xticks = []
        for tk in self.ax1.get_xticks():
            try:
                tklabels.append(mtpl.labeldict[tk])
                xticks.append(tk)
            except KeyError:
                pass
        self.ax1.set_xticks(xticks)
        self.ax1.set_xticklabels(tklabels, fontdict={'size': self.font_size})
        self.ax1.set_xlabel('Period (s)', fontdict=font_dict)
        self.ax1.set_ylim(ymin=-1.5 * self.ellipse_size,
                          ymax=1.5 * self.ellipse_size)

        self.ax1.grid(True,
                      alpha=.25,
                      which='major',
                      color=(.25, .25, .25),
                      lw=.25)

        plt.setp(self.ax1.get_yticklabels(), visible=False)
        # add colorbar for PT
        self.cbax = self.fig.add_axes(self.cb_position)
        if cmap == 'mt_seg_bl2wh2rd':
            # make a color list
            clist = [(cc, cc, 1) for cc in np.arange(0, 1 + 1. / (nseg), 1. / (nseg))] +\
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
        self.cbpt.set_ticks([ckmin, ckmax])
        self.cbpt.set_ticklabels(['{0:.0f}'.format(ckmin),
                                  '{0:.0f}'.format(ckmax)])
        self.cbpt.ax.yaxis.set_label_position('left')
        self.cbpt.ax.yaxis.set_label_coords(-1.05, .5)
        self.cbpt.ax.yaxis.tick_right()
        self.cbpt.ax.tick_params(axis='y', direction='in')
        self.cbpt.set_label(mtpl.ckdict[self.ellipse_colorby],
                            fontdict={'size': self.font_size, 'weight': 'bold'})

        #---------------plotStrikeAngle-----------------------------------
        self.ax2 = self.fig.add_subplot(3, 2, 3)
        az = self.pt.azimuth[0]
        az_err = self.pt.azimuth[1]

        # put the strike into a coordinate system that goes from -90 to 90
        az[np.where(az > 90)] -= 180
        az[np.where(az < -90)] += 180

        stlist = []
        stlabel = []

        # plot phase tensor strike
        ps2 = self.ax2.errorbar(self._mt.period,
                                az,
                                marker=self.strike_pt_marker,
                                ms=self.marker_size,
                                mfc=self.strike_pt_color,
                                mec=self.strike_pt_color,
                                mew=self.marker_lw,
                                ls='none',
                                yerr=az_err,
                                ecolor=self.strike_pt_color,
                                capsize=self.marker_size,
                                elinewidth=self.marker_lw)

        stlist.append(ps2[0])
        stlabel.append('PT')
        try:
            strike = self.zinv.strike
            strikeerr = np.nan_to_num(self.zinv.strike_err)
            # put the strike into a coordinate system that goes from -90 to 90
            strike[np.where(strike > 90)] = strike[np.where(strike > 90)] - 180
            strike[np.where(strike < -90)
                   ] = strike[np.where(strike < -90)] + 180

            # plot invariant strike
            erxy = self.ax2.errorbar(self._mt.period,
                                     strike,
                                     marker=self.strike_inv_marker,
                                     ms=self.marker_size,
                                     mfc=self.strike_inv_color,
                                     mec=self.strike_inv_color,
                                     mew=self.marker_lw,
                                     ls='none',
                                     yerr=strikeerr,
                                     ecolor=self.strike_inv_color,
                                     capsize=self.marker_size,
                                     elinewidth=self.marker_lw)

            stlist.append(erxy[0])
            stlabel.append('Z_inv')
        except AttributeError:
            print('Could not get z_invariants from pt, input z if desired.')

        if self._mt.tipper is not None:
            # strike from tipper
            tp = self._mt.Tipper
            s3 = tp.angle_real + 90

            # fold to go from -90 to 90
            s3[np.where(s3 > 90)] = s3[np.where(s3 > 90)] - 180
            s3[np.where(s3 < -90)] = s3[np.where(s3 < -90)] + 180

            # plot strike with error bars
            ps3 = self.ax2.errorbar(self._mt.period,
                                    s3,
                                    marker=self.strike_tp_marker,
                                    ms=self.marker_size,
                                    mfc=self.strike_tp_color,
                                    mec=self.strike_tp_color,
                                    mew=self.marker_lw,
                                    ls='none',
                                    yerr=np.zeros_like(s3),
                                    ecolor=self.strike_tp_color,
                                    capsize=self.marker_size,
                                    elinewidth=self.marker_lw)

            stlist.append(ps3[0])
            stlabel.append('Tipper')

        self.ax2.legend(stlist,
                        stlabel,
                        loc='lower left',
                        markerscale=.5 * self.marker_size,
                        borderaxespad=.01,
                        labelspacing=.1,
                        handletextpad=.2,
                        ncol=len(stlist),
                        borderpad=.1,
                        columnspacing=.1)

        leg = plt.gca().get_legend()
        ltext = leg.get_texts()  # all the text.Text instance in the legend
        plt.setp(ltext, fontsize=6)    # the legend text fontsize

        if self.strike_limits is None:
            self.strike_limits = (-89.99, 89.99)

        self.ax2.set_yscale('linear')
        self.ax2.set_xscale('log', nonposx='clip')
        self.ax2.set_xlim(xmax=10**xlimits[-1], xmin=10**xlimits[0])
        self.ax2.set_ylim(self.strike_limits)
        self.ax2.yaxis.set_major_locator(MultipleLocator(20))
        self.ax2.yaxis.set_minor_locator(MultipleLocator(5))
        self.ax2.grid(True, alpha=.25, which='both', color=(.25, .25, .25),
                      lw=.25)
        self.ax2.set_ylabel('Angle (deg)', fontdict=font_dict)
        self.ax2.set_title('Strike', fontdict=font_dictt)

        #---------plot Min & Max Phase-----------------------------------------
        minphi = self.pt.phimin[0]
        minphierr = self.pt.phimin[1]
        maxphi = self.pt.phimax[0]
        maxphierr = self.pt.phimax[1]

        self.ax3 = self.fig.add_subplot(3, 2, 4, sharex=self.ax2)

        ermin = self.ax3.errorbar(self._mt.period,
                                  minphi,
                                  marker=self.ptmin_marker,
                                  ms=self.marker_size,
                                  mfc='None',
                                  mec=self.ptmin_color,
                                  mew=self.marker_lw,
                                  ls='None',
                                  yerr=minphierr,
                                  ecolor=self.ptmin_color,
                                  capsize=self.marker_size,
                                  elinewidth=self.marker_lw)

        ermax = self.ax3.errorbar(self._mt.period,
                                  maxphi,
                                  marker=self.ptmax_marker,
                                  ms=self.marker_size,
                                  mfc='None',
                                  mec=self.ptmax_color,
                                  mew=self.marker_lw,
                                  ls='None',
                                  yerr=maxphierr,
                                  ecolor=self.ptmax_color,
                                  capsize=self.marker_size,
                                  elinewidth=self.marker_lw)

        if self.pt_limits is None:
            self.pt_limits = [min([self.pt.phimax[0].min(),
                                   self.pt.phimin[0].min()]) - 3,
                              max([self.pt.phimax[0].max(),
                                   self.pt.phimin[0].max()]) + 3]
            if self.pt_limits[0] < -10:
                self.pt_limits[0] = -9.9
            if self.pt_limits[1] > 100:
                self.pt_limits[1] = 99.99

        self.ax3.set_xscale('log', nonposx='clip')
        self.ax3.set_yscale('linear')

        self.ax3.legend((ermin[0], ermax[0]),
                        ('$\phi_{min}$', '$\phi_{max}$'),
                        loc='lower left',
                        markerscale=.5 * self.marker_size,
                        borderaxespad=.01,
                        labelspacing=.1,
                        handletextpad=.2,
                        ncol=2,
                        borderpad=.01,
                        columnspacing=.01)

        leg = plt.gca().get_legend()
        ltext = leg.get_texts()  # all the text.Text instance in the legend
        plt.setp(ltext, fontsize=6.5)    # the legend text fontsize

        self.ax3.set_ylim(self.pt_limits)
        self.ax3.grid(True, alpha=.25, which='both', color=(.25, .25, .25),
                      lw=.25)

        self.ax3.set_ylabel('Phase (deg)', fontdict=font_dict)
        self.ax3.set_title('$\mathbf{\phi_{min}}$ and $\mathbf{\phi_{max}}$',
                           fontdict=font_dictt)

        #-----------------------plotSkew---------------------------------------

        skew = self.pt.beta[0]
        skewerr = self.pt.beta[1]

        self.ax4 = self.fig.add_subplot(3, 2, 5, sharex=self.ax2)
        erskew = self.ax4.errorbar(self._mt.period,
                                   skew,
                                   marker=self.skew_marker,
                                   ms=self.marker_size,
                                   mfc='None',
                                   mec=self.skew_color,
                                   mew=self.marker_lw,
                                   ls='None',
                                   yerr=skewerr,
                                   ecolor=self.skew_color,
                                   capsize=self.marker_size,
                                   elinewidth=self.marker_lw)

        # plot lines indicating not 3d
        self.ax4.plot([10**xlimits[0], 10**xlimits[-1]],
                      [self.skew_cutoff, self.skew_cutoff],
                      ls='--',
                      color=self.skew_color,
                      lw=1)

        self.ax4.plot([10**xlimits[0], 10**xlimits[-1]],
                      [-self.skew_cutoff, -self.skew_cutoff],
                      ls='--',
                      color=self.skew_color,
                      lw=1)

        self.ax4.set_xscale('log', nonposx='clip')
        self.ax4.set_yscale('linear')
        self.ax4.yaxis.set_major_locator(MultipleLocator(ckstep))

        if self.skew_limits is None:
            self.skew_limits = (-10, 10)
        self.ax4.set_ylim(self.skew_limits)
        self.ax4.grid(True, alpha=.25, which='both', color=(.25, .25, .25),
                      lw=.25)
        self.ax4.set_xlabel('Period (s)', fontdict=font_dict)
        self.ax4.set_ylabel('Skew Angle (deg)', fontdict=font_dict)
        self.ax4.set_title('Skew Angle', fontdict=font_dictt)

        #----------------------plotEllipticity--------------------------------
        ellipticity = self.pt.ellipticity[0]
        ellipticityerr = self.pt.ellipticity[1]

        self.ax5 = self.fig.add_subplot(3, 2, 6, sharex=self.ax2)
        erskew = self.ax5.errorbar(self._mt.period,
                                   ellipticity,
                                   marker=self.ellip_marker,
                                   ms=self.marker_size,
                                   mfc='None',
                                   mec=self.ellip_color,
                                   mew=self.marker_lw,
                                   ls='None',
                                   yerr=ellipticityerr,
                                   ecolor=self.ellip_color,
                                   capsize=self.marker_size,
                                   elinewidth=self.marker_lw)

        # draw a line where the ellipticity is not 2d
        self.ax5.plot([10**xlimits[0], 10**xlimits[-1]],
                      [self.ellip_cutoff, self.ellip_cutoff],
                      ls='--',
                      color=self.ellip_color,
                      lw=1)

        self.ax5.set_xscale('log', nonposx='clip')
        self.ax5.set_yscale('linear')

        self.ax5.yaxis.set_major_locator(MultipleLocator(.1))

        self.ax5.set_ylim(ymin=0, ymax=1)
        self.ax5.grid(True, alpha=.25, which='both', color=(.25, .25, .25),
                      lw=.25)
        self.ax5.set_xlabel('Period (s)', fontdict=font_dict)
        self.ax5.set_ylabel('$\mathbf{\phi_{max}-\phi_{min}/\phi_{max}+\phi_{min}}$',
                            fontdict=font_dict)
        self.ax5.set_title('Ellipticity', fontdict=font_dictt)

        try:
            self.fig.suptitle('Phase Tensor Elements for: ' + self._mt.station,
                              fontdict={'size': self.font_size + 3, 'weight': 'bold'})
        except:
            self.fig.suptitle('Phase Tensor Elements for Station "unknown"',
                              fontdict={'size': self.font_size + 3, 'weight': 'bold'})

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
                            save_fn/station_name_PhaseTensor.file_format

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

            >>> # to save plot as jpg
            >>> import mtpy.imaging.mtpl.MTplottools as mtpl.MTplot
            >>> p1 = mtpl.MTplot.PlotPhaseTensor(r'/home/MT/mt01.edi')
            >>> p1.save_plot(r'/home/MT/figures', file_format='jpg')

        """

        if fig_dpi is None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation)
            plt.clf()
            plt.close(self.fig)

        else:
            station = self._mt.station
            if station is None:
                station = 'MT01'
            save_fn = os.path.join(save_fn, station + '_PhaseTensor.' +
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
            >>> import mtpy.imaging.mtpl.MTplottools as mtpl.MTplot
            >>> p1 = mtpl.MTplot.PlotResPhase(r'/home/MT/mt01.edi')
            >>> [ax.grid(True, which='major') for ax in [p1.axr,p1.axp]]
            >>> p1.update_plot()

        """

        self.fig.canvas.draw()

    def redraw_plot(self):
        """
        use this function if you updated some attributes and want to re-plot.

        :Example: ::

            >>> # change the color and marker of the xy components
            >>> import mtpy.imaging.mtpl.MTplottools as mtpl.MTplot
            >>> p1 = mtpl.MTplot.PlotResPhase(r'/home/MT/mt01.edi')
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

        return "Plots the phase tensor ellipses and other properties such\n" +\
               "strike angle, minimum and maximum phase, skew and ellipticity"
