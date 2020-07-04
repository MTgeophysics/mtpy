# -*- coding: utf-8 -*-
"""
Created on Thu May 30 18:28:24 2013

@author: jpeacock-pr
"""

#==============================================================================

import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
from mtpy.analysis.zinvariants import Zinvariants
import mtpy.imaging.mtplottools as mtpl

#==============================================================================


class PlotStrike(object):
    """
    PlotStrike will plot the strike estimated from the invariants, phase tensor
    and the tipper in either a rose diagram of xy plot


    plots the strike angle as determined by phase tensor azimuth (Caldwell et
    al. [2004]) and invariants of the impedance tensor (Weaver et al. [2003]).

    The data is split into decades where the histogram for each is plotted in
    the form of a rose diagram with a range of 0 to 180 degrees.
    Where 0 is North and 90 is East.   The median angle of the period band is
    set in polar diagram.  The top row is the strike estimated from
    the invariants of the impedance tensor.  The bottom row is the azimuth
    estimated from the phase tensor.  If tipper is 'y' then the 3rd row is the
    strike determined from the tipper, which is orthogonal to the induction
    arrow direction.

    Arguments
    ----------
        :param fn_list: full paths to .edi files to plot

        :param z_object: object of mtpy.core.z.  If this is input be sure the
                          attribute z.freq is filled.  *default* is None

        :param mt_object: object of mtpy.imaging.mtplot.MTplot
                         *default* is None

        :param fig_num: figure number to be plotted. *Default* is 1

        :param fs: font size for labels of plotting. *Default* is 10

        :param rot_z: angle of rotation clockwise positive. *Default* is 0

        :param period_tolerance: float
                   Tolerance level to match periods from different edi files.
                   *Default* is 0.05

        :param text_pad: padding of the angle label at the bottom of each
                         polar diagram.  *Default* is 1.65

        :param text_size: font size


        :param plot_range: [ 'data' | (period_min,period_max) ]
                    period range to estimate the strike angle. Options are:
                        * *'data'* for estimating the strike for all periods
                            in the data.
                        * (pmin,pmax) for period min and period max, input as
                          (log10(pmin),log10(pmax))

        :param plot_type: [ 1 | 2 ]
                        - *1* to plot individual decades in one plot
                        - *2* to plot all period ranges into one polar diagram
                              for each strike angle estimation

        :param plot_tipper: [ True | False ]
                          - True to plot the tipper strike
                          - False to not plot tipper strike

        :param pt_error_floor: Maximum error in degrees that is allowed to 
                               estimate strike. *Default* is None allowing all 
                               estimates to be used.

        :param fold: [ True | False ]
                    * True to plot only from 0 to 180
                    * False to plot from 0 to 360
                    
        :param plot_orthogonal: [ True | False]
                                * True to plot the orthogonal strike directions
                                * False to not
        
        :param color: [ True | False ]
                      * True to plot shade colors
                      * False to plot all in one color
                      
        :param color_inv: color of invariants plots
        
        :param color_pt: color of phase tensor plots
        
        :param color_tip: color of tipper plots
        
        :param ring_spacing: spacing of rings in polar plots
        
        :param ring_limits: (min count, max count) set each plot have these 
                            limits 
                            
        :param plot_orientation: [ 'h' | 'v' ] horizontal or vertical plots

    :Example: ::

        >>> import glob
        >>> import mtpy.imaging.mtplot as mtplot
        >>> edi_dir = r"/home/EDIFiles"
        >>> edi_list = glob.glob("{0}\*.edi".format(edi_dir)
        >>> #---plot rose plots in decades 
        >>> strike = mtplot.plot_strike(fn_list=edilist, plot_type=1)
        >>> #---Turn on Tipper
        >>> strike.plot_tipper = True
        >>> #---Plot only main directions
        >>> strike.plot_orthogonal = False
        >>> # Redraw plot
        >>> strike.redraw_plot()
        >>> # plot only from 0-180
        >>> strike.fold = True
        >>> strike.redraw_plot()
        >>> #---save the plot---
        >>> strike.save_plot(r"/home/Figures")

    """

    def __init__(self, **kwargs):

        fn_list = kwargs.pop('fn_list', None)
        z_object_list = kwargs.pop('z_object_list', None)
        tipper_object_list = kwargs.pop('tipper_object_list', None)
        mt_object_list = kwargs.pop('mt_object_list', None)
        self._rotation_angle = 0

        #------Set attributes of the class-----------------

        #--> get the inputs into a list of mt objects
        self.mt_list = mtpl.get_mtlist(fn_list=fn_list,
                                       z_object_list=z_object_list,
                                       tipper_object_list=tipper_object_list,
                                       mt_object_list=mt_object_list)

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

        #--> set plot properties
        self.fig_num = 1
        self.fig_dpi = 300
        self.fig_size = [7, 5]

        self.plot_type =  2
        self.plot_title = None
        self.plot_range = 'data'
        self.plot_tipper = False
        self.plot_orientation = 'h'

        self.period_tolerance = .05
        self.pt_error_floor = None
        self.fold = True
        self.bin_width = 5
        self.color = True
        self.color_inv = (.7, 0, .2)
        self.color_pt = (.2, 0, .7)
        self.color_tip = (.2, .65, .2)
        self.ring_spacing = 10
        self.ring_limits = None
        self.plot_orthogonal = True
        
        self.font_size = 7
        self.text_pad = 0.6
        self.text_size = self.font_size
            
        for key, value in kwargs.items():
            setattr(self, key, value)

        # make a dictionary for plotting titles
        self.title_dict = {}
        self.title_dict[-5] = '10$^{-5}$ - 10$^{-4}$ s'
        self.title_dict[-4] = '10$^{-4}$ - 10$^{-3}$ s'
        self.title_dict[-3] = '10$^{-3}$ - 10$^{-2}$ s'
        self.title_dict[-2] = '10$^{-2}$ - 10$^{-1}$ s'
        self.title_dict[-1] = '10$^{-1}$ - 10$^{0}$ s'
        self.title_dict[0] = '10$^{0}$ - 10$^{1}$ s'
        self.title_dict[1] = '10$^{1}$ - 10$^{2}$ s'
        self.title_dict[2] = '10$^{2}$ - 10$^{3}$ s'
        self.title_dict[3] = '10$^{3}$ - 10$^{4}$ s'
        self.title_dict[4] = '10$^{4}$ - 10$^{5}$ s'
        self.title_dict[5] = '10$^{5}$ - 10$^{6}$ s'
        self.title_dict[6] = '10$^{6}$ - 10$^{7}$ s'

        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()

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
            mt.rotation_angle = value
            
        self._rotation_angle = value
            
        self.make_strike_array()
    
    def make_strike_array(self):
        """
        make strike array
        """
        inv_list = []
        pt_list = []
        tip_list = []
        # initialize some parameters
        nc = len(self.mt_list)
        nt = 0
        kk = 0
        
        for dd, mt in enumerate(self.mt_list):
            
            if mt.period.size > nt:
                nt = mt.period.size
            #-----------get strike angle from invariants-----------------------
            zinv = Zinvariants(mt.Z)

            # add 90 degrees because invariants assume 0 is north, but plotting
            # assumes that 90 is north and measures clockwise, thus the negative
            # because the strike angle from invariants is measured
            # counter-clockwise

            zs = 90 - zinv.strike

            # fold so the angle goes from 0 to 180
            if self.fold == True:
                # for plotting put the NW angles into the SE quadrant
                zs[np.where(zs > 90)] -= 180
                zs[np.where(zs < -90)] += 180

            # leave as the total unit circle 0 to 360
            elif self.fold == False:
                zs %= 360
                
            # make a dictionary of strikes with keys as period
            mdictinv = dict([(ff, jj) for ff, jj in zip(mt.period, zs)])
            inv_list.append(mdictinv)

            #------------get strike from phase tensor strike angle-------------
            pt = mt.pt
            az = 90 - pt.azimuth
            az_err = pt.azimuth_err
            az[pt.phimax == 0] = np.nan

            # need to add 90 because pt assumes 0 is north and
            # negative because measures clockwise.
            # put an error max on the estimation of strike angle
            if self.pt_error_floor:
                az[np.where(az_err > self.pt_error_floor)] = 0.0

            # fold so the angle goes from 0 to 180
            if self.fold:
                az[np.where(az > 90)] -= 180
                az[np.where(az < -90)] += 180

            # leave as the total unit circle 0 to 360
            elif self.fold == False:
                az %= 360
                
            # make a dictionary of strikes with keys as period
            mdictpt = dict([(ff, jj) for ff, jj in zip(mt.period, az)])
            pt_list.append(mdictpt)

            #-----------get tipper strike------------------------------------
            tip = mt.Tipper
            if tip.tipper is None:
                tip.tipper = np.zeros((len(mt.period), 1, 2),
                                              dtype='complex')
                tip.compute_components()

            # needs to be negative because measures clockwise
            tipr = -tip.angle_real

            tipr[np.where(tipr == 180.)] = 0.0

            # fold so the angle goes from 0 to 180
            if self.fold == True:
                tipr[np.where(tipr > 90)] -= 180
                tipr[np.where(tipr < -90)] += 180

            # leave as the total unit circle 0 to 360
            elif self.fold == False:
                tipr %= 360
                tipr[np.where(tipr == 360.0)] = 0.0

            # make a dictionary of strikes with keys as period
            tiprdict = dict([(ff, jj) for ff, jj in zip(mt.period, tipr)])
            tip_list.append(tiprdict)

        #--> get min and max period
        self.max_per = np.amax([np.max(list(mm.keys())) for mm in inv_list], axis=0)
        self.min_per = np.amin([np.min(list(mm.keys())) for mm in pt_list], axis=0)

        # make empty arrays to put data into for easy manipulation
        medinv = np.zeros((nt, nc))
        medpt = np.zeros((nt, nc))
        medtipr = np.zeros((nt, nc))

        # make a list of periods from the longest period list
        self.period_arr = np.logspace(np.log10(self.min_per),
                                      np.log10(self.max_per),
                                      num=nt,
                                      base=10)
        self.period_dict = dict([(ii, jj) for jj, ii in 
                                 enumerate(self.period_arr)])

        # put data into arrays
        for ii, mm in enumerate(inv_list):
            mperiod = mm.keys()
            for jj, mp in enumerate(mperiod):
                for kk in self.period_dict.keys():
                    if mp > kk * (1 - self.period_tolerance) and \
                            mp < kk * (1 + self.period_tolerance):
                        ll = self.period_dict[kk]
                        medinv[ll, ii] = inv_list[ii][mp]
                        medpt[ll, ii] = pt_list[ii][mp]
                        medtipr[ll, ii] = tip_list[ii][mp]
                    else:
                        pass
        
        # make the arrays local variables
        self.med_inv = medinv
        self.med_pt = medpt
        self.med_tip = medtipr
        
        
    def get_mean(self, st_array):
        """
        get mean value
        """
        s_mean = 90 - np.mean(st_array[np.nonzero(st_array)])
        s_mean %= 360
        
        return s_mean
    
    def get_median(self, st_array):
        """
        get median value
        """
        s_median = 90 - np.median(st_array[np.nonzero(st_array)])
        s_median %= 360
        
        return s_median
        
    def get_mode(self, st_hist):
        """
        get mode from a historgram
        """
        s_mode = 90 - st_hist[1][np.where(st_hist[0] == st_hist[0].max())[0][0]]
        s_mode %= 360
        
        return s_mode
        
    def get_stats(self, st_array, st_hist, exponent=None):
        """
        print stats nicely
        """
        # print out the statistics of the strike angles
        s_mean = self.get_mean(st_array)
        s_median = self.get_median(st_array)
        s_mode = self.get_mode(st_hist)
        m = '-'*5
        if exponent is None:
            print('{0}All Periods{0}'.format(m))
        else:
            print('{0}Period Range {1:.3g} to {2:.3g} (s){0}'.format(m,
                                                                    10**exponent,
                                                                    10**(exponent + 1)))

        print('{0}median={1:.1f} mode={2:.1f} mean={3:.1f}'.format(' '*4,
                                                                   s_median,
                                                                   s_mode,
                                                                   s_mean))
        
        return s_median, s_mode, s_mean
    
    def get_plot_array(self, st_array):
        """
        get a plot array that has the min and max angles
        """
        
        # array goes from 
        st_array = st_array[np.nonzero(st_array)].flatten()
        st_array = st_array[np.isfinite(st_array)]
        plot_array = np.hstack([st_array, (st_array + 180) % 360])
        
        if self.plot_orthogonal:
            plot_array = np.hstack([plot_array, (plot_array + 90) % 360])
        
        if self.fold:
            plot_array %= 180
        
        return plot_array
    
    def plot(self, show=True):
        """
        plot Strike angles as rose plots
        
        """
        if not hasattr(self, 'med_inv'):
            self.make_strike_array()
            
        # font dictionary
        fd = {'size': self.font_size, 'weight': 'normal'}
        bw = self.bin_width
            
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = .07
        plt.rcParams['figure.subplot.right'] = .98
        plt.rcParams['figure.subplot.bottom'] = .09
        plt.rcParams['figure.subplot.top'] = .90
        plt.rcParams['figure.subplot.wspace'] = .2
        plt.rcParams['figure.subplot.hspace'] = .4

        if self.fold == True:
            histrange = (0, 180)
        elif self.fold == False:
            histrange = (0, 360)

        #-----Plot Histograms of the strike angles-----------------------------
        ### get the range in periods to plot
        if self.plot_range == 'data':
            self._bin_range = np.arange(np.floor(np.log10(self.min_per)),
                                        np.ceil(np.log10(self.max_per)), 1)
        else:
            self._bin_range = np.arange(np.floor(self.plot_range[0]),
                                        np.ceil(self.plot_range[1]), 1)

        #------------------plot indivdual decades------------------------------
        if self.plot_type == 1:
            nb = len(self._bin_range)
            # plot specs
            plt.rcParams['figure.subplot.hspace'] = .3
            plt.rcParams['figure.subplot.wspace'] = .3

            self.fig = plt.figure(self.fig_num, dpi=self.fig_dpi)
            plt.clf()
            
            n_subplots = 2
            if self.plot_tipper:
                n_subplots += 1
            for jj, bb in enumerate(self._bin_range, 1):
                # make subplots for invariants and phase tensor azimuths
                # dependent on vertical or horizontal orientation
                if 'h' in self.plot_orientation:
                    self.ax_inv = self.fig.add_subplot(n_subplots, nb, jj,
                                                        polar=True)
                    self.ax_pt = self.fig.add_subplot(n_subplots, nb, jj + nb,
                                                       polar=True)
                elif 'v' in self.plot_orientation:
                    self.ax_inv = self.fig.add_subplot(nb, n_subplots, 2*jj-1,
                                                       polar=True)
                    self.ax_pt = self.fig.add_subplot(nb, n_subplots, 2*jj,
                                                      polar=True)
                ax_list = [self.ax_inv, self.ax_pt]
                # vertical orientation
                if self.plot_tipper:
                    if 'h' in self.plot_orientation:
                        self.ax_tip = self.fig.add_subplot(n_subplots, nb,
                                                           jj + 2 * nb,
                                                           polar=True)
                    elif 'v' in self.plot_orientation:
                        self.ax_tip = self.fig.add_subplot(n_subplots, nb,
                                                           2*jj + 2,
                                                           polar=True)
                    ax_list.append(self.ax_tip)

                # make a list of indicies for each decades
                bin_list = []
                for ii, ff in enumerate(self.period_arr):
                    if ff > 10**bb and ff < 10**(bb + 1):
                        bin_list.append(ii)

                # extract just the subset for each decade
                plot_inv = self.get_plot_array(self.med_inv[bin_list, :])
                plot_pt = self.get_plot_array(self.med_pt[bin_list, :])
                
                if self.plot_tipper:
                    tr = self.get_plot_array(self.med_tip[bin_list, :])
                    # compute the historgram for the tipper strike
                    tr_hist = np.histogram(tr[np.nonzero(tr)].flatten(),
                                              bins=int(360/bw),
                                              range=histrange)

                    # make a bar graph with each bar being width of bw degrees
                    bar_tr = self.ax_tip.bar((tr_hist[1][:-1]) * np.pi / 180,
                                             tr_hist[0],
                                             width=bw * np.pi / 180)

                    # set color of the bars according to the number in that bin
                    # tipper goes from dark blue (low) to light blue (high)
                    for cc, bar in enumerate(bar_tr):
                        if self.color:
                            try:
                                fc = float(tr_hist[0][cc]) / tr_hist[0].max() * .9
                                bar.set_facecolor((.2, 1 - fc, .2))
                            except ZeroDivisionError:
                                pass
                        else:
                            bar.set_facecolor(self.color_tip)

                # estimate the histogram for the decade for invariants and pt
                inv_hist = np.histogram(plot_inv[np.nonzero(plot_inv)].flatten(),
                                        bins= int(360/bw),
                                        range=histrange)
                pt_hist = np.histogram(plot_pt,
                                       bins=int(360/bw),
                                       range=histrange)
        
                # plot the histograms
                self.bar_inv = self.ax_inv.bar((inv_hist[1][:-1]) * np.pi / 180,
                                                inv_hist[0],
                                                width=bw * np.pi / 180,
                                                zorder=10)

                self.bar_pt = self.ax_pt.bar((pt_hist[1][:-1]) * np.pi / 180,
                                              pt_hist[0],
                                              width=bw * np.pi / 180,
                                              zorder=10)

                # set the color of the bars according to the number in that bin
                # invariants go from purple (low) to red (high)
                for cc, bar in enumerate(self.bar_inv):
                    if self.color:
                        try:
                            fc = float(inv_hist[0][cc]) / inv_hist[0].max() * .8
                            bar.set_facecolor((.75, 1 - fc, 0))
                        except ZeroDivisionError:
                            pass
                    else:
                        bar.set_facecolor(self.color_inv)

                # pt goes from green (low) to orange (high)
                for cc, bar in enumerate(self.bar_pt):
                    if self.color:
                        try:
                            fc = float(pt_hist[0][cc]) / pt_hist[0].max() * .8
                            bar.set_facecolor((1-fc, 0, 1-fc))
                        except ZeroDivisionError:
                            pass
                    else:
                        bar.set_facecolor(self.color_pt)

                # make axis look correct with N to the top at 90.
                for aa, axh in enumerate(ax_list):
                    # set multiple locator to be every 15 degrees
                    axh.xaxis.set_major_locator(
                        MultipleLocator(30 * np.pi / 180))

                    ### set labels on the correct axis
                    axh.xaxis.set_ticklabels(['', '', '',
                                              '', '', '',
                                              '', '', '',
                                              '', '', ''])
                    ### set y limits if asked for
                    if self.ring_limits is not None:
                        axh.set_ylim(self.ring_limits)
                    axh.yaxis.set_major_locator(MultipleLocator(self.ring_spacing))

                    # make a light grid
                    axh.grid(alpha=.25, zorder=0)

                    # properties for the invariants
                    if aa == 0:
                        # limits need to be rotate 90 counter clockwise because
                        # we already rotated by 90 degrees so the range is
                        # from -90 to 270 with -90 being east
                        axh.set_xlim(-90 * np.pi / 180, 270 * np.pi / 180)

                        # label the plot with the mode value of strike
                        # need to subtract 90 again because the histogram is
                        # for ploting 0 east, 90 north measuring
                        # counter-clockwise
                        inv_median, inv_mode, inv_mean = self.get_stats(plot_inv, 
                                                                        inv_hist,
                                                                        bb)

                        ### place the estimated strike
                        axh.text(-np.pi/2, axh.get_ylim()[1] * self.text_pad,
                                 '{0:.1f}$^o$'.format(inv_mode),
                                 horizontalalignment='center',
                                 verticalalignment='baseline',
                                 fontdict={'size': self.text_size},
                                 bbox={'facecolor': self.color_inv,
                                       'alpha': .25})

                        #--> set title of subplot
                        if 'h' in self.plot_orientation:
                            axh.set_title(self.title_dict[bb], fontdict=fd,
                                          bbox={'facecolor': 'white', 'alpha': .25})

                            #--> set the title offset
                            axh.titleOffsetTrans._t = (0, .1)
                        elif 'v' in self.plot_orientation:
                            axh.set_ylabel(self.title_dict[bb], fontdict=fd,
                                           bbox={'facecolor': 'white', 'alpha': .25}, 
                                           rotation=0,
                                           labelpad=50)
                            axh.yaxis.set_label_position("right")

                    # set pt axes properties
                    elif aa == 1:
                        # limits go from -180 to 180 as that is how the angle
                        # is calculated
                        axh.set_xlim(-180 * np.pi / 180, 180 * np.pi / 180)
                        
                        pt_median, pt_mode, pt_mean = self.get_stats(plot_pt,
                                                                     pt_hist,
                                                                     bb)

                        ### put the estimated strike
                        axh.text(-np.pi/2, axh.get_ylim()[1] * self.text_pad,
                                 '{0:.1f}$^o$'.format(pt_mode),
                                 horizontalalignment='center',
                                 verticalalignment='baseline',
                                 fontdict={'size': self.text_size},
                                 bbox={'facecolor': self.color_pt, 
                                       'alpha': .25})

                    # set tipper axes properties
                    elif aa == 2:
                        # limits go from -180 to 180
                        axh.set_xlim(-180 * np.pi / 180, 180 * np.pi / 180)

                        tr_median, tr_mode, tr_mean = self.get_stats(tr,
                                                                     tr_hist,
                                                                     bb)
                        ### put the estimated strike 
                        axh.text(-np.pi/2, axh.get_ylim()[1] * self.text_pad,
                                 '{0:.1f}$^o$'.format(tr_mode),
                                 horizontalalignment='center',
                                 verticalalignment='baseline',
                                 fontdict={'size': self.text_size},
                                 bbox={'facecolor': self.color_tip, 
                                       'alpha': .25})
  
                    # set plot labels
                    if jj == 1 and 'h' in self.plot_orientation:
                        if aa == 0:
                            axh.set_ylabel('Strike (Z)', fontdict=fd,
                                           labelpad=self.font_size,
                                           bbox={'facecolor': self.color_inv,
                                                 'alpha': 0.25})
                        elif aa == 1:
                            axh.set_ylabel('PT Azimuth', fontdict=fd,
                                           labelpad=self.font_size,
                                           bbox={'facecolor': self.color_pt,
                                                 'alpha': .25})
                        elif aa == 2:
                            axh.set_ylabel('Tipper Strike', fd,
                                           labelpad=self.font_size,
                                           bbox={'facecolor': self.color_tip,
                                                 'alpha': 0.25})
                    # set plot labels
                    if jj == 1 and 'v' in self.plot_orientation:
                        if aa == 0:
                            axh.set_title('Strike (Z)', fontdict=fd,
                                           bbox={'facecolor': self.color_inv,
                                                 'alpha': 0.25})
                        elif aa == 1:
                            axh.set_title('PT Azimuth', fontdict=fd,
                                           bbox={'facecolor': self.color_pt,
                                                 'alpha': .25})
                        elif aa == 2:
                            axh.set_title('Tipper Strike', fd,
                                           bbox={'facecolor': self.color_tip,
                                                 'alpha': 0.25})

                    plt.setp(axh.yaxis.get_ticklabels(), visible=False)

            print('Note: North is assumed to be 0 and the strike angle is measured' +\
                  'clockwise positive.')

            if show:
                plt.show()

        #------------------Plot strike angles for all period ranges------------
        elif self.plot_type == 2:
            # plot specs
            plt.rcParams['figure.subplot.left'] = .07
            plt.rcParams['figure.subplot.right'] = .98
            plt.rcParams['figure.subplot.bottom'] = .100
            plt.rcParams['figure.subplot.top'] = .88
            plt.rcParams['figure.subplot.hspace'] = .3
            plt.rcParams['figure.subplot.wspace'] = .2

            self.fig = plt.figure(self.fig_num,
                                  self.fig_size,
                                  dpi=self.fig_dpi)
            plt.clf()
            # make subplots for invariants and phase tensor azimuths
            if not self.plot_tipper:
                self.ax_inv = self.fig.add_subplot(1, 2, 1, projection='polar')
                self.ax_pt = self.fig.add_subplot(1, 2, 2, projection='polar')
                ax_list = [self.ax_inv, self.ax_pt]
            else:
                self.ax_inv = self.fig.add_subplot(1, 3, 1, polar=True)
                self.ax_pt = self.fig.add_subplot(1, 3, 2, polar=True)
                self.ax_tip = self.fig.add_subplot(1, 3, 3, polar=True)
                ax_list = [self.ax_inv, self.ax_pt, self.ax_tip]

            # make a list of indicies for each decades
            bin_list = [self.period_dict[ff] for ff in self.period_arr
                       if ff > 10**self._bin_range.min() and 
                       ff < 10**self._bin_range.max()]

            # extract just the subset for each decade
            plot_inv = self.get_plot_array(self.med_inv[bin_list, :])
            plot_pt = self.get_plot_array(self.med_pt[bin_list, :])

            # estimate the histogram for the decade for invariants and pt
            inv_hist = np.histogram(plot_inv[np.nonzero(plot_inv)].flatten(),
                                   bins=int(360/bw),
                                   range=histrange)
                
            pt_hist = np.histogram(plot_pt,
                                  bins=int(360/bw),
                                  range=histrange)

            # plot the histograms
            self.bar_inv = self.ax_inv.bar((inv_hist[1][:-1]) * np.pi / 180,
                                          inv_hist[0],
                                          width=bw * np.pi / 180)

            self.bar_pt = self.ax_pt.bar((pt_hist[1][:-1]) * np.pi / 180,
                                        pt_hist[0],
                                        width=bw * np.pi / 180)

            # set color of invariants from purple (low) to red (high count)
            for cc, bar in enumerate(self.bar_inv):
                if self.color:
                    fc = float(inv_hist[0][cc]) / inv_hist[0].max() * .8
                    bar.set_facecolor((.75, 1-fc, 0))
                else:
                    bar.set_facecolor(self.color_inv)

            # set color of pt from green (low) to orange (high count)
            for cc, bar in enumerate(self.bar_pt):
                if self.color:
                    fc = float(pt_hist[0][cc]) / pt_hist[0].max() * .8
                    bar.set_facecolor((1-fc, 0, 1-fc))
                else:
                    bar.set_facecolor(self.color_pt)

            # plot tipper if desired
            if self.plot_tipper:
                tr = self.get_plot_array(self.med_tip[bin_list, :])

                tr_hist = np.histogram(tr[np.nonzero(tr)].flatten(),
                                      bins= int(360/bw),
                                      range=histrange)

                self.bar_tr = self.ax_tip.bar((tr_hist[1][:-1]) * np.pi / 180,
                                             tr_hist[0],
                                             width=bw * np.pi / 180)

                # set tipper color from dark blue (low) to light blue (high)
                for cc, bar in enumerate(self.bar_tr):
                    if self.color:
                        try:
                            fc = float(tr_hist[0][cc]) / tr_hist[0].max() * .8
                            bar.set_facecolor((.2, 1-fc, .2))
                        except ZeroDivisionError:
                            pass
                    else:
                        bar.set_facecolor(self.color_tip)
                    
            # make axis look correct with N to the top at 90.
            for aa, axh in enumerate(ax_list):
                # set major ticks to be every 30 degrees
                axh.xaxis.set_major_locator(MultipleLocator(2 * np.pi / 12))

                # set a light grid
                axh.grid(alpha=0.25)  # works in 2.0.2 not 2.1.0

                # set tick labels to be invisible
                plt.setp(axh.yaxis.get_ticklabels(), visible=False)  # works in 2.0.2 not 2.1.0

                # place the correct label at the cardinal directions
                axh.xaxis.set_ticklabels(['', '', '', '',
                                          '', '', '',
                                          '', '', '',
                                          '', '', ''])

                # set invariant axes propertiesz
                if aa == 0:
                    axh.set_ylim(0, inv_hist[0].max())
                    inv_median, inv_mode, inv_mean = self.get_stats(plot_inv,
                                                                    inv_hist)

                    axh.text(170 * np.pi / 180, axh.get_ylim()[1] * .65,
                             '{0:.1f}$^o$'.format(inv_mode),
                             horizontalalignment='center',
                             verticalalignment='baseline',
                             fontdict={'size': self.font_size},
                             bbox={'facecolor': self.color_inv, 'alpha': .25})

                    axh.set_title('Strike (Z)', fontdict=fd,
                                  bbox={'facecolor': self.color_inv,
                                        'alpha': 0.25})

                # set pt axes properties
                elif aa == 1:
                    axh.set_ylim(0, pt_hist[0].max())
                    pt_median, pt_mode, pt_mean = self.get_stats(plot_pt,
                                                                 pt_hist)

                    axh.text(170 * np.pi / 180, axh.get_ylim()[1] * .65,
                             '{0:.1f}$^o$'.format(pt_mode),
                             horizontalalignment='center',
                             verticalalignment='baseline',
                             fontdict={'size': self.text_size},
                             bbox={'facecolor': self.color_pt, 'alpha': 0.25})

                    axh.set_title('PT Azimuth', fontdict=fd,
                                  bbox={'facecolor': self.color_pt, 'alpha': 0.25})

                # set tipper axes properties
                elif aa == 2:
                    axh.set_ylim(0, tr_hist[0].max())
                    tr_median, tr_mode, tr_mean = self.get_stats(tr, tr_hist)

                    axh.text(170 * np.pi / 180, axh.get_ylim()[1] * .65,
                             '{0:.1f}$^o$'.format(tr_mode),
                             horizontalalignment='center',
                             verticalalignment='baseline',
                             fontdict={'size': self.text_size},
                             bbox={'facecolor': self.color_tip, 'alpha': 0.25})

                    axh.set_title('Tipper Strike', fontdict=fd,
                                  bbox={'facecolor': self.color_tip,
                                        'alpha': 0.25})

                # move title up a little to make room for labels
                axh.titleOffsetTrans._t = (0, .15)

            # remind the user what the assumptions of the strike angle are
            print('Note: North is assumed to be 0 and the strike angle is ' +\
                  'measured clockwise positive.')
            if show:
                plt.show()

    def save_plot(self, save_fn, file_format='pdf',
                  orientation='portrait', fig_dpi=None, close_plot='y'):
        """
        save_plot will save the figure to save_fn.

        Arguments
        -----------

            **save_fn: string
                          full path to save figure to, can be input as
                          - directory path -> the directory path to save to
                            in which the file will be saved as
                            save_fn/station_name_ResPhase.file_format

                          - full path -> file will be save to the given
                            path.  If you use this option then the format
                            will be assumed to be provided by the path

            **file_format: [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps...

            **orientation: [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait

            **fig_dpi: int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the dpi will be that at
                          which the figure was made.  I don't think that
                          it can be larger than dpi of the figure.

            **close_plot: [ y | n ]
                             -'y' will close the plot after saving.
                             -'n' will leave plot open

        Examples
        ---------

        :Example: ::

            >>> # to save plot as jpg
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotPhaseTensorMaps(edilist,freqspot=10)
            >>> p1.save_plot(r'/home/MT', file_format='jpg')
            'Figure saved to /home/MT/PTMaps/PTmap_phimin_10Hz.jpg'

        """
        # get rid of . in file format as it will be added later
        if file_format is not None:
            file_format = file_format.replace('.','')

        if fig_dpi is None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation)
            # plt.clf()
            # plt.close(self.fig)

        else:
            if not os.path.exists(save_fn):
                os.mkdir(save_fn)

            save_fn = os.path.join(save_fn, 'StrikeAnalysis.' + file_format)
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

        self.fig.clf()
        self.plot()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return "Plots phase tensor maps for one freq"

    def writeTextFiles(self, save_path=None):
        """
        Saves the strike information as a text file.
        """

        # check to see if the strikes have been calculated
        try:
            self.bin_width
        except AttributeError:
            self.plot()

        # get the path to save the file to
        if save_path is None:
            try:
                svpath = os.path.dirname(self.mt_list[0].fn)
            except TypeError:
                raise IOError('Need to input save_path, could not find path')

        else:
            svpath = save_path

        # set
        if self.fold == True:
            histrange = (-180, 180)

        elif self.fold == False:
            histrange = (0, 360)

        # set the bin width
        bw = self.bin_width

        slistinv = [['station']]
        slistpt = [['station']]
        slisttip = [['station']]

        # calculate the strikes for the different period bands
        for jj, bb in enumerate(self._bin_range):
            tstr = self.title_dict[bb].replace('$', '')
            tstr = tstr.replace('{', '').replace('}', '').replace('^', 'e')
            tstr = tstr.replace('s', '(s)')
            slistinv[0].append(tstr)
            slistpt[0].append(tstr)
            slisttip[0].append(tstr)

            # calculate the strike for the different period bands per station
            for kk, mt in enumerate(self.mt_list, 1):

                if jj == 0:
                    slistinv.append([mt.station])
                    slistpt.append([mt.station])
                    slisttip.append([mt.station])

                zinv = mt.Z.invariants
                pt = mt.pt
                tp = mt.Tipper

                bnlist = []
                for nn, per in enumerate(mt.period):
                    if per > 10**bb and per < 10**(bb + 1):
                        bnlist.append(nn)

                #---> strike from invariants
                zs = 90 - zinv.strike[bnlist]
                # fold so the angle goes from 0 to 180
                if self.fold == True:
                    # for plotting put the NW angles into the SE quadrant
                    zs[np.where(zs > 90)] = zs[np.where(zs > 90)] - 180
                    zs[np.where(zs < -90)] = zs[np.where(zs < -90)] + 180

                # leave as the total unit circle 0 to 360
                elif self.fold == False:
                    pass

                zshist = np.histogram(zs[np.nonzero(zs)].flatten(),
                                      bins=int(360/bw),
                                      range=histrange)

                #==============================================================
                # For putting the values into a useful text file
                # need to subtract 90 from the values to put them into
                # coordinates where north is 0 and east is 90, which is
                # different from plotting where east in the plotting function
                # is equal to 0 and north is 90, measuring counter-clockwise
                #==============================================================

                #==> compute mean
                invmean = 90 - zs.mean()
                if invmean < 0:
                    invmean += 360
                invmed = 90 - np.median(zs)

                #==> compute median
                if invmed < 0:
                    invmed += 360

                #==> compute mode
                invmode = 90 - zshist[1][np.where(
                    zshist[0] == zshist[0].max())[0][0]]
                if invmode < 0:
                    invmode += 360

                #==> append to list
                slistinv[kk].append((invmean,
                                     invmed,
                                     invmode))

                #---> strike from phase tensor
                az = pt.azimuth[0][bnlist]
                # fold so the angle goes from 0 to 180
                if self.fold == True:
                    az[np.where(az > 90)] = az[np.where(az > 90)] - 180
                    az[np.where(az < -90)] = az[np.where(az < -90)] + 180

                # leave as the total unit circle 0 to 360
                elif self.fold == False:
                    az[np.where(az < 0)] = az[np.where(az < 0)] + 360

                # == > compute mean
                ptmean1 = 90 - az.mean()
                if ptmean1 < 0:
                    ptmean1 += 360

                # == > compute median
                ptmed1 = 90 - np.median(az)
                if ptmed1 < 0:
                    ptmed1 += 360

                # == > compute mode
                azhist = np.histogram(az[np.nonzero(az)].flatten(),
                                      bins=int(360/bw),
                                      range=histrange)
                ptmode1 = 90 - azhist[1][np.where(
                    azhist[0] == azhist[0].max())[0][0]]
                if ptmode1 < 0:
                    ptmode1 += 360

                slistpt[kk].append((ptmean1,
                                    ptmed1,
                                    ptmode1))

                #---> strike from tipper
                # needs to be negative because measures clockwise
                if tp._Tipper.tipper is None:
                    tp._Tipper.tipper = np.zeros((len(mt.period), 1, 2),
                                                 dtype='complex')
                    tp.compute_components()

                tipr = -tp.angle_real[bnlist]

                # fold so the angle goes from 0 to 180
                if self.fold == True:
                    tipr[np.where(tipr > 90)] = tipr[np.where(tipr > 90)] - 180
                    tipr[np.where(tipr < -90)
                         ] = tipr[np.where(tipr < -90)] + 180

                # leave as the total unit circle 0 to 360
                elif self.fold == False:
                    tipr[np.where(tipr < 0)] = tipr[np.where(tipr < 0)] + 360

                tphist = np.histogram(tipr[np.nonzero(tipr)].flatten(),
                                      bins=int(360/bw),
                                      range=histrange)

                #==> compute mean
                tpmean1 = 90 - tipr.mean()
                if tpmean1 < 0:
                    tpmean1 += 360

                #==> compute median
                tpmed1 = 90 - np.median(tipr)
                if tpmed1 < 0:
                    tpmed1 += 360

                #==> compute mode
                tpmode1 = 90 - tphist[1][np.where(
                    tphist[0] == tphist[0].max())[0][0]]
                if tpmode1 < 0:
                    tpmode1 += 360

                #--> append statistics to list
                slisttip[kk].append((tpmean1,
                                     tpmed1,
                                     tpmode1))

            # make a list of indicies for each decades
            bin_list = []
            for ii, ff in enumerate(self._self.period_arr):
                if ff > 10**bb and ff < 10**(bb + 1):
                    bin_list.append(ii)

            # extract just the subset for each decade
            plot_inv = self._medinv[bin_list, :]
            plot_pt = self._medpt[bin_list, :]
            plot_pt = plot_pt[np.nonzero(plot_pt)].flatten()
            plot_pt = plot_pt[np.isfinite(plot_pt)]
            tr = self._medtp[bin_list, :]

            # estimate the histogram for the decade for invariants and pt
            inv_hist = np.histogram(plot_inv[np.nonzero(plot_inv)].flatten(),
                                   bins=int(360/bw),
                                   range=histrange)
            pt_hist = np.histogram(plot_pt,
                                  bins=int(360/bw),
                                  range=histrange)

            tr_hist = np.histogram(tr[np.nonzero(tr)].flatten(),
                                  bins=int(360/bw),
                                  range=histrange)

            #--> include the row for mean, median and mode for each parameter
            if jj == 0:
                slistinv.append(['mean'])
                slistinv.append(['median'])
                slistinv.append(['mode'])

                slistpt.append(['mean'])
                slistpt.append(['median'])
                slistpt.append(['mode'])

                slisttip.append(['mean'])
                slisttip.append(['median'])
                slisttip.append(['mode'])

            #--> compute mean, median and mode for invariants
            # == > mean
            imean = 90 - np.mean(plot_inv[np.nonzero(plot_inv)])
            if imean < 0:
                imean += 360

            # == > median
            imed = 90 - np.median(plot_inv[np.nonzero(plot_inv)])
            if imed < 0:
                imed += 360

            # == > mode
            imode = 90 - inv_hist[1][np.where(
                inv_hist[0] == inv_hist[0].max())[0][0]]
            if imode < 0:
                imode += 360

            #--> add them to the list of estimates
            slistinv[kk + 1].append(imean)
            slistinv[kk + 2].append(imed)
            slistinv[kk + 3].append(imode)

            #--> compute pt statistics
            # == > mean
            ptmean = 90 - np.mean(plot_pt)
            if ptmean < 0:
                ptmean = np.mean(plot_pt)

            # == > median
            ptmed = 90 - np.median(plot_pt)
            if ptmed < 0:
                ptmed += 360

            # == > mode
            ptmode = 90 - pt_hist[1][np.where(
                pt_hist[0] == pt_hist[0].max())[0][0]]
            if ptmode < 0:
                ptmode += 360

            #--> add the statistics to the parameter list
            slistpt[kk + 1].append(ptmean)
            slistpt[kk + 2].append(ptmed)
            slistpt[kk + 3].append(ptmode)

            #--> compute tipper statistics
            # == > mean
            tpmean = 90 - np.mean(tipr[np.nonzero(tipr)])
            if tpmean < 0:
                tpmean += 360

            # == > median
            tpmed = 90 - np.median(tipr[np.nonzero(tipr)])
            if tpmed < 0:
                tpmed += 360

            # == > mode
            tpmode = 90 - tr_hist[1][np.where(
                tr_hist[0] == tr_hist[0].max())[0][0]]
            if tpmode < 0:
                tpmode += 360

            #--> add the statistics to parameter list
            slisttip[kk + 1].append(tpmean)
            slisttip[kk + 2].append(tpmed)
            slisttip[kk + 3].append(tpmode)

        invfid = file(os.path.join(svpath, 'Strike.invariants'), 'w')
        ptfid = file(os.path.join(svpath, 'Strike.pt'), 'w')
        tpfid = file(os.path.join(svpath, 'Strike.tipper'), 'w')

        #---> write strike from the invariants
        # == > mean
        invfid.write('-' * 20 + 'MEAN' + '-' * 20 + '\n')
        for ii, l1 in enumerate(slistinv):
            for jj, l2 in enumerate(l1):
                if ii == 0:
                    invfid.write('{0:^16}'.format(l2))
                else:
                    if jj == 0:
                        invfid.write('{0:>16}'.format(l2 + ' ' * 6))
                    else:
                        try:
                            invfid.write('{0:^16}'.format(
                                         '{0: .2f}'.format(l2[0])))
                        except IndexError:
                            invfid.write('{0:^16}'.format(
                                         '{0: .2f}'.format(l2)))
            invfid.write('\n')

        # == > median
        invfid.write('-' * 20 + 'MEDIAN' + '-' * 20 + '\n')
        for ii, l1 in enumerate(slistinv):
            for jj, l2 in enumerate(l1):
                if ii == 0:
                    invfid.write('{0:^16}'.format(l2))
                else:
                    if jj == 0:
                        invfid.write('{0:>16}'.format(l2 + ' ' * 6))
                    else:
                        try:
                            invfid.write('{0:^16}'.format(
                                         '{0: .2f}'.format(l2[1])))
                        except IndexError:
                            invfid.write('{0:^16}'.format(
                                         '{0: .2f}'.format(l2)))
            invfid.write('\n')

        # == > mode
        invfid.write('-' * 20 + 'MODE' + '-' * 20 + '\n')
        for ii, l1 in enumerate(slistinv):
            for jj, l2 in enumerate(l1):
                if ii == 0:
                    invfid.write('{0:^16}'.format(l2))
                else:
                    if jj == 0:
                        invfid.write('{0:>16}'.format(l2 + ' ' * 6))
                    else:
                        try:
                            invfid.write('{0:^16}'.format(
                                         '{0: .2f}'.format(l2[2])))
                        except IndexError:
                            invfid.write('{0:^16}'.format(
                                         '{0: .2f}'.format(l2)))
            invfid.write('\n')
        invfid.close()

        #---> write the phase tensor text files
        ptfid.write('-' * 20 + 'MEAN' + '-' * 20 + '\n')
        for ii, l1 in enumerate(slistpt):
            for jj, l2 in enumerate(l1):
                if ii == 0:
                    ptfid.write('{0:^16}'.format(l2))
                else:
                    if jj == 0:
                        ptfid.write('{0:>16}'.format(l2 + ' ' * 6))
                    else:
                        try:
                            ptfid.write('{0:^16}'.format(
                                '{0: .2f}'.format(l2[0])))
                        except IndexError:
                            ptfid.write('{0:^16}'.format(
                                '{0: .2f}'.format(l2)))
            ptfid.write('\n')

        ptfid.write('-' * 20 + 'MEDIAN' + '-' * 20 + '\n')
        for ii, l1 in enumerate(slistpt):
            for jj, l2 in enumerate(l1):
                if ii == 0:
                    ptfid.write('{0:^16}'.format(l2))
                else:
                    if jj == 0:
                        ptfid.write('{0:>16}'.format(l2 + ' ' * 6))
                    else:
                        try:
                            ptfid.write('{0:^16}'.format(
                                '{0: .2f}'.format(l2[1])))
                        except IndexError:
                            ptfid.write('{0:^16}'.format(
                                '{0: .2f}'.format(l2)))
            ptfid.write('\n')

        ptfid.write('-' * 20 + 'MODE' + '-' * 20 + '\n')
        for ii, l1 in enumerate(slistpt):
            for jj, l2 in enumerate(l1):
                if ii == 0:
                    ptfid.write('{0:^16}'.format(l2))
                else:
                    if jj == 0:
                        ptfid.write('{0:>16}'.format(l2 + ' ' * 6))
                    else:
                        try:
                            ptfid.write('{0:^16}'.format(
                                '{0: .2f}'.format(l2[2])))
                        except IndexError:
                            ptfid.write('{0:^16}'.format(
                                '{0: .2f}'.format(l2)))
            ptfid.write('\n')
        ptfid.close()

        #---> write the tipper text files
        tpfid.write('-' * 20 + 'MEAN' + '-' * 20 + '\n')
        for ii, l1 in enumerate(slisttip):
            for jj, l2 in enumerate(l1):
                if ii == 0:
                    tpfid.write('{0:^16}'.format(l2))
                else:
                    if jj == 0:
                        tpfid.write('{0:>16}'.format(l2 + ' ' * 6))
                    else:
                        try:
                            tpfid.write('{0:^16}'.format(
                                '{0: .2f}'.format(l2[0])))
                        except IndexError:
                            tpfid.write('{0:^16}'.format(
                                '{0: .2f}'.format(l2)))
            tpfid.write('\n')

        tpfid.write('-' * 20 + 'MEDIAN' + '-' * 20 + '\n')
        for ii, l1 in enumerate(slisttip):
            for jj, l2 in enumerate(l1):
                if ii == 0:
                    tpfid.write('{0:^16}'.format(l2))
                else:
                    if jj == 0:
                        tpfid.write('{0:>16}'.format(l2 + ' ' * 6))
                    else:
                        try:
                            tpfid.write('{0:^16}'.format(
                                '{0: .2f}'.format(l2[1])))
                        except IndexError:
                            tpfid.write('{0:^16}'.format(
                                '{0: .2f}'.format(l2)))
            tpfid.write('\n')

        tpfid.write('-' * 20 + 'MODE' + '-' * 20 + '\n')
        for ii, l1 in enumerate(slisttip):
            for jj, l2 in enumerate(l1):
                if ii == 0:
                    tpfid.write('{0:^16}'.format(l2))
                else:
                    if jj == 0:
                        tpfid.write('{0:>16}'.format(l2 + ' ' * 6))
                    else:
                        try:
                            tpfid.write('{0:^16}'.format(
                                '{0: .2f}'.format(l2[2])))
                        except IndexError:
                            tpfid.write('{0:^16}'.format(
                                '{0: .2f}'.format(l2)))
            tpfid.write('\n')
        tpfid.close()
