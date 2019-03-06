# -*- coding: utf-8 -*-
"""
Created on Thu May 30 18:28:24 2013

@author: jpeacock-pr
"""

#==============================================================================

import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import MultipleLocator
import mtpy.imaging.mtplottools as mtpl
import mtpy.analysis.geometry as MTgy

#==============================================================================


class PlotStrike2D(object):
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

    Arguments:
    ----------
        **fn_list** : list of strings
                          full paths to .edi files to plot

        **z_object** : class mtpy.core.z.Z
                      object of mtpy.core.z.  If this is input be sure the
                      attribute z.freq is filled.  *default* is None

        **mt_object** : class mtpy.imaging.mtplot.MTplot
                        object of mtpy.imaging.mtplot.MTplot
                        *default* is None

        **fignum** : int
                     figure number to be plotted. *Default* is 1

        **fs** : float
                 font size for labels of plotting. *Default* is 10

        **dpi** : int
                  dots-per-inch resolution of figure, 300 is needed for
                  publications. *Default* is 300

        **thetar** : float
                     angle of rotation clockwise positive. *Default* is 0

        **ptol** : float
                   Tolerance level to match periods from different edi files.
                   *Default* is 0.05

        **text_dict** : dictionary
                  *'pad' : float
                           padding of the angle label at the bottom of each
                           polar diagram.  *Default* is 1.65

                  *'size' : float
                            font size


        **plot_range** : [ 'data' | (period_min,period_max) ]
                    period range to estimate the strike angle. Options are:
                        * *'data'* for estimating the strike for all periods
                            in the data.
                        * (pmin,pmax) for period min and period max, input as
                          (log10(pmin),log10(pmax))

        **plot_type** : [ 1 | 2 ]
                        -*1* to plot individual decades in one plot
                        -*2* to plot all period ranges into one polar diagram
                              for each strike angle estimation

        **plot_tipper** : [ 'y' | 'n' ]
                      -*'y'* to plot the tipper strike
                      -*'n'* to not plot tipper strike

        **pt_error_floor** : float
                    Maximum error in degrees that is allowed to estimate strike.
                    *Default* is None allowing all estimates to be used.

        **fold** : [ True | False ]
                    *True to plot only from 0 to 180
                    *False to plot from 0 to 360


    :Example: ::

        >>> import os
        >>> import mtpy.imaging.mtplot as mtplot
        >>> edipath = r"/home/EDIFiles"
        >>> edilist = [os.path.join(edipath,edi) for edi in os.listdir(edipath)
        >>> ...       if edi.find('.edi')>0]
        >>> #---plot rose plots in decades with tipper and an error floor on pt
        >>> strike = mtplot.plot_strike(fn_list=edilist, plot_type=1,\
                                        pt_error_floor=5)
        >>> #---plot all decades into one rose plot for each estimation---
        >>> strike.plot_type = 2
        >>> strike.redraw_plot()
        >>> #---save the plot---
        >>> strike.save_plot(r"/home/Figures")
        'Figure saved to /home/Figures/StrikeAnalysis_.pdf'

    Attributes:
    -----------

        -axhinv            matplotlib.axes instance for invariant strike
        -axhpt             matplotlib.axes instance for phase tensor strike
        -axhtip            matplotlib.axes instance for tipper strike

        -barinv            matplotlib.axes.bar instance for invariant strike
        -barpt             matplotlib.axes.bar instance for pt strike
        -bartr             matplotlib.axes.bar instance for tipper strike

        -bin_width         width of histogram bins in degrees

        -fig               matplotlib.figure instance of plot
        -fig_dpi           dots-per-inch resolution of figure
        -fig_num           number of figure being plotted
        -fig_size          size of figure in inches
        -fold              boolean to fold angles to range from [0,180] or
                           [0,360]

        -font_size         font size of axes tick labels

        -mt_list            list of mtplot.MTplot instances containing all
                           the important information for each station
        -period_tolerance  tolerance to look for periods being plotted

        -plot_range        range of periods to plot
        -plot_tipper       string to tell program to plot induction arrows
        -plot_type         string to tell program how to plot strike angles
        -plot_yn           plot strike on instance creation
        -pt_error_floor    error floor to plot phase tensor strike, anything
                           above this error will not be plotted
        -text_pad          padding between text and rose diagram
        -text_size         font size of text labeling the mode of the histogram
        -title_dict        title dictionary

    Methods:
    --------

        -plot                 plots the pseudo section
        -redraw_plot          on call redraws the plot from scratch
        -save_figure          saves figure to a file of given format
        -update_plot          updates the plot while still active
        -writeTextFiles       writes parameters of the phase tensor and tipper
                              to text files.

    """

    def __init__(self, **kwargs):

        fn_list = kwargs.pop('fn_list', None)
        z_object_list = kwargs.pop('z_object_list', None)
        tipper_object_list = kwargs.pop('tipper_object_list', None)
        mt_object_list = kwargs.pop('mt_object_list', None)

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
        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_dpi = kwargs.pop('fig_dpi', 300)
        self.fig_size = kwargs.pop('fig_size', [7, 5])

        self.plot_num = kwargs.pop('plot_num', 1)
        self.plot_type = kwargs.pop('plot_type', 2)
        self.plot_title = kwargs.pop('plot_title', None)
        self.plot_range = kwargs.pop('plot_range', 'data')
        self.plot_tipper = kwargs.pop('plot_tipper', 'n')

        self.period_tolerance = kwargs.pop('period_tolerance', .05)
        self.pt_error_floor = kwargs.pop('pt_error_floor', None)
        self.fold = kwargs.pop('fold', True)
        self.bin_width = kwargs.pop('bin_width', 5)
        self.skew_threshold = kwargs.pop('skew_threshold', 3)

        self.font_size = kwargs.pop('font_size', 7)

        text_dict = kwargs.pop('text_dict', {})
        try:
            self.text_pad = text_dict['pad']
        except KeyError:
            self.text_pad = 0.6

        try:
            self.text_size = text_dict['size']
        except KeyError:
            self.text_size = self.font_size

        # make a dictionary for plotting titles
        self.title_dict = {}
        self.title_dict[-5] = '10$^{-5}$--10$^{-4}$s'
        self.title_dict[-4] = '10$^{-4}$--10$^{-3}$s'
        self.title_dict[-3] = '10$^{-3}$--10$^{-2}$s'
        self.title_dict[-2] = '10$^{-2}$--10$^{-1}$s'
        self.title_dict[-1] = '10$^{-1}$--10$^{0}$s'
        self.title_dict[0] = '10$^{0}$--10$^{1}$s'
        self.title_dict[1] = '10$^{1}$--10$^{2}$s'
        self.title_dict[2] = '10$^{2}$--10$^{3}$s'
        self.title_dict[3] = '10$^{3}$--10$^{4}$s'
        self.title_dict[4] = '10$^{4}$--10$^{5}$s'
        self.title_dict[5] = '10$^{5}$--10$^{6}$s'

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

        for ii, mt in enumerate(self.mt_list):
            mt.rotation_angle = self._rot_z[ii]

    def _get_rot_z(self):
        return self._rot_z

    rot_z = property(fget=_get_rot_z, fset=_set_rot_z,
                     doc="""rotation angle(s)""")

    def plot(self):

        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = .07
        plt.rcParams['figure.subplot.right'] = .98
        plt.rcParams['figure.subplot.bottom'] = .09
        plt.rcParams['figure.subplot.top'] = .90
        plt.rcParams['figure.subplot.wspace'] = .2
        plt.rcParams['figure.subplot.hspace'] = .4

        bw = self.bin_width
        histrange = (0, 360)

        # set empty lists that will hold dictionaries with keys as the period
        ptlist = []
        tiprlist = []

        # initialize some parameters
        nc = len(self.mt_list)
        nt = 0
        kk = 0

        for dd, mt in enumerate(self.mt_list):

            #--> set the period
            period = mt.period

            # get maximum length of periods
            if len(period) > nt:
                nt = len(period)

            # estimate where only the 2D sections are
            dim_2d = MTgy.dimensionality(z_object=mt._Z,
                                         skew_threshold=self.skew_threshold)
            index_2d = np.where(dim_2d == 2)[0]
            #------------get strike from phase tensor strike angle-------------
            pt = mt.pt
            az = (90 - pt.azimuth[index_2d]) % 360
            az_err = pt.azimuth_err[index_2d]

            # need to add 90 because pt assumes 0 is north and
            # negative because measures clockwise.

            # put an error max on the estimation of strike angle
            if self.pt_error_floor:
                az[np.where(az_err > self.pt_error_floor)] = 0.0

            # make a dictionary of strikes with keys as period
            mdictpt = dict([(ff, jj)
                            for ff, jj in zip(mt.period[index_2d], az)])
            ptlist.append(mdictpt)

            #-----------get tipper strike------------------------------------
            tip = mt.Tipper
            if tip.tipper is None:
                tip.tipper = np.zeros((len(mt.period), 1, 2),
                                              dtype='complex')
                tip.compute_components()

            # needs to be negative because measures clockwise
            tipr = -tip.angle_real[index_2d]

            tipr[np.where(tipr == 180.)] = 0.0
            tipr[np.where(tipr == -180.)] = 0.0

            # make sure the angle is between 0 and 360
            tipr = tipr % 360

            # make a dictionary of strikes with keys as period
            tiprdict = dict([(ff, jj)
                             for ff, jj in zip(mt.period[index_2d], tipr)])
            tiprlist.append(tiprdict)

        #--> get min and max period
        maxper = np.max([np.max(list(mm.keys())) for mm in ptlist if list(mm.keys())])
        minper = np.min([np.min(list(mm.keys())) for mm in ptlist if list(mm.keys())])

        # make empty arrays to put data into for easy manipulation
        medpt = np.zeros((nt, nc))
        medtipr = np.zeros((nt, nc))

        # make a list of periods from the longest period list
        plist = np.logspace(
            np.log10(minper),
            np.log10(maxper),
            num=nt,
            base=10)
        pdict = dict([(ii, jj) for jj, ii in enumerate(plist)])

        self._plist = plist

        # put data into arrays
        for ii, mm in enumerate(ptlist):
            mperiod = list(mm.keys())
            for jj, mp in enumerate(mperiod):
                for kk in list(pdict.keys()):
                    if mp > kk * (1 - self.period_tolerance) and \
                            mp < kk * (1 + self.period_tolerance):
                        ll = pdict[kk]
                        medpt[ll, ii] = ptlist[ii][mp]
                        medtipr[ll, ii] = tiprlist[ii][mp]
                    else:
                        pass

        # make the arrays local variables
        self._medpt = medpt
        self._medtp = medtipr

        #-----Plot Histograms of the strike angles-----------------------------
        if self.plot_range == 'data':
            brange = np.arange(np.floor(np.log10(minper)),
                               np.ceil(np.log10(maxper)), 1)
        else:
            brange = np.arange(np.floor(self.plot_range[0]),
                               np.ceil(self.plot_range[1]), 1)

        self._brange = brange

        # font dictionary
        fd = {'size': self.font_size, 'weight': 'normal'}

        #------------------plot indivdual decades------------------------------
        if self.plot_type == 1:
            # plot specs
            plt.rcParams['figure.subplot.hspace'] = .3
            plt.rcParams['figure.subplot.wspace'] = .3

            self.fig = plt.figure(self.fig_num, dpi=self.fig_dpi)
            plt.clf()
            nb = len(brange)
            for jj, bb in enumerate(brange, 1):
                # make subplots for invariants and phase tensor azimuths
                if self.plot_tipper == 'n':
                    self.axhpt = self.fig.add_subplot(1, nb, jj, polar=True)
                    axlist = [self.axhpt]

                if self.plot_tipper == 'y':
                    self.axhpt = self.fig.add_subplot(2, nb, jj, polar=True)
                    self.axhtip = self.fig.add_subplot(2, nb, jj + nb,
                                                       polar=True)
                    axlist = [self.axhpt, self.axhtip]

                # make a list of indicies for each decades
                binlist = []
                for ii, ff in enumerate(plist):
                    if ff > 10**bb and ff < 10**(bb + 1):
                        binlist.append(ii)

                # extract just the subset for each decade
                gg = medpt[binlist, :]
                if self.plot_tipper == 'y':
                    tr = medtipr[binlist, :]

                    # compute the historgram for the tipper strike
                    trhist = np.histogram(tr[np.nonzero(tr)].flatten(),
                                          bins=int(360/bw),
                                          range=histrange)

                    # make a bar graph with each bar being width of bw degrees
                    bartr = self.axhtip.bar((trhist[1][:-1]) * np.pi / 180,
                                            trhist[0],
                                            width=bw * np.pi / 180)

                    # set color of the bars according to the number in that bin
                    # tipper goes from dark blue (low) to light blue (high)
                    for cc, bar in enumerate(bartr):
                        try:
                            fc = float(trhist[0][cc]) / trhist[0].max() * .9
                        except ZeroDivisionError:
                            fc = 1.0
                        bar.set_facecolor((0, 1 - fc / 2, fc))

                # estimate the histogram for the decade for invariants and pt
                pthist = np.histogram(gg[np.nonzero(gg)].flatten(),
                                      bins=int(360/bw),
                                      range=histrange)

                # plot the histograms
                self.barpt = self.axhpt.bar((pthist[1][:-1]) * np.pi / 180,
                                            pthist[0],
                                            width=bw * np.pi / 180)

                # set the color of the bars according to the number in that bin
                # pt goes from green (low) to orange (high)
                for cc, bar in enumerate(self.barpt):
                    try:
                        fc = float(pthist[0][cc]) / pthist[0].max() * .8
                    except ZeroDivisionError:
                        fc = 1.0
                    bar.set_facecolor((fc, 1 - fc, 0))

                # make axis look correct with N to the top at 90.
                for aa, axh in enumerate(axlist):
                    # set multiple locator to be every 15 degrees
                    axh.xaxis.set_major_locator(
                        MultipleLocator(30 * np.pi / 180))

                    # set labels on the correct axis
                    axh.xaxis.set_ticklabels(['', 'E', '', '',
                                              'N', '', '',
                                              'W', '', '',
                                              'S', '', ''])
                    # make a light grid
                    axh.grid(alpha=.25)

                    # set pt axes properties
                    if aa == 0:
                        # limits go from -180 to 180 as that is how the angle
                        # is calculated
                        axh.set_xlim(0, 2 * np.pi)

                        # label plot with the mode of the strike angle
                        ptmode = (90 - pthist[1][np.where(
                                  pthist[0] == pthist[0].max())[0][0]]) % 360

                        ptmedian = (90 - np.median(gg[np.nonzero(gg)])) % 360

                        ptmean = (90 - np.mean(gg[np.nonzero(gg)])) % 360

                        axh.text(np.pi, axh.get_ylim()[1] * self.text_pad,
                                 '{0:.1f}$^o$'.format(ptmode),
                                 horizontalalignment='center',
                                 verticalalignment='baseline',
                                 fontdict={'size': self.text_size},
                                 bbox={'facecolor': (.9, .9, 0), 'alpha': .25})

                        # print out the results for the strike angles
                        print('-----Period Range {0:.3g} to {1:.3g} (s)-----'.format(10**bb,
                                                                                     10**(bb + 1)))
                        print('   *PT Strike:     median={0:.1f} mode={1:.1f} mean={2:.1f}'.format(
                            ptmedian,
                            ptmode,
                            ptmean))

                        if self.plot_tipper != 'y':
                            print('\n')

                         #--> set title of subplot
                        axh.set_title(self.title_dict[bb], fontdict=fd,
                                      bbox={'facecolor': 'white', 'alpha': .25})

                        #--> set the title offset
                        axh.titleOffsetTrans._t = (0, .1)

                    # set tipper axes properties
                    elif aa == 1:
                        # limits go from -180 to 180
                        axh.set_xlim(0, 2 * np.pi)

                        # label plot with mode
                        tpmode = (90 - trhist[1][np.where(
                                  trhist[0] == trhist[0].max())[0][0]]) % 360

                        tpmedian = (90 - np.median(tr[np.nonzero(tr)])) % 360

                        tpmean = (90 - np.mean(tr[np.nonzero(tr)])) % 360

                        axh.text(np.pi, axh.get_ylim()[1] * self.text_pad,
                                 '{0:.1f}$^o$'.format(tpmode),
                                 horizontalalignment='center',
                                 verticalalignment='baseline',
                                 fontdict={'size': self.text_size},
                                 bbox={'facecolor': (0, .1, .9), 'alpha': .25})

                        # print out statistics for strike angle
                        print('   *Tipper Strike: median={0:.1f} mode={1:.1f} mean={2:.1f}'.format(
                            tpmedian,
                            tpmode,
                            tpmode))
                        print('\n')
                        if nb > 5:
                            axh.set_title(self.title_dict[bb], fontdict=fd,
                                          bbox={'facecolor': 'white', 'alpha': .25})

                    # set plot labels
                    if jj == 1:
                        if aa == 0:
                            axh.set_ylabel('PT Azimuth', fontdict=fd,
                                           labelpad=self.font_size,
                                           bbox={'facecolor': (.9, .9, 0),
                                                 'alpha': .25})
                        elif aa == 1:
                            axh.set_ylabel('Tipper Strike', fd,
                                           labelpad=self.font_size,
                                           bbox={'facecolor': (0, .1, .9),
                                                 'alpha': 0.25})

                    plt.setp(axh.yaxis.get_ticklabels(), visible=False)

            print('Note: North is assumed to be 0 and the strike angle is measured' +\
                  'clockwise positive.')

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
            if self.plot_tipper == 'n':
                self.axhpt = self.fig.add_subplot(1, 1, 1, polar=True)
                axlist = [self.axhpt]
            else:
                self.axhpt = self.fig.add_subplot(1, 2, 1, polar=True)
                self.axhtip = self.fig.add_subplot(1, 2, 2, polar=True)
                axlist = [self.axhpt, self.axhtip]

            # make a list of indicies for each decades
            binlist = [pdict[ff] for ff in plist
                       if ff > 10**brange.min() and ff < 10**brange.max()]

            # extract just the subset for each decade
            gg = medpt[binlist, :]

            # estimate the histogram for the decade for invariants and pt
            pthist = np.histogram(gg[np.nonzero(gg)].flatten(),
                                  bins=int(360/bw),
                                  range=histrange)

            # plot the histograms
            self.barpt = self.axhpt.bar((pthist[1][:-1]) * np.pi / 180,
                                        pthist[0],
                                        width=bw * np.pi / 180)

            # set color of pt from green (low) to orange (high count)
            for cc, bar in enumerate(self.barpt):
                fc = float(pthist[0][cc]) / pthist[0].max() * .8
                bar.set_facecolor((fc, 1 - fc, 0))

            # plot tipper if desired
            if self.plot_tipper == 'y':
                tr = self._medtp[binlist, :]

                trhist = np.histogram(tr[np.nonzero(tr)].flatten(),
                                      bins=int(360/bw),
                                      range=histrange)

                self.bartr = self.axhtip.bar((trhist[1][:-1]) * np.pi / 180,
                                             trhist[0],
                                             width=bw * np.pi / 180)

                # set tipper color from dark blue (low) to light blue (high)
                for cc, bar in enumerate(self.bartr):
                    try:
                        fc = float(trhist[0][cc]) / trhist[0].max() * .9
                        bar.set_facecolor((0, 1 - fc / 2, fc))
                    except ZeroDivisionError:
                        pass

            # make axis look correct with N to the top at 90.
            for aa, axh in enumerate(axlist):
                # set major ticks to be every 30 degrees
                axh.xaxis.set_major_locator(MultipleLocator(2 * np.pi / 12))

                # set a light grid
                axh.grid(alpha=0.25)

                # set tick labels to be invisible
                plt.setp(axh.yaxis.get_ticklabels(), visible=False)

                # place the correct label at the cardinal directions
                axh.xaxis.set_ticklabels(['', 'E', '', '',
                                          'N', '', '',
                                          'W', '', '',
                                          'S', '', ''])

                # set pt axes properties
                if aa == 0:
                    axh.set_ylim(0, pthist[0].max())

                    ptmode = (90 - pthist[1][np.where(
                              pthist[0] == pthist[0].max())[0][0]]) % 360

                    ptmedian = (90 - np.median(gg[np.nonzero(gg)])) % 360

                    ptmean = (90 - np.mean(gg[np.nonzero(gg)])) % 360

                    axh.text(170 * np.pi / 180, axh.get_ylim()[1] * .65,
                             '{0:.1f}$^o$'.format(ptmode),
                             horizontalalignment='center',
                             verticalalignment='baseline',
                             fontdict={'size': self.text_size},
                             bbox={'facecolor': (.9, .9, 0), 'alpha': 0.25})

                    # print results of strike analysis for pt
                    print('-----Period Range {0:.3g} to {1:.3g} (s)-----'.format(10**brange[0],
                                                                                 10**brange[-1]))
                    print('   *PT Strike:     median={0:.1f} mode={1:.1f} mean={2:.1f}'.format(
                        ptmedian,
                        ptmode,
                        ptmean))

                    if self.plot_tipper != 'y':
                        print('\n')

                    axh.set_title('PT Azimuth', fontdict=fd,
                                  bbox={'facecolor': (.9, .9, 0), 'alpha': 0.25})

                # set tipper axes properties
                elif aa == 2:
                    axh.set_ylim(0, trhist[0].max())

                    tpmode = (90 - trhist[1][np.where(
                              trhist[0] == trhist[0].max())[0][0]]) % 360

                    tpmedian = (90 - np.median(tr[np.nonzero(tr)])) % 360

                    tpmean = (90 - np.mean(tr[np.nonzero(tr)])) % 360

                    axh.text(170 * np.pi / 180, axh.get_ylim()[1] * .65,
                             '{0:.1f}$^o$'.format(tpmode),
                             horizontalalignment='center',
                             verticalalignment='baseline',
                             fontdict={'size': self.text_size},
                             bbox={'facecolor': (0, .1, .9), 'alpha': 0.25})

                    print('   *Tipper Stike:  median={0:.1f} mode={1:.1f} mean={2:.1f}\n'.format(
                        tpmedian,
                        tpmode,
                        tpmean))

                    axh.set_title('Tipper Strike', fontdict=fd,
                                  bbox={'facecolor': (0, .1, .9), 'alpha': 0.25})

                # move title up a little to make room for labels
                axh.titleOffsetTrans._t = (0, .15)

            # remind the user what the assumptions of the strike angle are
            print('Note: North is assumed to be 0 and the strike angle is ' +\
                  'measured clockwise positive.')

            plt.show()

    def save_plot(self, save_fn, file_format='pdf',
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

            >>> # to save plot as jpg
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotPhaseTensorMaps(edilist,freqspot=10)
            >>> p1.save_plot(r'/home/MT', file_format='jpg')
            'Figure saved to /home/MT/PTMaps/PTmap_phimin_10Hz.jpg'

        """

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

            save_fn = os.path.join(save_fn, 'StrikeAnalysis_' + file_format)
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
        for jj, bb in enumerate(self._brange):
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
            binlist = []
            for ii, ff in enumerate(self._plist):
                if ff > 10**bb and ff < 10**(bb + 1):
                    binlist.append(ii)

            # extract just the subset for each decade
            hh = self._medinv[binlist, :]
            gg = self._medpt[binlist, :]
            tr = self._medtp[binlist, :]

            # estimate the histogram for the decade for invariants and pt
            invhist = np.histogram(hh[np.nonzero(hh)].flatten(),
                                   bins=int(360/bw),
                                   range=histrange)
            pthist = np.histogram(gg[np.nonzero(gg)].flatten(),
                                  bins=int(360/bw),
                                  range=histrange)

            trhist = np.histogram(tr[np.nonzero(tr)].flatten(),
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
            imean = 90 - np.mean(hh[np.nonzero(hh)])
            if imean < 0:
                imean += 360

            # == > median
            imed = 90 - np.median(hh[np.nonzero(hh)])
            if imed < 0:
                imed += 360

            # == > mode
            imode = 90 - invhist[1][np.where(
                invhist[0] == invhist[0].max())[0][0]]
            if imode < 0:
                imode += 360

            #--> add them to the list of estimates
            slistinv[kk + 1].append(imean)
            slistinv[kk + 2].append(imed)
            slistinv[kk + 3].append(imode)

            #--> compute pt statistics
            # == > mean
            ptmean = 90 - np.mean(gg[np.nonzero(gg)])
            if ptmean < 0:
                ptmean = np.mean(gg[np.nonzero(gg)])

            # == > median
            ptmed = 90 - np.median(gg[np.nonzero(gg)])
            if ptmed < 0:
                ptmed += 360

            # == > mode
            ptmode = 90 - pthist[1][np.where(
                pthist[0] == pthist[0].max())[0][0]]
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
            tpmode = 90 - trhist[1][np.where(
                trhist[0] == trhist[0].max())[0][0]]
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
