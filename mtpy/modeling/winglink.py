# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 15:19:30 2011

deal with output files from winglink.

@author: jp
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
import matplotlib.colorbar as mcb
from matplotlib.colors import Normalize

#------------------------------------------------------------------------------


class WLInputError(Exception):
    pass


def read_output_file(output_fn):
    """
    Reads in an output file from winglink and returns the data
    in the form of a dictionary of structured arrays.

    Arguments:
    -----------
        **output_fn** : string
                        the full path to winglink outputfile

    Returns:
    ----------
        **wl_data** : dictionary with keys of station names
                      each station contains a structured array with keys
                      * 'station' --> station name
                      * 'period' --> periods to plot
                      * 'te_res' --> TE resistivity in linear scale
                      * 'tm_res' --> TM resistivity in linear scale
                      * 'te_phase' --> TE phase in deg
                      * 'tm_phase' --> TM phase in deg
                      * 're_tip' --> real tipper amplitude.
                      * 'im_tip' --> imaginary tipper amplitude
                      * 'rms' --> RMS for the station
                      * 'index' --> order from left to right of station number 

        .. note:: each data is an np.ndarray(2, num_periods) where the first
                  index is the data and the second index is the model response
    """
    if os.path.isfile(output_fn) is False:
        raise WLInputError('Cannot find {0}, check path'.format(output_fn))
    ofid = open(output_fn, 'r')
    lines = ofid.readlines()

    idict = {}
    stationlst = []

    # get title line
    titleline = lines[1].replace('"', '')
    titleline = titleline.rstrip().split(',')
    title = titleline[1].split(':')[1]
    profile = titleline[0].split(':')[1]
    inversiontype = lines[2].rstrip()

    dkeys = ['obs_tm_res', 'obs_tm_phase', 'mod_tm_res', 'mod_tm_phase',
             'obs_te_res', 'obs_te_phase', 'mod_te_res', 'mod_te_phase',
             'obs_re_tip', 'obs_im_tip', 'mod_re_tip', 'mod_im_tip', 'period']

    index = 0
    for line in lines[3:]:
        # get the beginning of the station block
        if line.find('Data for station') == 0:
            station = line.rstrip().split(':')[1][1:]
            idict[station] = {}
            stationlst.append(station)
            print('Reading in station: ', station)
            idict[station]['index'] = index
            index += 1
            for key in dkeys:
                idict[station][key] = []
        # get rms
        elif line.find('RMS') == 0:
            idict[station]['rms'] = float(line.strip().split(' = ')[1])
        # skip the divding line
        elif line.find('==') == 0:
            pass
        # get the data
        else:
            linelst = line.split()
            if len(linelst) == len(dkeys):
                for kk, key in enumerate(dkeys):
                    try:
                        if key.find('phase') >= 0:
                            idict[station][key].append(-1 * float(linelst[kk]))
                        else:
                            idict[station][key].append(float(linelst[kk]))
                    except ValueError:
                        idict[station][key].append(0)
            else:
                pass

    # get data into a more useful format that takes into account any masking of
    # data points.

    data = {}
    for st in list(idict.keys()):
        data[st] = {}
        data[st]['station'] = st
        data[st]['index'] = int(idict[st]['index'])
        data[st]['period'] = np.array(idict[st]['period'])
        data[st]['te_res'] = np.array([np.array(idict[st]['obs_te_res']),
                                       np.array(idict[st]['mod_te_res'])])
        data[st]['tm_res'] = np.array([np.array(idict[st]['obs_tm_res']),
                                       np.array(idict[st]['mod_tm_res'])])
        data[st]['te_phase'] = np.array([np.array(idict[st]['obs_te_phase']),
                                         np.array(idict[st]['mod_te_phase'])])
        data[st]['tm_phase'] = np.array([np.array(idict[st]['obs_tm_phase']),
                                         np.array(idict[st]['mod_tm_phase'])])
        data[st]['re_tip'] = np.array([np.array(idict[st]['obs_re_tip']),
                                       np.array(idict[st]['mod_re_tip'])])
        data[st]['im_tip'] = np.array([np.array(idict[st]['obs_im_tip']),
                                       np.array(idict[st]['mod_im_tip'])])
        data[st]['rms'] = float(idict[st]['rms'])

    return data


#------------------------------------------------------------------------------
def read_model_file(model_fn):
    """
    readModelFile reads in the XYZ txt file output by Winglink.    

    Inputs:
        modelfile = fullpath and filename to modelfile
        profiledirection = 'ew' for east-west predominantly, 'ns' for 
                            predominantly north-south.  This gives column to 
                            fix
    """

    mfid = open(model_fn, 'r')
    lines = mfid.readlines()
    nlines = len(lines)

    X = np.zeros(nlines)
    Y = np.zeros(nlines)
    Z = np.zeros(nlines)
    rho = np.zeros(nlines)
    # file starts from the bottom of the model grid in X Y Z Rho coordinates
    for ii, line in enumerate(lines):
        linestr = line.split()
        X[ii] = float(linestr[0])
        Y[ii] = float(linestr[1])
        Z[ii] = float(linestr[2])
        rho[ii] = float(linestr[3])

    X = X[np.nonzero(X)]
    Y = Y[np.nonzero[Y]]
    Z = Z[np.nonzero[Z]]
    rho = rho[np.nonzero[rho]]
    return X, Y, Z, rho


#==============================================================================
# plot the MT and model responses
#==============================================================================
class PlotResponse():
    """
    Helper class to deal with plotting the MT response and occam2d model.

    Arguments:
    -------------
        **data_fn** : string
                      full path to data file

        **resp_fn** : string or list
                      full path(s) to response file(s)   


    ==================== ======================================================
    Attributes/key words            description
    ==================== ======================================================
    ax_list              list of matplotlib.axes instances for use with
                         OccamPointPicker    
    color_mode           [ 'color' | 'bw' ] plot figures in color or 
                         black and white ('bw')
    cted                 color of Data TE marker and line
    ctem                 color of Model TE marker and line
    ctewl                color of Winglink Model TE marker and line
    ctmd                 color of Data TM marker and line
    ctmm                 color of Model TM marker and line
    ctmwl                color of Winglink Model TM marker and line
    e_capsize            size of error bar caps in points
    e_capthick           line thickness of error bar caps in points
    err_list             list of line properties of error bars for use with
                         OccamPointPicker
    fig_dpi              figure resolution in dots-per-inch 
    fig_list             list of dictionaries with key words
                         station --> station name
                         fig --> matplotlib.figure instance
                         axrte --> matplotlib.axes instance for TE app.res
                         axrtm --> matplotlib.axes instance for TM app.res
                         axpte --> matplotlib.axes instance for TE phase
                         axptm --> matplotlib.axes instance for TM phase

    fig_num              starting number of figure
    fig_size             size of figure in inches (width, height)
    font_size            size of axes ticklabel font in points
    line_list            list of matplotlib.Line instances for use with 
                         OccamPointPicker
    lw                   line width of lines in points
    ms                   marker size in points
    mted                 marker for Data TE mode
    mtem                 marker for Model TE mode
    mtewl                marker for Winglink Model TE
    mtmd                 marker for Data TM mode
    mtmm                 marker for Model TM mode
    mtmwl                marker for Winglink TM mode
    period               np.ndarray of periods to plot 
    phase_limits         limits on phase plots in degrees (min, max)
    plot_model_error     [ 'y' | 'n' ] *default* is 'y' to plot model errors 
    plot_num             [ 1 | 2 ] 
                         1 to plot both modes in a single plot
                         2 to plot modes in separate plots (default)
    plot_tipper          [ 'y' | 'n' ] plot tipper data if desired
    plot_type            [ '1' | station_list]
                         '1' --> to plot all stations in different figures
                         station_list --> to plot a few stations, give names
                         of stations ex. ['mt01', 'mt07']
    plot_yn              [ 'y' | 'n']
                         'y' --> to plot on instantiation
                         'n' --> to not plot on instantiation
    res_limits           limits on resistivity plot in log scale (min, max)
    rp_list               list of dictionaries from read2Ddata
    station_list          station_list list of stations in rp_list
    subplot_bottom       subplot spacing from bottom (relative coordinates) 
    subplot_hspace       vertical spacing between subplots
    subplot_left         subplot spacing from left  
    subplot_right        subplot spacing from right
    subplot_top          subplot spacing from top
    subplot_wspace       horizontal spacing between subplots
    wl_fn                Winglink file name (full path)
    ==================== ======================================================

    =================== =======================================================
    Methods             Description
    =================== =======================================================
    plot                plots the apparent resistiviy and phase of data and
                        model if given.  called on instantiation if plot_yn
                        is 'y'.
    redraw_plot         call redraw_plot to redraw the figures, 
                        if one of the attributes has been changed
    save_figures        save all the matplotlib.figure instances in fig_list
    =================== =======================================================


    :Example: ::
        >>> data_fn = r"/home/occam/line1/inv1/OccamDataFile.dat"
        >>> resp_list = [r"/home/occam/line1/inv1/test_{0:02}".format(ii) 
                         for ii in range(2, 8, 2)]
        >>> pr_obj = occam2d.PlotResponse(data_fn, resp_list, plot_tipper='y')

    """

    def __init__(self, wl_data_fn=None, resp_fn=None, **kwargs):

        self.wl_data_fn = wl_data_fn

        self.color_mode = kwargs.pop('color_mode', 'color')

        self.ms = kwargs.pop('ms', 1.5)
        self.lw = kwargs.pop('lw', .5)

        self.ax_list = []
        self.line_list = []
        self.err_list = []

        # color mode
        if self.color_mode == 'color':
            # color for data
            self.cted = kwargs.pop('cted', (0, 0, 1))
            self.ctmd = kwargs.pop('ctmd', (1, 0, 0))
            self.mted = kwargs.pop('mted', 's')
            self.mtmd = kwargs.pop('mtmd', 'o')

            # color for Winglink model
            self.ctewl = kwargs.pop('ctewl', (0, .6, .8))
            self.ctmwl = kwargs.pop('ctmwl', (.8, .7, 0))
            self.mtewl = kwargs.pop('mtewl', 'x')
            self.mtmwl = kwargs.pop('mtmwl', 'x')

            # color of tipper
            self.ctipr = kwargs.pop('ctipr', self.cted)
            self.ctipi = kwargs.pop('ctipi', self.ctmd)

        # black and white mode
        elif self.color_mode == 'bw':
            # color for data
            self.cted = kwargs.pop('cted', (0, 0, 0))
            self.ctmd = kwargs.pop('ctmd', (0, 0, 0))
            self.mted = kwargs.pop('mted', '*')
            self.mtmd = kwargs.pop('mtmd', 'v')

            # color for Winglink model
            self.ctewl = kwargs.pop('ctewl', (0.3, 0.3, 0.3))
            self.ctmwl = kwargs.pop('ctmwl', (0.3, 0.3, 0.3))
            self.mtewl = kwargs.pop('mtewl', '|')
            self.mtmwl = kwargs.pop('mtmwl', '_')

            self.ctipr = kwargs.pop('ctipr', self.cted)
            self.ctipi = kwargs.pop('ctipi', self.ctmd)

        self.phase_limits = kwargs.pop('phase_limits', (-5, 95))
        self.res_limits = kwargs.pop('res_limits', None)
        self.tip_limits = kwargs.pop('tip_limits', (-.5, .5))

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)

        self.subplot_wspace = .1
        self.subplot_hspace = .15
        self.subplot_right = .98
        self.subplot_left = .085
        self.subplot_top = .93
        self.subplot_bottom = .1

        self.font_size = kwargs.pop('font_size', 6)

        self.plot_type = kwargs.pop('plot_type', '1')
        self.plot_num = kwargs.pop('plot_num', 2)
        self.plot_tipper = kwargs.pop('plot_tipper', 'n')
        self.plot_yn = kwargs.pop('plot_yn', 'y')

        if self.plot_num == 1:
            self.ylabel_coord = kwargs.pop('ylabel_coords', (-.055, .5))
        elif self.plot_num == 2:
            self.ylabel_coord = kwargs.pop('ylabel_coords', (-.12, .5))

        self.fig_list = []

        if self.plot_yn == 'y':
            self.plot()

    def plot(self):
        """
        plot the data and model response, if given, in individual plots.

        """

        wl_data = read_output_file(self.wl_data_fn)

        # create station list
        self.station_list = list(wl_data.keys())

        #---------------plot each respones in a different figure---------------
        if self.plot_type == '1':
            pstation_list = list(range(len(self.station_list)))

        else:
            if type(self.plot_type) is not list:
                self.plot_type = [self.plot_type]

            pstation_list = []
            for ii, station in enumerate(self.station_list):
                for pstation in self.plot_type:
                    if station.find(pstation) >= 0:
                        pstation_list.append(ii)

        # set the grid of subplots
        if self.plot_tipper == 'y':
            gs = gridspec.GridSpec(3, 2,
                                   wspace=self.subplot_wspace,
                                   left=self.subplot_left,
                                   top=self.subplot_top,
                                   bottom=self.subplot_bottom,
                                   right=self.subplot_right,
                                   hspace=self.subplot_hspace,
                                   height_ratios=[2, 1.5, 1])
        else:
            gs = gridspec.GridSpec(2, 2,
                                   wspace=self.subplot_wspace,
                                   left=self.subplot_left,
                                   top=self.subplot_top,
                                   bottom=self.subplot_bottom,
                                   right=self.subplot_right,
                                   hspace=self.subplot_hspace,
                                   height_ratios=[2, 1.5])

        #--> set default font size
        plt.rcParams['font.size'] = self.font_size

        # loop over each station to plot
        for ii, jj in enumerate(pstation_list):
            fig = plt.figure(self.station_list[jj],
                             self.fig_size, dpi=self.fig_dpi)
            plt.clf()

            #--> set subplot instances
            #---plot both TE and TM in same subplot---
            if self.plot_num == 1:
                axrte = fig.add_subplot(gs[0, :])
                axrtm = axrte
                axpte = fig.add_subplot(gs[1, :], sharex=axrte)
                axptm = axpte
                if self.plot_tipper == 'y':
                    axtipre = fig.add_subplot(gs[2, :], sharex=axrte)
                    axtipim = axtipre

            #---plot TE and TM in separate subplots---
            elif self.plot_num == 2:
                axrte = fig.add_subplot(gs[0, 0])
                axrtm = fig.add_subplot(gs[0, 1])
                axpte = fig.add_subplot(gs[1, 0], sharex=axrte)
                axptm = fig.add_subplot(gs[1, 1], sharex=axrtm)
                if self.plot_tipper == 'y':
                    axtipre = fig.add_subplot(gs[2, 0], sharex=axrte)
                    axtipim = fig.add_subplot(gs[2, 1], sharex=axrtm)

            # plot the data, it should be the same for all response files
            # empty lists for legend marker and label
            rlistte = []
            llistte = []
            rlisttm = []
            llisttm = []

            #--------------add in winglink responses------------------------
            wlrms = wl_data[self.station_list[jj]]['rms']
            wl_dict = wl_data[self.station_list[jj]]

            zrxy = np.nonzero(wl_dict['te_res'][0] != 0)[0]
            zryx = np.nonzero(wl_dict['tm_res'][0] != 0)[0]

            #--> plot data
            if len(zrxy) > 0:
                r1 = axrte.loglog(wl_dict['period'][zrxy],
                                  wl_dict['te_res'][0, zrxy],
                                  ls='-.',
                                  marker=self.mted,
                                  ms=self.ms,
                                  color=self.cted,
                                  mfc=self.cted,
                                  mec=self.cted,
                                  lw=self.lw)
                axpte.semilogx(wl_dict['period'][zrxy],
                               wl_dict['te_phase'][0, zrxy],
                               ls='-.',
                               marker=self.mted,
                               ms=self.ms,
                               color=self.cted,
                               mfc=self.cted,
                               mec=self.cted,
                               lw=self.lw)
                rlistte.append(r1[0])
                llistte.append('$Obs_{TE}$ ')
            if len(zryx) > 0:
                r2 = axrtm.loglog(wl_dict['period'][zryx],
                                  wl_dict['tm_res'][0, zryx],
                                  ls='-.',
                                  marker=self.mtmd,
                                  ms=self.ms,
                                  color=self.ctmd,
                                  mfc=self.ctmd,
                                  lw=self.lw)

                # plot winglink phase
                axptm.semilogx(wl_dict['period'][zryx],
                               wl_dict['tm_phase'][0, zryx],
                               ls='-.',
                               marker=self.mtmwl,
                               ms=self.ms,
                               color=self.ctmd,
                               mfc=self.ctmd,
                               lw=self.lw)

                rlisttm.append(r2[0])
                llisttm.append('$Obs_{TM}$ ')

            if self.plot_tipper == 'y':
                t_list = []
                t_label = []
                txy = np.nonzero(wl_dict['re_tip'][0])[0]
                tyx = np.nonzero(wl_dict['im_tip'][0])[0]
                #--> real tipper  data
                if len(txy) > 0:
                    per_list_p = []
                    tpr_list_p = []
                    per_list_n = []
                    tpr_list_n = []
                    for per, tpr in zip(wl_dict['period'][txy],
                                        wl_dict['re_tip'][0, txy]):
                        if tpr >= 0:
                            per_list_p.append(per)
                            tpr_list_p.append(tpr)
                        else:
                            per_list_n.append(per)
                            tpr_list_n.append(tpr)
                    if len(per_list_p) > 0:
                        m_line, s_line, b_line = axtipre.stem(per_list_p,
                                                              tpr_list_p,
                                                              markerfmt='^',
                                                              basefmt='k')
                        plt.setp(m_line, 'markerfacecolor', self.cted)
                        plt.setp(m_line, 'markeredgecolor', self.cted)
                        plt.setp(m_line, 'markersize', self.ms)
                        plt.setp(s_line, 'linewidth', self.lw)
                        plt.setp(s_line, 'color', self.cted)
                        plt.setp(b_line, 'linewidth', .01)
                        t_list.append(m_line)
                        t_label.append('Real')
                    if len(per_list_n) > 0:
                        m_line, s_line, b_line = axtipre.stem(per_list_n,
                                                              tpr_list_n,
                                                              markerfmt='v',
                                                              basefmt='k')
                        plt.setp(m_line, 'markerfacecolor', self.cted)
                        plt.setp(m_line, 'markeredgecolor', self.cted)
                        plt.setp(m_line, 'markersize', self.ms)
                        plt.setp(s_line, 'linewidth', self.lw)
                        plt.setp(s_line, 'color', self.cted)
                        plt.setp(b_line, 'linewidth', .01)
                        if len(t_list) == 0:
                            t_list.append(m_line)
                            t_label.append('Real')

                else:
                    pass
                if len(tyx) > 0:
                    per_list_p = []
                    tpi_list_p = []
                    per_list_n = []
                    tpi_list_n = []
                    for per, tpi in zip(wl_dict['period'][tyx],
                                        wl_dict['im_tip'][0, tyx]):
                        if tpi >= 0:
                            per_list_p.append(per)
                            tpi_list_p.append(tpi)
                        else:
                            per_list_n.append(per)
                            tpi_list_n.append(tpi)
                    if len(per_list_p) > 0:
                        m_line, s_line, b_line = axtipim.stem(per_list_p,
                                                              tpi_list_p,
                                                              markerfmt='^',
                                                              basefmt='k')
                        plt.setp(m_line, 'markerfacecolor', self.ctmd)
                        plt.setp(m_line, 'markeredgecolor', self.ctmd)
                        plt.setp(m_line, 'markersize', self.ms)
                        plt.setp(s_line, 'linewidth', self.lw)
                        plt.setp(s_line, 'color', self.ctmd)
                        plt.setp(b_line, 'linewidth', .01)
                        t_list.append(m_line)
                        t_label.append('Imag')
                    if len(per_list_n) > 0:
                        m_line, s_line, b_line = axtipim.stem(per_list_n,
                                                              tpi_list_n,
                                                              markerfmt='v',
                                                              basefmt='k')
                        plt.setp(m_line, 'markerfacecolor', self.ctmd)
                        plt.setp(m_line, 'markeredgecolor', self.ctmd)
                        plt.setp(m_line, 'markersize', self.ms)
                        plt.setp(s_line, 'linewidth', self.lw)
                        plt.setp(s_line, 'color', self.ctmd)
                        plt.setp(b_line, 'linewidth', .01)
                        if len(t_list) <= 1:
                            t_list.append(m_line)
                            t_label.append('Imag')

            #--> plot model response
            if len(zrxy) > 0:
                r5 = axrte.loglog(wl_dict['period'][zrxy],
                                  wl_dict['te_res'][1, zrxy],
                                  ls='-.',
                                  marker=self.mtewl,
                                  ms=self.ms,
                                  color=self.ctewl,
                                  mfc=self.ctewl,
                                  lw=self.lw)
                axpte.semilogx(wl_dict['period'][zrxy],
                               wl_dict['te_phase'][1, zrxy],
                               ls='-.',
                               marker=self.mtewl,
                               ms=self.ms,
                               color=self.ctewl,
                               mfc=self.ctewl,
                               lw=self.lw)
                rlistte.append(r5[0])
                llistte.append('$Mod_{TE}$ ' + '{0:.2f}'.format(wlrms))
            if len(zryx) > 0:
                r6 = axrtm.loglog(wl_dict['period'][zryx],
                                  wl_dict['tm_res'][1, zryx],
                                  ls='-.',
                                  marker=self.mtmwl,
                                  ms=self.ms,
                                  color=self.ctmwl,
                                  mfc=self.ctmwl,
                                  lw=self.lw)

                # plot winglink phase
                axptm.semilogx(wl_dict['period'][zryx],
                               wl_dict['tm_phase'][1, zryx],
                               ls='-.',
                               marker=self.mtmwl,
                               ms=self.ms,
                               color=self.ctmwl,
                               mfc=self.ctmwl,
                               lw=self.lw)

                rlisttm.append(r6[0])
                llisttm.append('$Mod_{TM}$ ' + '{0:.2f}'.format(wlrms))

            if self.plot_tipper == 'y':
                txy = np.nonzero(wl_dict['re_tip'][0])[0]
                tyx = np.nonzero(wl_dict['im_tip'][0])[0]
                #--> real tipper  data
                if len(txy) > 0:
                    per_list_p = []
                    tpr_list_p = []
                    per_list_n = []
                    tpr_list_n = []
                    for per, tpr in zip(wl_dict['period'][txy],
                                        wl_dict['re_tip'][1, txy]):
                        if tpr >= 0:
                            per_list_p.append(per)
                            tpr_list_p.append(tpr)
                        else:
                            per_list_n.append(per)
                            tpr_list_n.append(tpr)
                    if len(per_list_p) > 0:
                        m_line, s_line, b_line = axtipre.stem(per_list_p,
                                                              tpr_list_p,
                                                              markerfmt='^',
                                                              basefmt='k')
                        plt.setp(m_line, 'markerfacecolor', self.ctewl)
                        plt.setp(m_line, 'markeredgecolor', self.ctewl)
                        plt.setp(m_line, 'markersize', self.ms)
                        plt.setp(s_line, 'linewidth', self.lw)
                        plt.setp(s_line, 'color', self.ctewl)
                        plt.setp(b_line, 'linewidth', .01)
                    if len(per_list_n) > 0:
                        m_line, s_line, b_line = axtipre.stem(per_list_n,
                                                              tpr_list_n,
                                                              markerfmt='v',
                                                              basefmt='k')
                        plt.setp(m_line, 'markerfacecolor', self.ctewl)
                        plt.setp(m_line, 'markeredgecolor', self.ctewl)
                        plt.setp(m_line, 'markersize', self.ms)
                        plt.setp(s_line, 'linewidth', self.lw)
                        plt.setp(s_line, 'color', self.ctewl)
                        plt.setp(b_line, 'linewidth', .01)

                else:
                    pass
                if len(tyx) > 0:
                    per_list_p = []
                    tpi_list_p = []
                    per_list_n = []
                    tpi_list_n = []
                    for per, tpi in zip(wl_dict['period'][tyx],
                                        wl_dict['im_tip'][1, tyx]):
                        if tpi >= 0:
                            per_list_p.append(per)
                            tpi_list_p.append(tpi)
                        else:
                            per_list_n.append(per)
                            tpi_list_n.append(tpi)
                    if len(per_list_p) > 0:
                        m_line, s_line, b_line = axtipim.stem(per_list_p,
                                                              tpi_list_p,
                                                              markerfmt='^',
                                                              basefmt='k')
                        plt.setp(m_line, 'markerfacecolor', self.ctmwl)
                        plt.setp(m_line, 'markeredgecolor', self.ctmwl)
                        plt.setp(m_line, 'markersize', self.ms)
                        plt.setp(s_line, 'linewidth', self.lw)
                        plt.setp(s_line, 'color', self.ctmwl)
                        plt.setp(b_line, 'linewidth', .01)
                    if len(per_list_n) > 0:
                        m_line, s_line, b_line = axtipim.stem(per_list_n,
                                                              tpi_list_n,
                                                              markerfmt='v',
                                                              basefmt='k')
                        plt.setp(m_line, 'markerfacecolor', self.ctmwl)
                        plt.setp(m_line, 'markeredgecolor', self.ctmwl)
                        plt.setp(m_line, 'markersize', self.ms)
                        plt.setp(s_line, 'linewidth', self.lw)
                        plt.setp(s_line, 'color', self.ctmwl)
                        plt.setp(b_line, 'linewidth', .01)

                else:
                    pass
            else:
                if self.plot_num == 1:
                    axrte.set_title(self.station_list[jj],
                                    fontdict={'size': self.font_size + 2,
                                              'weight': 'bold'})
                elif self.plot_num == 2:
                    fig.suptitle(self.station_list[jj],
                                 fontdict={'size': self.font_size + 2,
                                           'weight': 'bold'})

            # set the axis properties
            ax_list = [axrte, axrtm]
            for aa, axr in enumerate(ax_list):
                # set both axes to logarithmic scale
                axr.set_xscale('log', nonposx='clip')

                try:
                    axr.set_yscale('log', nonposy='clip')
                except ValueError:
                    pass

                # put on a grid
                axr.grid(True, alpha=.3, which='both', lw=.5 * self.lw)
                axr.yaxis.set_label_coords(self.ylabel_coord[0],
                                           self.ylabel_coord[1])

                # set resistivity limits if desired
                if self.res_limits != None:
                    axr.set_ylim(10**self.res_limits[0],
                                 10**self.res_limits[1])

                # set the tick labels to invisible
                plt.setp(axr.xaxis.get_ticklabels(), visible=False)
                if aa == 0:
                    axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                                   fontdict={'size': self.font_size + 2,
                                             'weight': 'bold'})

                # set legend based on the plot type
                if self.plot_num == 1:
                    if aa == 0:
                        axr.legend(rlistte + rlisttm, llistte + llisttm,
                                   loc=2, markerscale=1,
                                   borderaxespad=.05,
                                   labelspacing=.08,
                                   handletextpad=.15,
                                   borderpad=.05,
                                   prop={'size': self.font_size + 1})
                elif self.plot_num == 2:
                    if aa == 0:
                        axr.legend(rlistte,
                                   llistte,
                                   loc=2, markerscale=1,
                                   borderaxespad=.05,
                                   labelspacing=.08,
                                   handletextpad=.15,
                                   borderpad=.05,
                                   prop={'size': self.font_size + 1})

                    if aa == 1:
                        axr.legend(rlisttm,
                                   llisttm,
                                   loc=2, markerscale=1,
                                   borderaxespad=.05,
                                   labelspacing=.08,
                                   handletextpad=.15,
                                   borderpad=.05,
                                   prop={'size': self.font_size + 1})

            # set Properties for the phase axes
            for aa, axp in enumerate([axpte, axptm]):
                # set the x-axis to log scale
                axp.set_xscale('log', nonposx='clip')

                # set the phase limits
                axp.set_ylim(self.phase_limits)

                # put a grid on the subplot
                axp.grid(True, alpha=.3, which='both', lw=.5 * self.lw)

                # set the tick locations
                axp.yaxis.set_major_locator(MultipleLocator(10))
                axp.yaxis.set_minor_locator(MultipleLocator(2))

                # set the x axis label
                if self.plot_tipper == 'y':
                    plt.setp(axp.get_xticklabels(), visible=False)
                else:
                    axp.set_xlabel('Period (s)',
                                   fontdict={'size': self.font_size + 2,
                                             'weight': 'bold'})

                # put the y label on the far left plot
                axp.yaxis.set_label_coords(self.ylabel_coord[0],
                                           self.ylabel_coord[1])
                if aa == 0:
                    axp.set_ylabel('Phase (deg)',
                                   fontdict={'size': self.font_size + 2,
                                             'weight': 'bold'})

            # set axes properties of tipper axis
            if self.plot_tipper == 'y':
                for aa, axt in enumerate([axtipre, axtipim]):
                    axt.set_xscale('log', nonposx='clip')

                    # set tipper limits
                    axt.set_ylim(self.tip_limits)

                    # put a grid on the subplot
                    axt.grid(True, alpha=.3, which='both', lw=.5 * self.lw)

                    # set the tick locations
                    axt.yaxis.set_major_locator(MultipleLocator(.2))
                    axt.yaxis.set_minor_locator(MultipleLocator(.1))

                    # set the x axis label
                    axt.set_xlabel('Period (s)',
                                   fontdict={'size': self.font_size + 2,
                                             'weight': 'bold'})

                    axt.set_xlim(10**np.floor(np.log10(wl_dict['period'].min())),
                                 10**np.ceil(np.log10(wl_dict['period'].max())))

                    # put the y label on the far left plot
                    axt.yaxis.set_label_coords(self.ylabel_coord[0],
                                               self.ylabel_coord[1])
                    if aa == 0:
                        axt.set_ylabel('Tipper',
                                       fontdict={'size': self.font_size + 2,
                                                 'weight': 'bold'})
                        if self.plot_num == 2:
                            axt.text(axt.get_xlim()[0] * 1.25,
                                     self.tip_limits[1] * .9,
                                     'Real', horizontalalignment='left',
                                     verticalalignment='top',
                                     bbox={'facecolor': 'white'},
                                     fontdict={'size': self.font_size + 1})
                        else:
                            axt.legend(t_list, t_label,
                                       loc=2, markerscale=1,
                                       borderaxespad=.05,
                                       labelspacing=.08,
                                       handletextpad=.15,
                                       borderpad=.05,
                                       prop={'size': self.font_size + 1})
                    if aa == 1:
                        if self.plot_num == 2:
                            axt.text(axt.get_xlim()[0] * 1.25,
                                     self.tip_limits[1] * .9,
                                     'Imag', horizontalalignment='left',
                                     verticalalignment='top',
                                     bbox={'facecolor': 'white'},
                                     fontdict={'size': self.font_size + 1})

            # make sure the axis and figure are accessible to the user
            self.fig_list.append({'station': self.station_list[jj],
                                  'fig': fig, 'axrte': axrte, 'axrtm': axrtm,
                                  'axpte': axpte, 'axptm': axptm})

        # set the plot to be full screen well at least try
        plt.show()

    def redraw_plot(self):
        """
        redraw plot if parameters were changed

        use this function if you updated some attributes and want to re-plot.

        :Example: ::

            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plot2DResponses()
            >>> #change color of te markers to a gray-blue
            >>> p1.cted = (.5, .5, .7)
            >>> p1.redraw_plot()
        """

        plt.close('all')
        self.plot()

    def save_figures(self, save_path, fig_fmt='pdf', fig_dpi=None,
                     close_fig='y'):
        """
        save all the figure that are in self.fig_list

        :Example: ::

            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plot2DResponses()
            >>> p1.save_figures(r"/home/occam2d/Figures", fig_fmt='jpg')
        """

        if not os.path.exists(save_path):
            os.mkdir(save_path)

        for fdict in self.fig_list:
            svfn = '{0}_resp.{1}'.format(fdict['station'], fig_fmt)
            fdict['fig'].savefig(os.path.join(save_path, svfn),
                                 dpi=self.fig_dpi)
            if close_fig == 'y':
                plt.close(fdict['fig'])

            print("saved figure to {0}".format(os.path.join(save_path, svfn)))


#==============================================================================
# plot pseudo section of data and model response
#==============================================================================
class PlotPseudoSection(object):
    """
    plot a pseudo section of the data and response if given


    Arguments:
    -------------
        **wl_data_fn** : string
                         full path to winglink output data file.

    ==================== ======================================================
    key words            description
    ==================== ======================================================
    axmpte               matplotlib.axes instance for TE model phase
    axmptm               matplotlib.axes instance for TM model phase
    axmrte               matplotlib.axes instance for TE model app. res 
    axmrtm               matplotlib.axes instance for TM model app. res 
    axpte                matplotlib.axes instance for TE data phase 
    axptm                matplotlib.axes instance for TM data phase
    axrte                matplotlib.axes instance for TE data app. res.
    axrtm                matplotlib.axes instance for TM data app. res.
    cb_pad               padding between colorbar and axes
    cb_shrink            percentage to shrink the colorbar to
    fig                  matplotlib.figure instance
    fig_dpi              resolution of figure in dots per inch
    fig_num              number of figure instance
    fig_size             size of figure in inches (width, height)
    font_size            size of font in points
    label_list            list to label plots
    ml                   factor to label stations if 2 every other station
                         is labeled on the x-axis
    period               np.array of periods to plot
    phase_cmap           color map name of phase
    phase_limits_te      limits for te phase in degrees (min, max)
    phase_limits_tm      limits for tm phase in degrees (min, max)            
    plot_resp            [ 'y' | 'n' ] to plot response
    plot_tipper          [ 'y' | 'n' ] to plot tipper
    plot_yn              [ 'y' | 'n' ] 'y' to plot on instantiation
    res_cmap             color map name for resistivity
    res_limits_te        limits for te resistivity in log scale (min, max)
    res_limits_tm        limits for tm resistivity in log scale (min, max)
    rp_list               list of dictionaries as made from read2Dresp
    station_id           index to get station name (min, max)
    station_list          station list got from rp_list
    subplot_bottom       subplot spacing from bottom (relative coordinates) 
    subplot_hspace       vertical spacing between subplots
    subplot_left         subplot spacing from left  
    subplot_right        subplot spacing from right
    subplot_top          subplot spacing from top
    subplot_wspace       horizontal spacing between subplots
    ==================== ======================================================

    =================== =======================================================
    Methods             Description
    =================== =======================================================
    plot                plots a pseudo-section of apparent resistiviy and phase
                        of data and model if given.  called on instantiation 
                        if plot_yn is 'y'.
    redraw_plot         call redraw_plot to redraw the figures, 
                        if one of the attributes has been changed
    save_figure         saves the matplotlib.figure instance to desired 
                        location and format
    =================== =======================================================

   :Example: ::

        >>> import mtpy.modeling.winglink as winglink
        >>> d_fn = r"/home/winglink/Line1/Inv1/DataRW.txt"
        >>> ps_plot = winglink.PlotPseudoSection(d_fn) 

    """

    def __init__(self, wl_data_fn=None, **kwargs):

        self.wl_data_fn = wl_data_fn

        self.plot_resp = kwargs.pop('plot_resp', 'y')

        self.label_list = [r'$\rho_{TE-Data}$', r'$\rho_{TE-Model}$',
                           r'$\rho_{TM-Data}$', r'$\rho_{TM-Model}$',
                           '$\phi_{TE-Data}$', '$\phi_{TE-Model}$',
                           '$\phi_{TM-Data}$', '$\phi_{TM-Model}$',
                           '$\Re e\{T_{Data}\}$', '$\Re e\{T_{Model}\}$',
                           '$\Im m\{T_{Data}\}$', '$\Im m\{T_{Model}\}$']

        self.phase_limits_te = kwargs.pop('phase_limits_te', (-5, 95))
        self.phase_limits_tm = kwargs.pop('phase_limits_tm', (-5, 95))
        self.res_limits_te = kwargs.pop('res_limits_te', (0, 3))
        self.res_limits_tm = kwargs.pop('res_limits_tm', (0, 3))
        self.tip_limits_re = kwargs.pop('tip_limits_re', (-1, 1))
        self.tip_limits_im = kwargs.pop('tip_limits_im', (-1, 1))

        self.phase_cmap = kwargs.pop('phase_cmap', 'jet')
        self.res_cmap = kwargs.pop('res_cmap', 'jet_r')
        self.tip_cmap = kwargs.pop('res_cmap', 'Spectral_r')

        self.ml = kwargs.pop('ml', 2)
        self.station_id = kwargs.pop('station_id', [0, 4])

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)

        self.subplot_wspace = .025
        self.subplot_hspace = .0
        self.subplot_right = .95
        self.subplot_left = .085
        self.subplot_top = .97
        self.subplot_bottom = .1

        self.font_size = kwargs.pop('font_size', 6)

        self.plot_type = kwargs.pop('plot_type', '1')
        self.plot_num = kwargs.pop('plot_num', 2)
        self.plot_tipper = kwargs.pop('plot_tipper', 'n')
        self.plot_yn = kwargs.pop('plot_yn', 'y')

        self.cb_shrink = .7
        self.cb_pad = .015

        self.axrte = None
        self.axrtm = None
        self.axpte = None
        self.axptm = None
        self.axmrte = None
        self.axmrtm = None
        self.axmpte = None
        self.axmptm = None
        self.axtpr = None
        self.axtpi = None
        self.axmtpr = None
        self.axmtpi = None

        self.te_res_arr = None
        self.tm_res_arr = None
        self.te_phase_arr = None
        self.tm_phase_arr = None
        self.tip_real_arr = None
        self.tip_imag_arr = None

        self.fig = None

        if self.plot_yn == 'y':
            self.plot()

    def plot(self):
        """
        plot pseudo section of data and response if given

        """
        if self.plot_resp == 'y':
            nr = 2
        else:
            nr = 1

        wl_data = read_output_file(self.wl_data_fn)
        stations = list(wl_data.keys())

        #--> need to sort the stations to be in order
        slst = np.array([(ss, wl_data[ss]['index']) for ss in stations],
                        dtype=[('station', '|S20'), ('index', np.int)])
        slst.sort(order='index')
        stations = slst['station']

        ns = len(list(wl_data.keys()))
        #--> get periods in decreasing order for easier plotting
        periods = []
        for ss in stations:
            periods.extend(list(wl_data[ss]['period']))

        periods = np.array(sorted(set(periods), reverse=True))
        # need to make a dictionary of where periods go cause they are not all
        # the same from winglink
        p_dict = dict([(per, index) for index, per in enumerate(periods)])
        nf = periods.shape[0]
        # set limits of y-axis in plot
        ylimits = (periods.max(), periods.min())

        # make a grid for pcolormesh so you can have a log scale
        # get things into arrays for plotting
        te_res_arr = np.ones((nf, ns, nr))
        tm_res_arr = np.ones((nf, ns, nr))
        te_phase_arr = np.zeros((nf, ns, nr))
        tm_phase_arr = np.zeros((nf, ns, nr))
        tip_real_arr = np.zeros((nf, ns, nr))
        tip_imag_arr = np.zeros((nf, ns, nr))

        for ss in stations:
            d_dict = wl_data[ss]
            ii = d_dict['index']
            for jj, per in enumerate(d_dict['period']):
                p_index = p_dict[per]
                te_res_arr[p_index, ii, 0] = d_dict['te_res'][0, jj]
                tm_res_arr[p_index, ii, 0] = d_dict['tm_res'][0, jj]
                te_phase_arr[p_index, ii, 0] = d_dict['te_phase'][0, jj]
                tm_phase_arr[p_index, ii, 0] = d_dict['tm_phase'][0, jj]
                tip_real_arr[p_index, ii, 0] = d_dict['re_tip'][0, jj]
                tip_imag_arr[p_index, ii, 0] = d_dict['im_tip'][0, jj]

                # read in response data
                if self.plot_resp == 'y':
                    te_res_arr[p_index, ii, 1] = d_dict['te_res'][1, jj]
                    tm_res_arr[p_index, ii, 1] = d_dict['tm_res'][1, jj]
                    te_phase_arr[p_index, ii, 1] = d_dict['te_phase'][1, jj]
                    tm_phase_arr[p_index, ii, 1] = d_dict['tm_phase'][1, jj]
                    tip_real_arr[p_index, ii, 1] = d_dict['re_tip'][1, jj]
                    tip_imag_arr[p_index, ii, 1] = d_dict['im_tip'][1, jj]

        # need to make any zeros 1 for taking log10
        te_res_arr[np.where(te_res_arr == 0)] = 1.0
        tm_res_arr[np.where(tm_res_arr == 0)] = 1.0

        self.te_res_arr = te_res_arr
        self.tm_res_arr = tm_res_arr
        self.te_phase_arr = te_phase_arr
        self.tm_phase_arr = tm_phase_arr
        self.tip_real_arr = tip_real_arr
        self.tip_imag_arr = tip_imag_arr

        # need to extend the last grid cell because meshgrid expects n+1 cells
        offset_list = np.arange(ns + 1)
        # make a meshgrid for plotting
        # flip frequency so bottom corner is long period
        dgrid, fgrid = np.meshgrid(offset_list, periods)

        # make list for station labels
        sindex_1 = self.station_id[0]
        sindex_2 = self.station_id[1]
        slabel = [stations[ss][sindex_1:sindex_2]
                  for ss in range(0, ns, self.ml)]

        xloc = offset_list[0] + abs(offset_list[0] - offset_list[1]) / 5
        yloc = 1.10 * periods[-2]

        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.left'] = self.subplot_left

        log_labels_te = ['10$^{0}$'.format('{' + str(nn) + '}')
                         for nn in np.arange(int(self.res_limits_te[0]),
                                             int(self.res_limits_te[1]) + 1)]
        log_labels_tm = ['10$^{0}$'.format('{' + str(nn) + '}')
                         for nn in np.arange(int(self.res_limits_tm[0]),
                                             int(self.res_limits_tm[1]) + 1)]

        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()

        if self.plot_resp == 'y':
            if self.plot_tipper == 'y':
                gs1 = gridspec.GridSpec(1, 3,
                                        left=self.subplot_left,
                                        right=self.subplot_right,
                                        wspace=self.subplot_wspace)
                gs4 = gridspec.GridSpecFromSubplotSpec(2, 2,
                                                       hspace=self.subplot_hspace,
                                                       wspace=0,
                                                       subplot_spec=gs1[2])
            else:
                gs1 = gridspec.GridSpec(1, 2,
                                        left=self.subplot_left,
                                        right=self.subplot_right,
                                        wspace=self.subplot_wspace)
            gs2 = gridspec.GridSpecFromSubplotSpec(2, 2,
                                                   hspace=self.subplot_hspace,
                                                   wspace=0,
                                                   subplot_spec=gs1[0])
            gs3 = gridspec.GridSpecFromSubplotSpec(2, 2,
                                                   hspace=self.subplot_hspace,
                                                   wspace=0,
                                                   subplot_spec=gs1[1])

            # plot TE resistivity data
            self.axrte = plt.Subplot(self.fig, gs2[0, 0])
            self.fig.add_subplot(self.axrte)
            self.axrte.pcolormesh(dgrid,
                                  fgrid,
                                  np.log10(te_res_arr[:, :, 0]),
                                  cmap=self.res_cmap,
                                  vmin=self.res_limits_te[0],
                                  vmax=self.res_limits_te[1])

            # plot TE resistivity model
            self.axmrte = plt.Subplot(self.fig, gs2[0, 1])
            self.fig.add_subplot(self.axmrte)
            self.axmrte.pcolormesh(dgrid,
                                   fgrid,
                                   np.log10(te_res_arr[:, :, 1]),
                                   cmap=self.res_cmap,
                                   vmin=self.res_limits_te[0],
                                   vmax=self.res_limits_te[1])

            # plot TM resistivity data
            self.axrtm = plt.Subplot(self.fig, gs3[0, 0])
            self.fig.add_subplot(self.axrtm)
            self.axrtm.pcolormesh(dgrid,
                                  fgrid,
                                  np.log10(tm_res_arr[:, :, 0]),
                                  cmap=self.res_cmap,
                                  vmin=self.res_limits_tm[0],
                                  vmax=self.res_limits_tm[1])

            # plot TM resistivity model
            self.axmrtm = plt.Subplot(self.fig, gs3[0, 1])
            self.fig.add_subplot(self.axmrtm)
            self.axmrtm.pcolormesh(dgrid,
                                   fgrid,
                                   np.log10(tm_res_arr[:, :, 1]),
                                   cmap=self.res_cmap,
                                   vmin=self.res_limits_tm[0],
                                   vmax=self.res_limits_tm[1])

            # plot TE phase data
            self.axpte = plt.Subplot(self.fig, gs2[1, 0])
            self.fig.add_subplot(self.axpte)
            self.axpte.pcolormesh(dgrid,
                                  fgrid,
                                  te_phase_arr[:, :, 0],
                                  cmap=self.phase_cmap,
                                  vmin=self.phase_limits_te[0],
                                  vmax=self.phase_limits_te[1])

            # plot TE phase model
            self.axmpte = plt.Subplot(self.fig, gs2[1, 1])
            self.fig.add_subplot(self.axmpte)
            self.axmpte.pcolormesh(dgrid,
                                   fgrid,
                                   te_phase_arr[:, :, 1],
                                   cmap=self.phase_cmap,
                                   vmin=self.phase_limits_te[0],
                                   vmax=self.phase_limits_te[1])

            # plot TM phase data
            self.axptm = plt.Subplot(self.fig, gs3[1, 0])
            self.fig.add_subplot(self.axptm)
            self.axptm.pcolormesh(dgrid,
                                  fgrid,
                                  tm_phase_arr[:, :, 0],
                                  cmap=self.phase_cmap,
                                  vmin=self.phase_limits_tm[0],
                                  vmax=self.phase_limits_tm[1])

            # plot TM phase model
            self.axmptm = plt.Subplot(self.fig, gs3[1, 1])
            self.fig.add_subplot(self.axmptm)
            self.axmptm.pcolormesh(dgrid,
                                   fgrid,
                                   tm_phase_arr[:, :, 1],
                                   cmap=self.phase_cmap,
                                   vmin=self.phase_limits_tm[0],
                                   vmax=self.phase_limits_tm[1])

            ax_list = [self.axrte, self.axmrte, self.axrtm, self.axmrtm,
                       self.axpte, self.axmpte, self.axptm, self.axmptm]

            if self.plot_tipper == 'y':
                # plot real tipper  data
                self.axtpr = plt.Subplot(self.fig, gs4[0, 0])
                self.fig.add_subplot(self.axtpr)
                self.axtpr.pcolormesh(dgrid,
                                      fgrid,
                                      tip_real_arr[:, :, 0],
                                      cmap=self.tip_cmap,
                                      vmin=self.tip_limits_re[0],
                                      vmax=self.tip_limits_re[1])
                # plot real tipper  model
                self.axmtpr = plt.Subplot(self.fig, gs4[0, 1])
                self.fig.add_subplot(self.axmtpr)
                self.axmtpr.pcolormesh(dgrid,
                                       fgrid,
                                       tip_real_arr[:, :, 1],
                                       cmap=self.tip_cmap,
                                       vmin=self.tip_limits_re[0],
                                       vmax=self.tip_limits_re[1])

                # plot imag tipper  data
                self.axtpi = plt.Subplot(self.fig, gs4[1, 0])
                self.fig.add_subplot(self.axtpi)
                self.axtpi.pcolormesh(dgrid,
                                      fgrid,
                                      tip_imag_arr[:, :, 0],
                                      cmap=self.tip_cmap,
                                      vmin=self.tip_limits_re[0],
                                      vmax=self.tip_limits_re[1])
                # plot imag tipper  model
                self.axmtpi = plt.Subplot(self.fig, gs4[1, 1])
                self.fig.add_subplot(self.axmtpi)
                self.axmtpi.pcolormesh(dgrid,
                                       fgrid,
                                       tip_imag_arr[:, :, 1],
                                       cmap=self.tip_cmap,
                                       vmin=self.tip_limits_re[0],
                                       vmax=self.tip_limits_re[1])

                ax_list.append(self.axtpr)
                ax_list.append(self.axmtpr)
                ax_list.append(self.axtpi)
                ax_list.append(self.axmtpi)

            # make everthing look tidy
            for xx, ax in enumerate(ax_list):
                ax.semilogy()
                ax.set_ylim(ylimits)
                ax.xaxis.set_ticks(offset_list[np.arange(0, ns, self.ml)])
                ax.xaxis.set_ticks(offset_list, minor=True)
                ax.xaxis.set_ticklabels(slabel)
                ax.set_xlim(offset_list.min(), offset_list.max())
                if np.remainder(xx, 2.0) == 1:
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    cbx = mcb.make_axes(ax,
                                        shrink=self.cb_shrink,
                                        pad=self.cb_pad)
                if xx == 2 or xx == 6 or xx == 8 or xx == 10:
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)

                if xx < 4:
                    plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                    if xx == 1:
                        cb = mcb.ColorbarBase(cbx[0], cmap=self.res_cmap,
                                              norm=Normalize(vmin=self.res_limits_te[0],
                                                             vmax=self.res_limits_te[1]))
                        cb.set_ticks(np.arange(int(self.res_limits_te[0]),
                                               int(self.res_limits_te[1]) + 1))
                        cb.set_ticklabels(log_labels_te)
                    if xx == 3:
                        cb = mcb.ColorbarBase(cbx[0], cmap=self.res_cmap,
                                              norm=Normalize(vmin=self.res_limits_tm[0],
                                                             vmax=self.res_limits_tm[1]))
                        cb.set_label('App. Res. ($\Omega \cdot$m)',
                                     fontdict={'size': self.font_size + 1,
                                               'weight': 'bold'})
                        cb.set_label('Resistivity ($\Omega \cdot$m)',
                                     fontdict={'size': self.font_size + 1,
                                               'weight': 'bold'})
                        cb.set_ticks(np.arange(int(self.res_limits_tm[0]),
                                               int(self.res_limits_tm[1]) + 1))
                        cb.set_ticklabels(log_labels_tm)
                else:
                    # color bar TE phase
                    if xx == 5:
                        cb = mcb.ColorbarBase(cbx[0], cmap=self.phase_cmap,
                                              norm=Normalize(vmin=self.phase_limits_te[0],
                                                             vmax=self.phase_limits_te[1]))
                    # color bar TM phase
                    if xx == 7:
                        cb = mcb.ColorbarBase(cbx[0], cmap=self.phase_cmap,
                                              norm=Normalize(vmin=self.phase_limits_tm[0],
                                                             vmax=self.phase_limits_tm[1]))
                        cb.set_label('Phase (deg)',
                                     fontdict={'size': self.font_size + 1,
                                               'weight': 'bold'})
                    # color bar tipper Imag
                    if xx == 9:
                        cb = mcb.ColorbarBase(cbx[0], cmap=self.tip_cmap,
                                              norm=Normalize(vmin=self.tip_limits_re[0],
                                                             vmax=self.tip_limits_re[1]))
                        cb.set_label('Re{T}',
                                     fontdict={'size': self.font_size + 1,
                                               'weight': 'bold'})
                    if xx == 11:
                        cb = mcb.ColorbarBase(cbx[0], cmap=self.tip_cmap,
                                              norm=Normalize(vmin=self.tip_limits_im[0],
                                                             vmax=self.tip_limits_im[1]))
                        cb.set_label('Im{T}',
                                     fontdict={'size': self.font_size + 1,
                                               'weight': 'bold'})

                ax.text(xloc, yloc, self.label_list[xx],
                        fontdict={'size': self.font_size + 1},
                        bbox={'facecolor': 'white'},
                        horizontalalignment='left',
                        verticalalignment='top')
                if xx == 0 or xx == 4:
                    ax.set_ylabel('Period (s)',
                                  fontdict={'size': self.font_size + 2,
                                            'weight': 'bold'})
                if xx > 3:
                    ax.set_xlabel('Station', fontdict={'size': self.font_size + 2,
                                                       'weight': 'bold'})

            plt.show()

        else:
            if self.plot_tipper == 'y':
                gs1 = gridspec.GridSpec(2, 3,
                                        left=self.subplot_left,
                                        right=self.subplot_right,
                                        hspace=self.subplot_hspace,
                                        wspace=self.subplot_wspace)

            else:
                gs1 = gridspec.GridSpec(2, 2,
                                        left=self.subplot_left,
                                        right=self.subplot_right,
                                        hspace=self.subplot_hspace,
                                        wspace=self.subplot_wspace)

            # plot TE resistivity data
            self.axrte = self.fig.add_subplot(gs1[0, 0])
            self.axrte.pcolormesh(dgrid,
                                  fgrid,
                                  np.log10(te_res_arr[:, :, 0]),
                                  cmap=self.res_cmap,
                                  vmin=self.res_limits_te[0],
                                  vmax=self.res_limits_te[1])

            # plot TM resistivity data
            self.axrtm = self.fig.add_subplot(gs1[0, 1])
            self.axrtm.pcolormesh(dgrid,
                                  fgrid,
                                  np.log10(tm_res_arr[:, :, 0]),
                                  cmap=self.res_cmap,
                                  vmin=self.res_limits_tm[0],
                                  vmax=self.res_limits_tm[1])

            # plot TE phase data
            self.axpte = self.fig.add_subplot(gs1[1, 0])
            self.axpte.pcolormesh(dgrid,
                                  fgrid,
                                  te_phase_arr[:, :, 0],
                                  cmap=self.phase_cmap,
                                  vmin=self.phase_limits_te[0],
                                  vmax=self.phase_limits_te[1])

            # plot TM phase data
            self.axptm = self.fig.add_subplot(gs1[1, 1])
            self.axptm.pcolormesh(dgrid,
                                  fgrid,
                                  tm_phase_arr[:, :, 0],
                                  cmap=self.phase_cmap,
                                  vmin=self.phase_limits_tm[0],
                                  vmax=self.phase_limits_tm[1])
            ax_list = [self.axrte, self.axrtm, self.axpte, self.axptm]
            if self.plot_tipper == 'y':
                # plot real tipper  data
                self.axtpr = plt.Subplot(self.fig, gs1[0, 2])
                self.fig.add_subplot(self.axtpr)
                self.axtpr.pcolormesh(dgrid,
                                      fgrid,
                                      tip_real_arr[:, :, 0],
                                      cmap=self.tip_cmap,
                                      vmin=self.tip_limits_re[0],
                                      vmax=self.tip_limits_re[1])
                # plot real tipper  data
                self.axtpi = plt.Subplot(self.fig, gs1[1, 2])
                self.fig.add_subplot(self.axtpi)
                self.axtpi.pcolormesh(dgrid,
                                      fgrid,
                                      tip_imag_arr[:, :, 0],
                                      cmap=self.tip_cmap,
                                      vmin=self.tip_limits_re[0],
                                      vmax=self.tip_limits_re[1])
                ax_list.append(self.axtpr)
                ax_list.append(self.axtpi)

            # make everything look tidy
            for xx, ax in enumerate(ax_list):
                ax.semilogy()
                ax.set_ylim(ylimits)
                ax.xaxis.set_ticks(offset_list[np.arange(0, ns, self.ml)])
                ax.xaxis.set_ticks(offset_list, minor=True)
                ax.xaxis.set_ticklabels(slabel)
                ax.grid(True, alpha=.25)
                ax.set_xlim(offset_list.min(), offset_list.max())
                cbx = mcb.make_axes(ax,
                                    shrink=self.cb_shrink,
                                    pad=self.cb_pad)
                if xx == 0:
                    plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                    cb = mcb.ColorbarBase(cbx[0], cmap=self.res_cmap,
                                          norm=Normalize(vmin=self.res_limits_te[0],
                                                         vmax=self.res_limits_te[1]))
                    cb.set_ticks(np.arange(self.res_limits_te[0],
                                           self.res_limits_te[1] + 1))
                    cb.set_ticklabels(log_labels_te)
                elif xx == 1:
                    plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)

                    cb = mcb.ColorbarBase(cbx[0], cmap=self.res_cmap,
                                          norm=Normalize(vmin=self.res_limits_tm[0],
                                                         vmax=self.res_limits_tm[1]))
                    cb.set_label('App. Res. ($\Omega \cdot$m)',
                                 fontdict={'size': self.font_size + 1,
                                           'weight': 'bold'})
                    cb.set_ticks(np.arange(self.res_limits_tm[0],
                                           self.res_limits_tm[1] + 1))
                    cb.set_ticklabels(log_labels_tm)
                elif xx == 2:
                    cb = mcb.ColorbarBase(cbx[0], cmap=self.phase_cmap,
                                          norm=Normalize(vmin=self.phase_limits_te[0],
                                                         vmax=self.phase_limits_te[1]))
                    cb.set_ticks(np.arange(self.phase_limits_te[0],
                                           self.phase_limits_te[1] + 1, 15))
                elif xx == 3:
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    cb = mcb.ColorbarBase(cbx[0], cmap=self.phase_cmap,
                                          norm=Normalize(vmin=self.phase_limits_tm[0],
                                                         vmax=self.phase_limits_tm[1]))
                    cb.set_label('Phase (deg)',
                                 fontdict={'size': self.font_size + 1,
                                           'weight': 'bold'})
                    cb.set_ticks(np.arange(self.phase_limits_te[0],
                                           self.phase_limits_te[1] + 1, 15))

                # real tipper
                elif xx == 4:
                    plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    cb = mcb.ColorbarBase(cbx[0], cmap=self.tip_cmap,
                                          norm=Normalize(vmin=self.tip_limits_re[0],
                                                         vmax=self.tip_limits_re[1]))
                    cb.set_label('Re{T}',
                                 fontdict={'size': self.font_size + 1,
                                           'weight': 'bold'})
                # imag tipper
                elif xx == 5:
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    cb = mcb.ColorbarBase(cbx[0], cmap=self.tip_cmap,
                                          norm=Normalize(vmin=self.tip_limits_im[0],
                                                         vmax=self.tip_limits_im[1]))
                    cb.set_label('Im{T}',
                                 fontdict={'size': self.font_size + 1,
                                           'weight': 'bold'})

                ax.text(xloc, yloc, self.label_list[2 * xx],
                        fontdict={'size': self.font_size + 1},
                        bbox={'facecolor': 'white'},
                        horizontalalignment='left',
                        verticalalignment='top')
                if xx == 0 or xx == 2:
                    ax.set_ylabel('Period (s)',
                                  fontdict={'size': self.font_size + 2,
                                            'weight': 'bold'})
                if xx > 1:
                    ax.set_xlabel('Station', fontdict={'size': self.font_size + 2,
                                                       'weight': 'bold'})

            plt.show()

    def redraw_plot(self):
        """
        redraw plot if parameters were changed

        use this function if you updated some attributes and want to re-plot.

        :Example: ::

            >>> # plot tipper and change station id
            >>> import mtpy.modeling.winglink as winglink
            >>> ps_plot = winglink.PlotPseudosection(wl_fn)
            >>> ps_plot.plot_tipper = 'y'
            >>> ps_plot.station_id = [2, 5]
            >>> #label only every 3rd station
            >>> ps_plot.ml = 3
            >>> ps_plot.redraw_plot()
        """

        plt.close(self.fig)
        self.plot()

    def save_figure(self, save_fn, file_format='pdf', orientation='portrait',
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
            >>> ps_plot.save_plot(r'/home/MT/figures', file_format='jpg')

        """

        if fig_dpi == None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        else:
            save_fn = os.path.join(save_fn, 'WLPseudoSection.' +
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

    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.  

        Use this if you change an of the .fig or axes properties

        :Example: ::

            >>> # to change the grid lines to only be on the major ticks
            >>> [ax.grid(True, which='major') for ax in [ps_plot.axrte]]
            >>> ps_plot.update_plot()

        """

        self.fig.canvas.draw()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return ("Plots a pseudo section of TE and TM modes for data and "
                "response if given.")
#==============================================================================
# plot misfits as a pseudo-section
#==============================================================================


class PlotMisfitPseudoSection(object):
    """
    plot a pseudo section misfit of the data and response if given

    .. note:: the output file from winglink does not contain errors, so to get
              a normalized error, you need to input the error for each 
              component as a percent for resistivity and a value for phase
              and tipper.  If you used the data errors, unfortunately, you
              have to input those as arrays.  


    Arguments:
    -------------
        **wl_data_fn** : string
                         full path to output data file from winglink


    ==================== ==================================================
    key words            description
    ==================== ==================================================
    axmpte               matplotlib.axes instance for TE model phase
    axmptm               matplotlib.axes instance for TM model phase
    axmrte               matplotlib.axes instance for TE model app. res 
    axmrtm               matplotlib.axes instance for TM model app. res 
    axpte                matplotlib.axes instance for TE data phase 
    axptm                matplotlib.axes instance for TM data phase
    axrte                matplotlib.axes instance for TE data app. res.
    axrtm                matplotlib.axes instance for TM data app. res.
    cb_pad               padding between colorbar and axes
    cb_shrink            percentage to shrink the colorbar to
    fig                  matplotlib.figure instance
    fig_dpi              resolution of figure in dots per inch
    fig_num              number of figure instance
    fig_size             size of figure in inches (width, height)
    font_size            size of font in points
    label_list            list to label plots
    ml                   factor to label stations if 2 every other station
                         is labeled on the x-axis
    period               np.array of periods to plot
    phase_cmap           color map name of phase
    phase_limits_te      limits for te phase in degrees (min, max)
    phase_limits_tm      limits for tm phase in degrees (min, max)            
    plot_resp            [ 'y' | 'n' ] to plot response
    plot_yn              [ 'y' | 'n' ] 'y' to plot on instantiation

    res_cmap             color map name for resistivity
    res_limits_te        limits for te resistivity in log scale (min, max)
    res_limits_tm        limits for tm resistivity in log scale (min, max)
    rp_list               list of dictionaries as made from read2Dresp
    station_id           index to get station name (min, max)
    station_list          station list got from rp_list
    subplot_bottom       subplot spacing from bottom (relative coordinates) 
    subplot_hspace       vertical spacing between subplots
    subplot_left         subplot spacing from left  
    subplot_right        subplot spacing from right
    subplot_top          subplot spacing from top
    subplot_wspace       horizontal spacing between subplots
    ==================== ==================================================

    =================== =======================================================
    Methods             Description
    =================== =======================================================
    plot                plots a pseudo-section of apparent resistiviy and phase
                        of data and model if given.  called on instantiation 
                        if plot_yn is 'y'.
    redraw_plot         call redraw_plot to redraw the figures, 
                        if one of the attributes has been changed
    save_figure         saves the matplotlib.figure instance to desired 
                        location and format
    =================== =======================================================

   :Example: ::

        >>> import mtpy.modeling.occam2d as occam2d
        >>> ocd = occam2d.Occam2DData()
        >>> rfile = r"/home/Occam2D/Line1/Inv1/Test_15.resp"
        >>> ocd.data_fn = r"/home/Occam2D/Line1/Inv1/DataRW.dat"
        >>> ps1 = ocd.plot2PseudoSection(resp_fn=rfile) 

    """

    def __init__(self, data_fn, resp_fn, **kwargs):

        self.data_fn = data_fn
        self.resp_fn = resp_fn

        self.label_list = [r'$\rho_{TE}$', r'$\rho_{TM}$',
                           '$\phi_{TE}$', '$\phi_{TM}$',
                           '$\Re e\{T\}$', '$\Im m\{T\}$']

        self.phase_limits_te = kwargs.pop('phase_limits_te', (-10, 10))
        self.phase_limits_tm = kwargs.pop('phase_limits_tm', (-10, 10))
        self.res_limits_te = kwargs.pop('res_limits_te', (-2, 2))
        self.res_limits_tm = kwargs.pop('res_limits_tm', (-2, 2))
        self.tip_limits_re = kwargs.pop('tip_limits_re', (-.2, .2))
        self.tip_limits_im = kwargs.pop('tip_limits_im', (-.2, .2))

        self.phase_cmap = kwargs.pop('phase_cmap', 'BrBG')
        self.res_cmap = kwargs.pop('res_cmap', 'BrBG_r')
        self.tip_cmap = kwargs.pop('tip_cmap', 'PuOr')
        self.plot_tipper = kwargs.pop('plot_tipper', 'n')

        self.ml = kwargs.pop('ml', 2)
        self.station_id = kwargs.pop('station_id', [0, 4])

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)

        self.subplot_wspace = .0025
        self.subplot_hspace = .0
        self.subplot_right = .95
        self.subplot_left = .085
        self.subplot_top = .97
        self.subplot_bottom = .1

        self.font_size = kwargs.pop('font_size', 6)
        self.plot_yn = kwargs.pop('plot_yn', 'y')

        self.cb_shrink = .7
        self.cb_pad = .015

        self.axrte = None
        self.axrtm = None
        self.axpte = None
        self.axptm = None
        self.axtpr = None
        self.axtpi = None

        self.misfit_te_res = None
        self.misfit_te_phase = None
        self.misfit_tm_res = None
        self.misfit_tm_phase = None
        self.misfit_tip_real = None
        self.misfit_tip_imag = None

        self.fig = None
        self._data_obj = None

        if self.plot_yn == 'y':
            self.plot()

    def get_misfit(self):
        """
        compute misfit of MT response found from the model and the data.

        Need to normalize correctly
        """
        data_obj = Data()
        data_obj.read_data_file(self.data_fn)
        self._data_obj = data_obj

        resp_obj = Response()
        resp_obj.read_response_file(self.resp_fn)

        n_stations = len(data_obj.station_list)
        n_periods = len(data_obj.freq)

        self.misfit_te_res = np.zeros((n_periods, n_stations))
        self.misfit_te_phase = np.zeros((n_periods, n_stations))
        self.misfit_tm_res = np.zeros((n_periods, n_stations))
        self.misfit_tm_phase = np.zeros((n_periods, n_stations))
        self.misfit_tip_real = np.zeros((n_periods, n_stations))
        self.misfit_tip_imag = np.zeros((n_periods, n_stations))

        for rr, r_dict in zip(list(range(n_stations)), resp_obj.resp):
            self.misfit_te_res[:, rr] = r_dict['te_res'][1]
            self.misfit_tm_res[:, rr] = r_dict['tm_res'][1]
            self.misfit_te_phase[:, rr] = r_dict['te_phase'][1]
            self.misfit_tm_phase[:, rr] = r_dict['tm_phase'][1]
            self.misfit_tip_real[:, rr] = r_dict['re_tip'][1]
            self.misfit_tip_imag[:, rr] = r_dict['im_tip'][1]

        self.misfit_te_res = np.nan_to_num(self.misfit_te_res)
        self.misfit_te_phase = np.nan_to_num(self.misfit_te_phase)
        self.misfit_tm_res = np.nan_to_num(self.misfit_tm_res)
        self.misfit_tm_phase = np.nan_to_num(self.misfit_tm_phase)
        self.misfit_tip_real = np.nan_to_num(self.misfit_tip_real)
        self.misfit_tip_imag = np.nan_to_num(self.misfit_tip_imag)

    def plot(self):
        """
        plot pseudo section of data and response if given

        """

        self.get_misfit()

        ylimits = (self._data_obj.period.max(), self._data_obj.period.min())

        offset_list = np.append(self._data_obj.station_locations,
                                self._data_obj.station_locations[-1] * 1.15)

        # make a meshgrid for plotting
        # flip frequency so bottom corner is long period
        dgrid, fgrid = np.meshgrid(offset_list, self._data_obj.period[::-1])

        # make list for station labels
        ns = len(self._data_obj.station_list)
        sindex_1 = self.station_id[0]
        sindex_2 = self.station_id[1]
        slabel = [self._data_obj.station_list[ss][sindex_1:sindex_2]
                  for ss in range(0, ns, self.ml)]

        xloc = offset_list[0] + abs(offset_list[0] - offset_list[1]) / 5
        yloc = 1.10 * self._data_obj.period[1]

        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.hspace'] = self.subplot_hspace
        plt.rcParams['figure.subplot.wspace'] = self.subplot_wspace

        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()

        if self.plot_tipper != 'y':
            self.axrte = self.fig.add_subplot(2, 2, 1)
            self.axrtm = self.fig.add_subplot(2, 2, 2, sharex=self.axrte)
            self.axpte = self.fig.add_subplot(2, 2, 3, sharex=self.axrte)
            self.axptm = self.fig.add_subplot(2, 2, 4, sharex=self.axrte)

        else:
            self.axrte = self.fig.add_subplot(2, 3, 1)
            self.axrtm = self.fig.add_subplot(2, 3, 2, sharex=self.axrte)
            self.axpte = self.fig.add_subplot(2, 3, 4, sharex=self.axrte)
            self.axptm = self.fig.add_subplot(2, 3, 5, sharex=self.axrte)
            self.axtpr = self.fig.add_subplot(2, 3, 3, sharex=self.axrte)
            self.axtpi = self.fig.add_subplot(2, 3, 6, sharex=self.axrte)

        #--> TE Resistivity
        self.axrte.pcolormesh(dgrid,
                              fgrid,
                              np.flipud(self.misfit_te_res),
                              cmap=self.res_cmap,
                              vmin=self.res_limits_te[0],
                              vmax=self.res_limits_te[1])
        #--> TM Resistivity
        self.axrtm.pcolormesh(dgrid,
                              fgrid,
                              np.flipud(self.misfit_tm_res),
                              cmap=self.res_cmap,
                              vmin=self.res_limits_tm[0],
                              vmax=self.res_limits_tm[1])
        #--> TE Phase
        self.axpte.pcolormesh(dgrid,
                              fgrid,
                              np.flipud(self.misfit_te_phase),
                              cmap=self.phase_cmap,
                              vmin=self.phase_limits_te[0],
                              vmax=self.phase_limits_te[1])
        #--> TM Phase
        self.axptm.pcolormesh(dgrid,
                              fgrid,
                              np.flipud(self.misfit_tm_phase),
                              cmap=self.phase_cmap,
                              vmin=self.phase_limits_tm[0],
                              vmax=self.phase_limits_tm[1])

        ax_list = [self.axrte, self.axrtm, self.axpte, self.axptm]

        if self.plot_tipper == 'y':
            self.axtpr.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(self.misfit_tip_real),
                                  cmap=self.tip_cmap,
                                  vmin=self.tip_limits_re[0],
                                  vmax=self.tip_limits_re[1])
            self.axtpi.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(self.misfit_tip_imag),
                                  cmap=self.tip_cmap,
                                  vmin=self.tip_limits_im[0],
                                  vmax=self.tip_limits_im[1])

            ax_list.append(self.axtpr)
            ax_list.append(self.axtpi)
         # make everthing look tidy
        for xx, ax in enumerate(ax_list):
            ax.semilogy()
            ax.set_ylim(ylimits)
            ax.xaxis.set_ticks(offset_list[np.arange(0, ns, self.ml)])
            ax.xaxis.set_ticks(offset_list, minor=True)
            ax.xaxis.set_ticklabels(slabel)
            ax.set_xlim(offset_list.min(), offset_list.max())
            cbx = mcb.make_axes(ax,
                                shrink=self.cb_shrink,
                                pad=self.cb_pad)

            # te res
            if xx == 0:
                plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                cb = mcb.ColorbarBase(cbx[0], cmap=self.res_cmap,
                                      norm=Normalize(vmin=self.res_limits_te[0],
                                                     vmax=self.res_limits_te[1]))
            # tm res
            elif xx == 1:
                plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                cb = mcb.ColorbarBase(cbx[0], cmap=self.res_cmap,
                                      norm=Normalize(vmin=self.res_limits_tm[0],
                                                     vmax=self.res_limits_tm[1]))
                cb.set_label('Log$_{10}$ App. Res. ($\Omega \cdot$m)',
                             fontdict={'size': self.font_size + 1,
                                       'weight': 'bold'})
            # te phase
            elif xx == 2:
                cb = mcb.ColorbarBase(cbx[0], cmap=self.phase_cmap,
                                      norm=Normalize(vmin=self.phase_limits_te[0],
                                                     vmax=self.phase_limits_te[1]))
            # tm phase
            elif xx == 3:
                plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                cb = mcb.ColorbarBase(cbx[0], cmap=self.phase_cmap,
                                      norm=Normalize(vmin=self.phase_limits_tm[0],
                                                     vmax=self.phase_limits_tm[1]))
                cb.set_label('Phase (deg)',
                             fontdict={'size': self.font_size + 1,
                                       'weight': 'bold'})

            # real tipper
            elif xx == 4:
                plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                cb = mcb.ColorbarBase(cbx[0], cmap=self.tip_cmap,
                                      norm=Normalize(vmin=self.tip_limits_re[0],
                                                     vmax=self.tip_limits_re[1]))
                cb.set_label('Re{Tip}',
                             fontdict={'size': self.font_size + 1,
                                       'weight': 'bold'})
            # imag tipper
            elif xx == 5:
                plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                cb = mcb.ColorbarBase(cbx[0], cmap=self.tip_cmap,
                                      norm=Normalize(vmin=self.tip_limits_im[0],
                                                     vmax=self.tip_limits_im[1]))
                cb.set_label('Im{Tip}',
                             fontdict={'size': self.font_size + 1,
                                       'weight': 'bold'})

            # make label for plot
            ax.text(xloc, yloc, self.label_list[xx],
                    fontdict={'size': self.font_size + 2},
                    bbox={'facecolor': 'white'},
                    horizontalalignment='left',
                    verticalalignment='top')

            if xx == 0 or xx == 2:
                ax.set_ylabel('Period (s)',
                              fontdict={'size': self.font_size + 2,
                                        'weight': 'bold'})
            if xx > 1:
                ax.set_xlabel('Station', fontdict={'size': self.font_size + 2,
                                                   'weight': 'bold'})

        plt.show()

    def redraw_plot(self):
        """
        redraw plot if parameters were changed

        use this function if you updated some attributes and want to re-plot.

        :Example: ::

            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotPseudoSection()
            >>> #change color of te markers to a gray-blue
            >>> p1.res_cmap = 'seismic_r'
            >>> p1.redraw_plot()
        """

        plt.close(self.fig)
        self.plot()

    def save_figure(self, save_fn, file_format='pdf', orientation='portrait',
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
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> ps1.save_plot(r'/home/MT/figures', file_format='jpg')

        """

        if fig_dpi == None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        else:
            save_fn = os.path.join(save_fn, 'OccamMisfitPseudoSection.' +
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

    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.  

        Use this if you change an of the .fig or axes properties

        :Example: ::

            >>> # to change the grid lines to only be on the major ticks
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> [ax.grid(True, which='major') for ax in [ps1.axrte,ps1.axtep]]
            >>> ps1.update_plot()

        """

        self.fig.canvas.draw()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return ("Plots a pseudo section of TE and TM modes for data and "
                "response if given.")
