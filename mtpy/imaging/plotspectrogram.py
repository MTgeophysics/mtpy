# -*- coding: utf-8 -*-
"""
==================
plotft
==================

    *PlotTF --> will plot a time frequency distribution of your choice


Created on Mon Aug 19 16:24:29 2013

@author: jpeacock
"""

#=================================================================

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import mtpy.processing.tf as mttf
import os

import  mtpy.utils.exceptions as mtex

#=================================================================


class PlotTF(object):
    """
    class to plot Time-Frequency


    """

    def __init__(self, time_series, tf_type='smethod', **kwargs):

        self.time_series = time_series
        self.tf_type = tf_type

        self.tf_array = None
        self.time_list = None
        self.freq_list = None

        self.tf_nh = kwargs.pop('nh', None)
        self.tf_ng = kwargs.pop('ng', None)
        self.tf_tstep = kwargs.pop('tstep', 2**5)
        self.tf_nfbins = kwargs.pop('nfbins', 2**9)
        self.tf_L = kwargs.pop('L', 11)
        self.tf_beta = kwargs.pop('beta', 0.2)
        self.tf_alpha = kwargs.pop('alpha', None)
        self.tf_sigmat = kwargs.pop('sigmat', None)
        self.tf_sigmaf = kwargs.pop('sigmaf', None)
        self.tf_sigmaL = kwargs.pop('sigmaL', None)
        self.tf_ngwv = kwargs.pop('ngwv', None)
        self.tf_nhwv = kwargs.pop('nhwv', None)
        self.tf_thresh = kwargs.pop('thresh', None)
        self.tf_robust_type = kwargs.pop('robusttype', 'median')


        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_dpi = kwargs.pop('fig_dpi', 300)
        self.fig_size = kwargs.pop('fig_size', [6,6])

        self.font_size = kwargs.pop('font_size', 7)

        self.df = kwargs.pop('df', 1.0)

        self.start_time = kwargs.pop('start_time', 0)
        self.time_units = kwargs.pop('time_units', 'hrs')

        self.tf_scale = kwargs.pop('tf_scale', 'log')

        self.freq_scale = kwargs.pop('freq_scale', 'log')
        self.freq_units = kwargs.pop('freq_units', 'hz')

        self.cmap = kwargs.pop('cmap', 'jet')
        self.climits = kwargs.pop('climits', None)

        self.plot_title = kwargs.pop('title', None)
        self.plot_interpolation = kwargs.pop('plot_interpolation', 'gaussian')
        self.plot_aspect_ratio = kwargs.pop('plot_aspect_ratio', 'auto')
        self.plot_type = kwargs.pop('plot_type', 'tf')
        self.plot_normalize = kwargs.pop('plot_normalize', 'n')

        self.cb_orientation = kwargs.pop('cb_orientation', 'vertical')
        self.cb_shrink = kwargs.pop('cb_shrink', .8)
        self.cb_aspect_ratio = kwargs.pop('cb_aspect_ratio', 20)
        self.cb_pad = kwargs.pop('cb_pad', .05)

        self.cb = None
        self.fig = None
        self.axtf = None
        self.axts = None
        self.axps = None
        self.fig_fn = None

        self.x_major_tick = None
        self.x_minor_tick = None

        self.subplot_left = 0.12
        self.subplot_right = 0.99
        self.subplot_bottom = 0.12
        self.subplot_top = 0.96
        self.subplot_wspace = 0.25
        self.subplot_hspace = 0.20

        self.lw = 0.75
        self.line_color_ts = 'k'
        self.line_color_ps = 'k'

        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()

    def _get_tf(self):
        """
        get the specified time frequency distribution
        """
        #--> short time fourier transform
        if self.tf_type == 'stft':
            if self.tf_nh == None:
                self.tf_nh = 2**8
            if self.tf_ng == None:
                self.tf_ng = 1

            kwargs = {'nh':self.tf_nh,
                      'tstep': self.tf_tstep,
                      'ng':self.tf_ng,
                      'nfbins':self.tf_nfbins,
                      'df':self.df}

            tf_tuple = mttf.stft(self.time_series, **kwargs)
            self.tf_array = tf_tuple[0]
            self.time_list = tf_tuple[1]
            self.freq_list = tf_tuple[2]

       #--> reassigned stft
        elif self.tf_type == 'reassigned_stft':
            if self.tf_nh == None:
                self.tf_nh = 2**6-1
            if self.tf_alpha == None:
                self.tf_alpha = 4.0
            if self.tf_thresh == None:
                pass

            kwargs = {'nh':self.tf_nh,
                      'tstep': self.tf_tstep,
                      'nfbins':self.tf_nfbins,
                      'alpha':self.tf_alpha,
                      'threshold':self.tf_thresh,
                      'df':self.df}

            tf_tuple = mttf.reassigned_stft(self.time_series, **kwargs)
            self.tf_array = tf_tuple[0]
            self.time_list = tf_tuple[1]
            self.freq_list = tf_tuple[2]

        #--> Wigner-ville distribution
        elif self.tf_type == 'wvd':
            if self.tf_nh == None:
                self.tf_nh = 2**8-1

            kwargs = {'nh':self.tf_nh,
                      'tstep': self.tf_tstep,
                      'nfbins':self.tf_nfbins,
                      'df':self.df}

            tf_tuple = mttf.wvd(self.time_series, **kwargs)
            self.tf_array = tf_tuple[0]
            self.time_list = tf_tuple[1]
            self.freq_list = tf_tuple[2]

        #--> smoothe pseudo wigner-ville distribution
        elif self.tf_type == 'spwvd':
            kwargs = {'nh':self.tf_nh,
                      'ng':self.tf_ng,
                      'sigmat':self.tf_sigmat,
                      'sigmaf':self.tf_sigmaf,
                      'tstep': self.tf_tstep,
                      'nfbins':self.tf_nfbins,
                      'df':self.df}

            tf_tuple = mttf.spwvd(self.time_series, **kwargs)
            self.tf_array = tf_tuple[0]
            self.time_list = tf_tuple[1]
            self.freq_list = tf_tuple[2]

        #--> robust wigner ville-distribution
        elif self.tf_type == 'robust_wvd':
            if self.tf_nh == None:
                self.tf_nh = 2**7-1
            if self.tf_ng == None:
                self.tf_ng = 2**4-1

            kwargs = {'nh':self.tf_nh,
                      'ng':self.tf_ng,
                      'sigmat':self.tf_sigmat,
                      'sigmaf':self.tf_sigmaf,
                      'tstep': self.tf_tstep,
                      'nfbins':self.tf_nfbins,
                      'df':self.df}
            tf_tuple = mttf.robust_wvd(self.time_series, **kwargs)
            self.tf_array = tf_tuple[0]
            self.time_list = tf_tuple[1]
            self.freq_list = tf_tuple[2]

        #--> robust wigner ville-distribution
        elif self.tf_type == 'specwv':
            if self.tf_nh == None:
                self.tf_nh = 2**8
            if self.tf_nhwv == None:
                self.tf_nhwv = 2**9-1
            if self.tf_ngwv == None:
                self.tf_ngwv = 2**3-1

            kwargs = {'nhs':self.tf_nh,
                      'nhwv':self.tf_nh,
                      'ngwv':self.tf_ng,
                      'sigmat':self.tf_sigmat,
                      'sigmaf':self.tf_sigmaf,
                      'tstep': self.tf_tstep,
                      'nfbins':self.tf_nfbins,
                      'df':self.df}

            tf_tuple = mttf.specwv(self.time_series, **kwargs)
            self.tf_array = tf_tuple[0]
            self.time_list = tf_tuple[1]
            self.freq_list = tf_tuple[2]

        #--> modified b
        elif self.tf_type == 'modifiedb':
            if self.tf_nh == None:
                self.tf_nh = 2**8-1

            kwargs = {'nh':self.tf_nh,
                      'beta':self.tf_beta,
                      'tstep': self.tf_tstep,
                      'nfbins':self.tf_nfbins,
                      'df':self.df}

            tf_tuple = mttf.modifiedb(self.time_series, **kwargs)
            self.tf_array = tf_tuple[0]
            self.time_list = tf_tuple[1]
            self.freq_list = tf_tuple[2]

        #--> robust stft with vector median filter
        elif self.tf_type == 'robust_stft_median':
            if self.tf_nh == None:
                self.tf_nh = 2**8

            kwargs = {'nh':self.tf_nh,
                      'tstep': self.tf_tstep,
                      'nfbins':self.tf_nfbins,
                      'df':self.df}

            tf_tuple = mttf.robust_stft_median(self.time_series, **kwargs)
            self.tf_array = tf_tuple[0]
            self.time_list = tf_tuple[1]
            self.freq_list = tf_tuple[2]

        #--> robust stft with L-distribution
        elif self.tf_type == 'robust_stft_L':
            if self.tf_nh == None:
                self.tf_nh = 2**8
            if self.tf_alpha == None:
                self.tf_alpha = 0.325

            kwargs = {'nh':self.tf_nh,
                      'alpha':self.tf_alpha,
                      'tstep': self.tf_tstep,
                      'nfbins':self.tf_nfbins,
                      'df':self.df}

            tf_tuple = mttf.robust_stft_L(self.time_series, **kwargs)
            self.tf_array = tf_tuple[0]
            self.time_list = tf_tuple[1]
            self.freq_list = tf_tuple[2]

        #--> smethod
        elif self.tf_type == 'smethod':
            if self.tf_nh == None:
                self.tf_nh = 2**8
            if self.tf_ng == None:
                self.tf_ng = 1
            if self.tf_alpha == None:
                self.tf_alpha = 0.325

            kwargs = {'nh':self.tf_nh,
                      'L':self.tf_L,
                      'tstep': self.tf_tstep,
                      'nfbins':self.tf_nfbins,
                      'sigmaL':self.tf_sigmaL,
                      'df':self.df}

            tf_tuple = mttf.smethod(self.time_series, **kwargs)
            self.tf_array = tf_tuple[0]
            self.time_list = tf_tuple[1]
            self.freq_list = tf_tuple[2]

        #--> robust smethod
        elif self.tf_type == 'robust_smethod':
            if self.tf_nh == None:
                self.tf_nh = 2**8
            if self.tf_ng == None:
                self.tf_ng = 1
            if self.tf_alpha == None:
                self.tf_alpha = 0.325


            kwargs = {'nh':self.tf_nh,
                      'L':self.tf_L,
                      'alpha':self.tf_alpha,
                      'tstep': self.tf_tstep,
                      'nfbins':self.tf_nfbins,
                      'sigmaL':self.tf_sigmaL,
                      'robusttype':self.tf_robust_type,
                      'df':self.df}

            tf_tuple = mttf.robust_smethod(self.time_series, **kwargs)
            self.tf_array = tf_tuple[0]
            self.time_list = tf_tuple[1]
            self.freq_list = tf_tuple[2]

        #--> reassigned smethod
        elif self.tf_type == 'reassigned_smethod':
            if self.tf_nh == None:
                self.tf_nh = 2**8-1
            if self.tf_ng == None:
                self.tf_ng = 1
            if self.tf_alpha == None:
                self.tf_alpha = 4.0
            if self.tf_thresh == None:
                self.tf_thresh = .01


            kwargs = {'nh':self.tf_nh,
                      'L':self.tf_L,
                      'alpha':self.tf_alpha,
                      'tstep': self.tf_tstep,
                      'nfbins':self.tf_nfbins,
                      'threshold':self.tf_thresh,
                      'robusttype':self.tf_robust_type,
                      'df':self.df}

            tf_tuple = mttf.reassigned_smethod(self.time_series, **kwargs)
            self.tf_array = tf_tuple[0]
            self.time_list = tf_tuple[1]
            self.freq_list = tf_tuple[2]

        else:
            raise mtex.MTpyError_inputarguments('{0}'.format(self.tf_type)+
                        ' is not definded see mtpy.processing.tf for options')

        #print information for user
        print('{0} tf parameters {0}'.format('-'*5))
        for kw in sorted(kwargs.keys()):
            print('{0}{1} = {2}'.format(' '*4, kw, kwargs[kw]))

    def plot(self):
        """
        plot the time frequency distribution

        """
        #get the requested time-frequency distribution
        self._get_tf()

        #time increment
        if self.time_units == 'hrs':
            tinc = 3600*self.df
            if self.x_major_tick == None:
                x_major_tick = 1
            if self.x_minor_tick == None:
                x_minor_tick = .15
        elif self.time_units == 'min':
            tinc = 60*self.df
            if self.x_major_tick == None:
                x_major_tick = 5
            if self.x_minor_tick == None:
                x_minor_tick = 1
        elif self.time_units == 'sec':
            tinc = 1*self.df
            if self.x_major_tick == None:
                x_major_tick = 60
            if self.x_minor_tick == None:
                x_minor_tick = 15
        else:
            raise mtex.MTpyError_inputarguments('{0} is not defined'.format(
                                                self.time_units))

        #scale time-frequency
        if self.tf_scale == 'log':
            self.tf_array[np.where(abs(self.tf_array)==0)] = 1.0
            if self.plot_normalize == 'y':
                plottfarray = 10*np.log10(abs(self.tf_array/
                                                np.max(abs(self.tf_array))))
            else:
                plottfarray = 10*np.log10(abs(self.tf_array))
        elif self.tf_scale == 'linear':
            if self.plot_normalize == 'y':
                plottfarray = abs(self.tf_array/np.max(abs(self.tf_array)))
            else:
                plottfarray = abs(self.tf_array)

        #period or frequency
        if self.freq_units == 'y':
            self.freq_list[1:] = 1./self.freq_list[1:]
            self.freq_list[0] = 2*self.freq_list[1]
        elif self.freq_units == 'n':
            pass

        #set properties for the plot
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        plt.rcParams['figure.subplot.wspace'] = self.subplot_wspace
        plt.rcParams['figure.subplot.hspace'] = self.subplot_hspace

        #set the font dictionary
        fdict={'size':self.font_size+2, 'weight':'bold'}

        #make a meshgrid if yscale is logarithmic
        if self.freq_scale == 'log':
            logt, logf = np.meshgrid(self.time_list/tinc, self.freq_list)

        #make figure
        self.fig = plt.figure(self.fig_num,self.fig_size, dpi=self.fig_dpi)
        self.fig.clf()

        if self.plot_type == 'all':
            self.axps = self.fig.add_axes([.05, .25, .1, .7])
            self.axts = self.fig.add_axes([.25, .05, .60, .1])
            self.axtf = self.fig.add_axes([.25, .25, .75, .7])

            #plot time series
            st = self.start_time
            time_array = np.arange(st,
                                   st+self.time_series.size/self.df,
                                   1./self.df)

            self.axts.plot(time_array,
                           self.time_series,
                           color=self.line_color_ts,
                           lw=self.lw)
            self.axts.axis('tight')

            FX = np.fft.fft(mttf.padzeros(self.time_series))
            FXfreq = np.fft.fftfreq(len(FX), 1./self.df)

            #plot power spectra
            if self.freq_scale == 'log':
                self.axps.loglog(abs(FX[0:len(FX)/2]/max(abs(FX))),
                                   FXfreq[0:len(FX)/2],
                                          color=self.line_color_ps,
                                          lw=self.lw)
            else:
                self.axps.semilogx(abs(FX[0:len(FX)/2]/max(abs(FX))),
                                   FXfreq[0:len(FX)/2],
                                   color=self.line_color_ps,
                                   lw=self.lw)
            self.axps.axis('tight')
            self.axps.set_ylim(self.freq_list[1],self.freq_list[-1])
        else:
            self.axtf = self.fig.add_subplot(1, 1, 1,
                                         aspect=self.plot_aspect_ratio)

        #--> get color limits
        if self.climits != None:
            vmin = self.climits[0]
            vmax = self.climits[1]
        else:
            vmin = plottfarray.min()
            vmax = plottfarray.max()


        #add in log yscale
        if self.freq_scale == 'log':
            #need to flip the matrix so that origin is bottom right
            cbp = self.axtf.pcolormesh(logt,
                                       logf,
                                       np.flipud(plottfarray),
                                       cmap=self.cmap,
                                       vmin=vmin,
                                       vmax=vmax)
            self.axtf.semilogy()
            self.axtf.set_ylim(self.freq_list[1],self.freq_list[-1])

            self.axtf.set_xlim(logt.min(), logt.max())
            self.cb = plt.colorbar(cbp,
                                   orientation=self.cb_orientation,
                                   shrink=self.cb_shrink,
                                   pad=self.cb_pad,
                                   aspect=self.cb_aspect_ratio,
                                   use_gridspec=True)
        else:
            cbp = self.axtf.imshow(plottfarray,
                                 extent=(self.time_list[0]/tinc+self.start_time,
                                         self.time_list[-1]/tinc+self.start_time,
                                         self.freq_list[1],self.freq_list[-1]),
                                aspect=self.plot_aspect_ratio,
                                vmin=vmin,
                                vmax=vmax,
                                cmap=self.cmap,
                                interpolation=self.plot_interpolation)

            self.cb = plt.colorbar(orientation=self.cb_orientation,
                                   shrink=self.cb_shrink,
                                   pad=self.cb_pad,
                                   aspect=self.cb_aspect_ratio,
                                   use_gridspec=True)

        #--> make the plot look nice
        self.axtf.set_xlabel('time({0})'.format(self.time_units),
                           fontdict=fdict)
        self.axtf.xaxis.set_major_locator(MultipleLocator(x_major_tick))
        self.axtf.xaxis.set_minor_locator(MultipleLocator(x_minor_tick))

        if self.freq_units == 's':
            self.axtf.set_ylabel('period (s)', fontdict=fdict)
        else:
            self.axtf.set_ylabel('frequency (Hz)', fontdict=fdict)
        if self.plot_title != None:
            self.axtf.set_title(self.plot_title, fontdict=fdict)

        plt.show()

    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.

        Use this if you change an of the .fig or axes properties

        :Example: ::

            >>> # to change the grid lines to be on the major ticks and gray
            >>> tf1.ax.grid(True, which='major', color=(.5,.5,.5))
            >>> tf1.update_plot()

        """

        self.fig.canvas.draw()

    def redraw_plot(self):
        """
        use this function if you updated some attributes and want to re-plot.

        :Example: ::
            >>> tf1.plot_type = 'tf'
            >>> tf1.redraw_plot()
        """

        plt.close(self.fig)
        self.plot()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return "Plots time frequency distribution"

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
                            save_fn/TF_tftype.file_format

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
            >>> tf1.save_plot(r'/home/MT/figures', file_format='jpg')

        """

        if fig_dpi == None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        else:
            save_fn = os.path.join(save_fn, 'TF_{0}.'.format(self.tf_type)+
                                   file_format)
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        if close_plot == 'y':
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print('Saved figure to: {0}'.format(self.fig_fn))
