"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch

"""
import os

import numpy as np
from matplotlib import colors as colors, pyplot as plt, colorbar as mcb
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from mtpy.modeling.modem import Data, Residual

__all__ = ['PlotRMSMaps']


class PlotRMSMaps(object):
    """
    plots the RMS as (data-model)/(error) in map view for all components
    of the data file.  Gets this infomration from the .res file output
    by ModEM.

    Arguments:
    ------------------

        **residual_fn** : string
                          full path to .res file

    =================== =======================================================
    Attributes                   Description
    =================== =======================================================
    fig                 matplotlib.figure instance for a single plot
    fig_dpi             dots-per-inch resolution of figure *default* is 200
    fig_num             number of fig instance *default* is 1
    fig_size            size of figure in inches [width, height]
                        *default* is [7,6]
    font_size           font size of tick labels, axis labels are +2
                        *default* is 8
    marker              marker style for station rms,
                        see matplotlib.line for options,
                        *default* is 's' --> square
    marker_size         size of marker in points. *default* is 10
    pad_x               padding in map units from edge of the axis to stations
                        at the extremeties in longitude.
                        *default* is 1/2 tick_locator
    pad_y               padding in map units from edge of the axis to stations
                        at the extremeties in latitude.
                        *default* is 1/2 tick_locator
    period_index        index of the period you want to plot according to
                        self.residual.period_list. *default* is 1
    plot_yn             [ 'y' | 'n' ] default is 'y' to plot on instantiation
    plot_z_list         internal variable for plotting
    residual            modem.Data instance that holds all the information
                        from the residual_fn given
    residual_fn         full path to .res file
    rms_cmap            matplotlib.cm object for coloring the markers
    rms_cmap_dict       dictionary of color values for rms_cmap
    rms_max             maximum rms to plot. *default* is 5.0
    rms_min             minimum rms to plot. *default* is 1.0
    save_path           path to save figures to. *default* is directory of
                        residual_fn
    subplot_bottom      spacing from axis to bottom of figure canvas.
                        *default* is .1
    subplot_hspace      horizontal spacing between subplots.
                        *default* is .1
    subplot_left        spacing from axis to left of figure canvas.
                        *default* is .1
    subplot_right       spacing from axis to right of figure canvas.
                        *default* is .9
    subplot_top         spacing from axis to top of figure canvas.
                        *default* is .95
    subplot_vspace      vertical spacing between subplots.
                        *default* is .01
    tick_locator        increment for x and y major ticks. *default* is
                        limits/5
    =================== =======================================================

    =================== =======================================================
    Methods             Description
    =================== =======================================================
    plot                plot rms maps for a single period
    plot_loop           loop over all frequencies and save figures to save_path
    read_residual_fn    read in residual_fn
    redraw_plot         after updating attributes call redraw_plot to
                        well redraw the plot
    save_figure         save the figure to a file
    =================== =======================================================


    :Example: ::

        >>> import mtpy.modeling.modem as modem
        >>> rms_plot = PlotRMSMaps(r"/home/ModEM/Inv1/mb_NLCG_030.res")
        >>> # change some attributes
        >>> rms_plot.fig_size = [6, 4]
        >>> rms_plot.rms_max = 3
        >>> rms_plot.redraw_plot()
        >>> # happy with the look now loop over all periods
        >>> rms_plot.plot_loop()
    """

    def __init__(self, residual_fn, **kwargs):
        self.residual_fn = residual_fn
        self.residual = None
        self.save_path = kwargs.pop('save_path', os.path.dirname(self.residual_fn))

        self.period_index = kwargs.pop('period_index', 0)

        self.subplot_left = kwargs.pop('subplot_left', .1)
        self.subplot_right = kwargs.pop('subplot_right', .9)
        self.subplot_top = kwargs.pop('subplot_top', .95)
        self.subplot_bottom = kwargs.pop('subplot_bottom', .1)
        self.subplot_hspace = kwargs.pop('subplot_hspace', .1)
        self.subplot_vspace = kwargs.pop('subplot_vspace', .01)

        self.font_size = kwargs.pop('font_size', 8)

        self.fig_size = kwargs.pop('fig_size', [7.75, 6.75])
        self.fig_dpi = kwargs.pop('fig_dpi', 200)
        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig = None

        self.marker = kwargs.pop('marker', 's')
        self.marker_size = kwargs.pop('marker_size', 10)

        self.rms_max = kwargs.pop('rms_max', 5)
        self.rms_min = kwargs.pop('rms_min', 0)

        self.tick_locator = kwargs.pop('tick_locator', None)
        self.pad_x = kwargs.pop('pad_x', None)
        self.pad_y = kwargs.pop('pad_y', None)

        self.plot_yn = kwargs.pop('plot_yn', 'y')

        # colormap for rms, goes white to black from 0 to rms max and
        # red below 1 to show where the data is being over fit
        self.rms_cmap_dict = {'red': ((0.0, 1.0, 1.0),
                                      (0.2, 1.0, 1.0),
                                      (1.0, 0.0, 0.0)),
                              'green': ((0.0, 0.0, 0.0),
                                        (0.2, 1.0, 1.0),
                                        (1.0, 0.0, 0.0)),
                              'blue': ((0.0, 0.0, 0.0),
                                       (0.2, 1.0, 1.0),
                                       (1.0, 0.0, 0.0))}

        self.rms_cmap = colors.LinearSegmentedColormap('rms_cmap',
                                                       self.rms_cmap_dict,
                                                       256)

        self.plot_z_list = [{'label': r'$Z_{xx}$', 'index': (0, 0), 'plot_num': 1},
                            {'label': r'$Z_{xy}$', 'index': (0, 1), 'plot_num': 2},
                            {'label': r'$Z_{yx}$', 'index': (1, 0), 'plot_num': 3},
                            {'label': r'$Z_{yy}$', 'index': (1, 1), 'plot_num': 4},
                            {'label': r'$T_{x}$', 'index': (0, 0), 'plot_num': 5},
                            {'label': r'$T_{y}$', 'index': (0, 1), 'plot_num': 6}]

        if self.plot_yn == 'y':
            self.plot()

    def read_residual_fn(self):
        if self.residual is None:
            self.residual = Residual(residual_fn=self.residual_fn)
#            self.residual.read_data_file(self.residual_fn)
            self.residual.read_residual_file()
            self.residual.get_rms()
        else:
            pass

    def plot(self):
        """
        plot rms in map view
        """

        self.read_residual_fn()

        font_dict = {'size': self.font_size + 2, 'weight': 'bold'}
        rms_1 = 1. / self.rms_max

        if self.tick_locator is None:
            x_locator = np.round((self.residual.residual_array['lon'].max() -
                                  self.residual.residual_array['lon'].min()) / 5, 2)
            y_locator = np.round((self.residual.residual_array['lat'].max() -
                                  self.residual.residual_array['lat'].min()) / 5, 2)

            if x_locator > y_locator:
                self.tick_locator = x_locator

            elif x_locator < y_locator:
                self.tick_locator = y_locator

        if self.pad_x is None:
            self.pad_x = self.tick_locator / 2
        if self.pad_y is None:
            self.pad_y = self.tick_locator / 2

        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        plt.rcParams['figure.subplot.wspace'] = self.subplot_hspace
        plt.rcParams['figure.subplot.hspace'] = self.subplot_vspace
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)

        for p_dict in self.plot_z_list:
            ax = self.fig.add_subplot(3, 2, p_dict['plot_num'], aspect='equal')

            ii = p_dict['index'][0]
            jj = p_dict['index'][0]

#            for r_arr in self.residual.residual_array:
            for ridx in range(len(self.residual.residual_array)):
                
                if self.period_index == 'all':
                    r_arr = self.residual.rms_array[ridx]
                    if p_dict['plot_num'] < 5:
                        rms = r_arr['rms_z']
                    else:
                        rms = r_arr['rms_tip']
                else:
                    r_arr = self.residual.residual_array[ridx]
                    
                    # calulate the rms self.residual/error
                    if p_dict['plot_num'] < 5:
                        rms = r_arr['z'][self.period_index, ii, jj].__abs__() / \
                              r_arr['z_err'][self.period_index, ii, jj].real
    
                    else:
                        rms = r_arr['tip'][self.period_index, ii, jj].__abs__() / \
                              r_arr['tip_err'][self.period_index, ii, jj].real

                # color appropriately
                if np.nan_to_num(rms) == 0.0:
                    marker_color = (1, 1, 1)
                    marker = '.'
                    marker_size = .1
                    marker_edge_color = (1, 1, 1)
                if rms > self.rms_max:
                    marker_color = (0, 0, 0)
                    marker = self.marker
                    marker_size = self.marker_size
                    marker_edge_color = (0, 0, 0)

                elif 1 <= rms <= self.rms_max:
                    r_color = 1 - rms / self.rms_max + rms_1
                    marker_color = (r_color, r_color, r_color)
                    marker = self.marker
                    marker_size = self.marker_size
                    marker_edge_color = (0, 0, 0)

                elif rms < 1:
                    r_color = 1 - rms / self.rms_max
                    marker_color = (1, r_color, r_color)
                    marker = self.marker
                    marker_size = self.marker_size
                    marker_edge_color = (0, 0, 0)

                ax.plot(r_arr['lon'], r_arr['lat'],
                        marker=marker,
                        ms=marker_size,
                        mec=marker_edge_color,
                        mfc=marker_color,
                        zorder=3)

            if p_dict['plot_num'] == 1 or p_dict['plot_num'] == 3:
                ax.set_ylabel('Latitude (deg)', fontdict=font_dict)
                plt.setp(ax.get_xticklabels(), visible=False)

            elif p_dict['plot_num'] == 2 or p_dict['plot_num'] == 4:
                plt.setp(ax.get_xticklabels(), visible=False)
                plt.setp(ax.get_yticklabels(), visible=False)

            elif p_dict['plot_num'] == 6:
                plt.setp(ax.get_yticklabels(), visible=False)
                ax.set_xlabel('Longitude (deg)', fontdict=font_dict)

            else:
                ax.set_xlabel('Longitude (deg)', fontdict=font_dict)
                ax.set_ylabel('Latitude (deg)', fontdict=font_dict)

            ax.text(self.residual.residual_array['lon'].min() + .005 - self.pad_x,
                    self.residual.residual_array['lat'].max() - .005 + self.pad_y,
                    p_dict['label'],
                    verticalalignment='top',
                    horizontalalignment='left',
                    bbox={'facecolor': 'white'},
                    zorder=3)

            ax.tick_params(direction='out')
            ax.grid(zorder=0, color=(.75, .75, .75))

            # [line.set_zorder(3) for line in ax.lines]

            ax.set_xlim(self.residual.residual_array['lon'].min() - self.pad_x,
                        self.residual.residual_array['lon'].max() + self.pad_x)

            ax.set_ylim(self.residual.residual_array['lat'].min() - self.pad_y,
                        self.residual.residual_array['lat'].max() + self.pad_y)

            ax.xaxis.set_major_locator(MultipleLocator(self.tick_locator))
            ax.yaxis.set_major_locator(MultipleLocator(self.tick_locator))
            ax.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%2.2f'))

        # cb_ax = mcb.make_axes(ax, orientation='vertical', fraction=.1)
        cb_ax = self.fig.add_axes([self.subplot_right + .02, .225, .02, .45])
        color_bar = mcb.ColorbarBase(cb_ax,
                                     cmap=self.rms_cmap,
                                     norm=colors.Normalize(vmin=self.rms_min,
                                                           vmax=self.rms_max),
                                     orientation='vertical')

        color_bar.set_label('RMS', fontdict=font_dict)
        if self.period_index == 'all':
            self.fig.suptitle('all periods',
                              fontdict={'size': self.font_size + 3, 'weight': 'bold'})
        else:            
            self.fig.suptitle('period = {0:.5g} (s)'.format(self.residual.period_list[self.period_index]),
                              fontdict={'size': self.font_size + 3, 'weight': 'bold'})
        self.fig.show()

    def redraw_plot(self):
        plt.close(self.fig)
        self.plot()

    def save_figure(self, save_path=None, save_fn_basename=None,
                    save_fig_dpi=None, fig_format='png', fig_close=True):
        """
        save figure in the desired format
        """
        if save_path is not None:
            self.save_path = save_path

        if save_fn_basename is not None:
            pass
        else:
            if self.period_index == 'all':
                save_fn_basename = 'RMS_AllPeriods.{}'.format(fig_format)
            else:
                save_fn_basename = '{0:02}_RMS_{1:.5g}_s.{2}'.format(self.period_index,
                                                                 self.residual.period_list[self.period_index],
                                                                 fig_format)
        save_fn = os.path.join(self.save_path, save_fn_basename)

        if save_fig_dpi is not None:
            self.fig_dpi = save_fig_dpi

        self.fig.savefig(save_fn, dpi=self.fig_dpi)
        print('saved file to {0}'.format(save_fn))

        if fig_close:
            plt.close(self.fig)

    def plot_loop(self, fig_format='png'):
        """
        loop over all periods and save figures accordingly
        """
        self.read_residual_fn()

        for f_index in range(self.residual.period_list.size):
            self.period_index = f_index
            self.plot()
            self.save_figure(fig_format=fig_format)

# ==================================================================================
# FZ: add example usage code
# Justdo>  python mtpy/modeling/modem/plot_rms_maps.py
# ==================================================================================
if __name__ == "__main__":

    from mtpy.mtpy_globals import *

    # directory where files are located
    wd = os.path.join(SAMPLE_DIR, 'ModEM')

    # file stem for inversion result
    filestem = 'Modular_MPI_NLCG_004'

    # directory to save to
    save_path = NEW_TEMP_DIR


    # period index to plot (0 plots the first (shortest) period, 1 for the second, etc)
    period_index = 0

    # plot map
    rmsmap = PlotRMSMaps(residual_fn=os.path.join(wd, filestem + '.res'), period_index=period_index,
                         xminorticks=50000, yminorticks=50000, save_plots='y', plot_yn='n')
    rmsmap.plot()

    rmsmap.save_figure(save_path, fig_close=False)  # this will save a file to
