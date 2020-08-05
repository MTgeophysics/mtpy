"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch

Revision History:
    brenainn.moushall@ga.gov.au 31-03-2020 13:38:10 AEDT:
        - Add ability to plot on background geotiff
        - Add ability to write RMS map as shapefile
        - Add 'plot_elements' attribute for selecting whether to plot
          impedance, tippers or both
        - Allow selection of period by providing period in seconds
"""
import os
import logging

import numpy as np
import geopandas as gpd
from shapely.geometry import Point, Polygon
from matplotlib import colors as colors, pyplot as plt, colorbar as mcb, cm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from mtpy.utils import basemap_tools
from mtpy.utils.plot_geotiff_imshow import plot_geotiff_on_axes
from mtpy.utils.mtpylog import MtPyLog
from mtpy.utils.gis_tools import epsg_project
from mtpy.modeling.modem import Data, Residual

__all__ = ['PlotRMSMaps']
_logger = MtPyLog.get_mtpy_logger(__name__)

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
    bimg                path to a geotiff to display as background of
                        plotted maps
    bimg_band           band of bimg to plot. *default* is None, which 
                        will plot all available bands
    bimg_cmap           cmap for bimg. *default* is 'viridis'. Ignored 
                        if bimg is RBG/A
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
        self._residual_fn = None
        self.residual = None
        self.residual_fn = residual_fn
        self.model_epsg = kwargs.pop('model_epsg', None)
        self.read_residual_fn()

        self.save_path = kwargs.pop('save_path', os.path.dirname(self.residual_fn))

        self.period = kwargs.pop('period', None)
        if self.period is not None:
            # Get period index closest to provided period
            index = np.argmin(np.fabs(self.residual.period_list - self.period))
            _logger.info("Plotting nearest available period ({}s) for selected period ({}s)"
                         .format(self.residual.period_list[index], self.period))
            self.period_index = index
        else:
            self.period_index = kwargs.pop('period_index', 0)

        self.plot_elements = kwargs.pop('plot_elements', 'both')

        self.subplot_left = kwargs.pop('subplot_left', .1)
        self.subplot_right = kwargs.pop('subplot_right', .9)
        self.subplot_top = kwargs.pop('subplot_top', .95)
        self.subplot_bottom = kwargs.pop('subplot_bottom', .1)
        self.subplot_hspace = kwargs.pop('subplot_hspace', .1)
        self.subplot_vspace = kwargs.pop('subplot_vspace', .01)

        self.font_size = kwargs.pop('font_size', 8)

        self.fig = None
        self.fig_size = kwargs.pop('fig_size', [7.75, 6.75])
        self.fig_dpi = kwargs.pop('fig_dpi', 200)
        self.fig_num = kwargs.pop('fig_num', 1)
        self.font_dict = {'size': self.font_size + 2, 'weight': 'bold'}

        self.marker = kwargs.pop('marker', 's')
        self.marker_size = kwargs.pop('marker_size', 10)

        self.rms_max = kwargs.pop('rms_max', 5)
        self.rms_min = kwargs.pop('rms_min', 0)

        self.tick_locator = kwargs.pop('tick_locator', None)
        self.pad_x = kwargs.pop('pad_x', None)
        self.pad_y = kwargs.pop('pad_y', None)

        self.plot_yn = kwargs.pop('plot_yn', 'y')

        self.bimg = kwargs.pop('bimg', None)
        if self.bimg and self.model_epsg is None:
            _logger.warning("You have provided a geotiff as a background image but model_epsg is "
                            "not set. It's assumed that the CRS of the model and the CRS of the "
                            "geotiff are the same. If this is not the case, please provide "
                            "model_epsg to PlotRMSMaps.")
        self.bimg_band = kwargs.pop('bimg_band', None)
        self.bimg_cmap = kwargs.pop('bimg_cmap', 'viridis')

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

        self.rms_cmap = None
        if 'rms_cmap' in list(kwargs.keys()):
            # check if it is a valid matplotlib color stretch
            if kwargs['rms_cmap'] in dir(cm):
                self.rms_cmap = cm.get_cmap(kwargs['rms_cmap'])
            else:
                print("provided rms_cmap invalid, using default colormap")

        if self.rms_cmap is None:
            self.rms_cmap = colors.LinearSegmentedColormap('rms_cmap',
                                                           self.rms_cmap_dict,
                                                           256)

        if self.plot_elements == 'both':
            self.plot_z_list = [{'label': r'$Z_{xx}$', 'index': (0, 0), 'plot_num': 1},
                                {'label': r'$Z_{xy}$', 'index': (0, 1), 'plot_num': 2},
                                {'label': r'$Z_{yx}$', 'index': (1, 0), 'plot_num': 3},
                                {'label': r'$Z_{yy}$', 'index': (1, 1), 'plot_num': 4},
                                {'label': r'$T_{x}$', 'index': (0, 0), 'plot_num': 5},
                                {'label': r'$T_{y}$', 'index': (0, 1), 'plot_num': 6}]
        elif self.plot_elements == 'impedance':
            self.plot_z_list = [{'label': r'$Z_{xx}$', 'index': (0, 0), 'plot_num': 1},
                                {'label': r'$Z_{xy}$', 'index': (0, 1), 'plot_num': 2},
                                {'label': r'$Z_{yx}$', 'index': (1, 0), 'plot_num': 3},
                                {'label': r'$Z_{yy}$', 'index': (1, 1), 'plot_num': 4}]
        elif self.plot_elements == 'tippers':
            self.plot_z_list = [{'label': r'$T_{x}$', 'index': (0, 0), 'plot_num': 1},
                                {'label': r'$T_{y}$', 'index': (0, 1), 'plot_num': 2}]
        else:
            raise ValueError("'plot_elements' value '{}' is not recognised. Please set "
                             "'plot_elements' to 'impedance', 'tippers' or 'both'.")

        if self.plot_yn == 'y':
            self.plot()

    def _fig_title(self, font_size, font_weight):
        if self.period_index == 'all':
            title = 'All periods'
        else:
            title = 'period = {0:.5g} (s)'.format(self.residual.period_list[self.period_index])
        self.fig.suptitle(title, fontdict={'size': font_size, 'weight': font_weight})

    def _calculate_rms(self, plot_dict):
        ii = plot_dict['index'][0]
        jj = plot_dict['index'][1]

        rms = np.zeros(self.residual.residual_array.shape[0])
        self.residual.get_rms()
        if plot_dict['label'].startswith('$Z'):
            rms = self.residual.rms_array['rms_z_component_period'][:, self.period_index, ii, jj]
        elif plot_dict['label'].startswith('$T'):
            rms = self.residual.rms_array['rms_tip_component_period'][:, self.period_index, ii, jj]
        
        # for ridx in range(len(self.residual.residual_array)):

        #     if self.period_index == 'all':
        #         r_arr = self.residual.rms_array[ridx]
        #         if plot_dict['label'].startswith('$Z'):
        #             rms[ridx] = r_arr['rms_z']
        #         else:
        #             rms[ridx] = r_arr['rms_tip']
        #     else:
        #         r_arr = self.residual.residual_array[ridx]
        #         # calulate the rms self.residual/error
        #         if plot_dict['label'].startswith('$Z'):
        #             rms[ridx] = r_arr['z'][self.period_index, ii, jj].__abs__() / \
        #                 r_arr['z_err'][self.period_index, ii, jj].real

        #         else:
        #             rms[ridx] = r_arr['tip'][self.period_index, ii, jj].__abs__() / \
        #                 r_arr['tip_err'][self.period_index, ii, jj].real

        filt = np.nan_to_num(rms).astype(bool)

        if len(rms[filt]) == 0:
            _logger.warning("No RMS available for component {}"
                            .format(self._normalize_label(plot_dict['label'])))

        return rms, filt

    @staticmethod
    def _normalize_label(label):
        return label.replace('$', '').replace('{', '').replace('}', '').replace('_', '')

    def read_residual_fn(self):
        if self.residual is None:
            self.residual = Residual(residual_fn=self.residual_fn,
                                     model_epsg=self.model_epsg)
            self.residual.read_residual_file()
            self.residual.get_rms()
        else:
            pass

    def create_shapefiles(self, dst_epsg, save_path=None):
        """
        Creates RMS map elements as shapefiles which can displayed in a
        GIS viewer. Intended to be called as part of the 'plot' 
        function.

        The points to plot defined by `lons` and `lats` are the centre 
        of the rectangular markers.

        If `model_epsg` hasn't been set on class, then 4326 is assumed.

        Parameters
        ----------
        dst_epsg : int
            EPSG code of the CRS that Shapefiles will be projected to.
            Make this the same as the CRS of the geotiff you intend to
            display on.
        marker_width : float
             Radius of the circular markers. Units are defined by
             `model_epsg`.
        """
        if save_path is None:
            save_path = self.save_path
        lon = self.residual.residual_array['lon']
        lat = self.residual.residual_array['lat']
        if self.model_epsg is None:
            _logger.warning("model_epsg has not been provided. Model EPSG is assumed to be 4326. "
                            "If this is not correct, please provide model_epsg to PlotRMSMaps. "
                            "Otherwise, shapefiles may have projection errors.")
            src_epsg = 4326
        else:
            src_epsg = self.model_epsg
        src_epsg = {'init': 'epsg:{}'.format(src_epsg)}
        for p_dict in self.plot_z_list:
            rms, _ = self._calculate_rms(p_dict)
            markers = []
            for x, y in zip(lon, lat):
                markers.append(Point(x, y))

            df = gpd.GeoDataFrame({'lon': lon, 'lat': lat, 'rms': rms},
                                  crs=src_epsg, geometry=markers)
            df.to_crs(epsg=dst_epsg, inplace=True)

            if self.period_index == 'all':
                period = 'all'
            else:
                period = self.residual.period_list[self.period_index]
            filename = '{}_EPSG_{}_Period_{}.shp'.format(self._normalize_label(p_dict['label']),
                                                         dst_epsg, period)
            directory = os.path.join(self.save_path, 'shapefiles_for_period_{}s'.format(period))
            if not os.path.exists(directory):
                os.mkdir(directory)
            outpath = os.path.join(directory, filename)
            df.to_file(outpath, driver='ESRI Shapefile')
            print("Saved shapefiles to %s", outpath)

    def plot(self):
        """
        plot rms in map view
        """
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

        # Get number of rows based on what is being plotted.
        sp_rows, sp_cols = len(self.plot_z_list) / 2, 2
        # Adjust dimensions based on number of rows.
        # Hardcoded - having issues getting the right spacing between
        # labels and subplots.
        if sp_rows == 1:
            self.fig_size[1] = 3.
        elif sp_rows == 2:
            self.fig_size[1] = 5.6

        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        plt.rcParams['figure.subplot.wspace'] = self.subplot_hspace
        plt.rcParams['figure.subplot.hspace'] = self.subplot_vspace
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)

        lon = self.residual.residual_array['lon']
        lat = self.residual.residual_array['lat']

        for p_dict in self.plot_z_list:
            rms, filt = self._calculate_rms(p_dict)

            ax = self.fig.add_subplot(sp_rows, sp_cols, 
                                      p_dict['plot_num'],
                                      aspect='equal')

            plt.scatter(lon[filt],
                        lat[filt],
                        c=rms[filt],
                        marker=self.marker,
                        edgecolors=(0, 0, 0),
                        cmap=self.rms_cmap,
                        norm=colors.Normalize(vmin=self.rms_min,
                                              vmax=self.rms_max),
                        )

            if not np.all(filt):
                filt2 = (1 - filt).astype(bool)
                plt.plot(lon[filt2],
                         lat[filt2],
                         '.',
                         ms=0.1,
                         mec=(0, 0, 0),
                         mfc=(1, 1, 1)
                         )

            # Hide y-ticks on subplots in column 2.
            if p_dict['plot_num'] in (2, 4, 6):
                plt.setp(ax.get_yticklabels(), visible=False)
            else:
                ax.set_ylabel('Latitude (deg)', fontdict=self.font_dict)

            # Only show x-ticks in final row.
            if p_dict['plot_num'] in (sp_rows * 2 - 1, sp_rows * 2):
                ax.set_xlabel('Longitude (deg)', fontdict=self.font_dict)
            else:
                plt.setp(ax.get_xticklabels(), visible=False)

            ax.text(self.residual.residual_array['lon'].min() + .005 - self.pad_x,
                    self.residual.residual_array['lat'].max() - .005 + self.pad_y,
                    p_dict['label'],
                    verticalalignment='top',
                    horizontalalignment='left',
                    bbox={'facecolor': 'white'},
                    zorder=3)

            ax.tick_params(direction='out')
            ax.grid(zorder=0, color=(.75, .75, .75))

            ax.set_xlim(self.residual.residual_array['lon'].min() - self.pad_x,
                        self.residual.residual_array['lon'].max() + self.pad_x)

            ax.set_ylim(self.residual.residual_array['lat'].min() - self.pad_y,
                        self.residual.residual_array['lat'].max() + self.pad_y)

            if self.bimg:
                plot_geotiff_on_axes(self.bimg, ax, epsg_code=self.model_epsg,
                                     band_number=self.bimg_band, cmap=self.bimg_cmap)

            ax.xaxis.set_major_locator(MultipleLocator(self.tick_locator))
            ax.yaxis.set_major_locator(MultipleLocator(self.tick_locator))
            ax.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%2.2f'))

        cb_ax = self.fig.add_axes([self.subplot_right + .02, .225, .02, .45])
        color_bar = mcb.ColorbarBase(cb_ax,
                                     cmap=self.rms_cmap,
                                     norm=colors.Normalize(vmin=self.rms_min,
                                                           vmax=self.rms_max),
                                     orientation='vertical')

        color_bar.set_label('RMS', fontdict=self.font_dict)

        self._fig_title(font_size=self.font_size + 3, font_weight='bold')
        self.fig.show()

    # BM: Is this still in use? `Residual` has no attribute `data_array`
    # which breaks this function.
    def plot_map(self):
        """
        plot the misfit as a map instead of points
        """
        rms_1 = 1. / self.rms_max

        if self.tick_locator is None:
            x_locator = np.round((self.residual.data_array['lon'].max() -
                                  self.residual.data_array['lon'].min()) / 5, 2)
            y_locator = np.round((self.residual.data_array['lat'].max() -
                                  self.residual.data_array['lat'].min()) / 5, 2)

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

        lat_arr = self.residual.data_array['lat']
        lon_arr = self.residual.data_array['lon']
        data_points = np.array([lon_arr, lat_arr])

        interp_lat = np.linspace(lat_arr.min(),
                                 lat_arr.max(),
                                 3 * self.residual.data_array.size)

        interp_lon = np.linspace(lon_arr.min(),
                                 lon_arr.max(),
                                 3 * self.residual.data_array.size)

        grid_x, grid_y = np.meshgrid(interp_lon, interp_lat)

        # calculate rms
        z_err = self.residual.data_array['z_err'].copy()
        z_err[np.where(z_err == 0.0)] = 1.0
        z_rms = np.abs(self.residual.data_array['z']) / z_err.real

        t_err = self.residual.data_array['tip_err'].copy()
        t_err[np.where(t_err == 0.0)] = 1.0
        t_rms = np.abs(self.residual.data_array['tip']) / t_err.real

        # --> plot maps
        for p_dict in self.plot_z_list:
            ax = self.fig.add_subplot(3, 2, p_dict['plot_num'], aspect='equal')

            if p_dict['plot_num'] == 1 or p_dict['plot_num'] == 3:
                ax.set_ylabel('Latitude (deg)', fontdict=self.font_dict)
                plt.setp(ax.get_xticklabels(), visible=False)

            elif p_dict['plot_num'] == 2 or p_dict['plot_num'] == 4:
                plt.setp(ax.get_xticklabels(), visible=False)
                plt.setp(ax.get_yticklabels(), visible=False)

            elif p_dict['plot_num'] == 6:
                plt.setp(ax.get_yticklabels(), visible=False)
                ax.set_xlabel('Longitude (deg)', fontdict=self.font_dict)

            else:
                ax.set_xlabel('Longitude (deg)', fontdict=self.font_dict)
                ax.set_ylabel('Latitude (deg)', fontdict=self.font_dict)

            ax.text(self.residual.data_array['lon'].min() + .005 - self.pad_x,
                    self.residual.data_array['lat'].max() - .005 + self.pad_y,
                    p_dict['label'],
                    verticalalignment='top',
                    horizontalalignment='left',
                    bbox={'facecolor': 'white'},
                    zorder=3)

            ax.tick_params(direction='out')
            ax.grid(zorder=0, color=(.75, .75, .75), lw=.75)

            # [line.set_zorder(3) for line in ax.lines]

            ax.set_xlim(self.residual.data_array['lon'].min() - self.pad_x,
                        self.residual.data_array['lon'].max() + self.pad_x)

            ax.set_ylim(self.residual.data_array['lat'].min() - self.pad_y,
                        self.residual.data_array['lat'].max() + self.pad_y)

            ax.xaxis.set_major_locator(MultipleLocator(self.tick_locator))
            ax.yaxis.set_major_locator(MultipleLocator(self.tick_locator))
            ax.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%2.2f'))

            # -----------------------------
            ii = p_dict['index'][0]
            jj = p_dict['index'][1]

            # calulate the rms self.residual/error
            if p_dict['plot_num'] < 5:
                rms = z_rms[:, self.period_index, ii, jj]
            else:
                rms = t_rms[:, self.period_index, ii, jj]

            # check for non zeros
            nz = np.nonzero(rms)
            data_points = np.array([lon_arr[nz], lat_arr[nz]])
            rms = rms[nz]
            if len(rms) < 5:
                continue

            # interpolate onto a grid
            rms_map = interpolate.griddata(data_points.T,
                                           rms,
                                           (grid_x, grid_y),
                                           method='cubic')
            # plot the grid
            im = ax.pcolormesh(grid_x,
                               grid_y,
                               rms_map,
                               cmap=self.rms_map_cmap,
                               vmin=self.rms_min,
                               vmax=self.rms_max,
                               zorder=3)
            ax.grid(zorder=0, color=(.75, .75, .75), lw=.75)

        # cb_ax = mcb.make_axes(ax, orientation='vertical', fraction=.1)
        cb_ax = self.fig.add_axes([self.subplot_right + .02, .225, .02, .45])
        color_bar = mcb.ColorbarBase(cb_ax,
                                     cmap=self.rms_map_cmap,
                                     norm=colors.Normalize(vmin=self.rms_min,
                                                           vmax=self.rms_max),
                                     orientation='vertical')

        color_bar.set_label('RMS', fontdict=self.font_dict)

        self._fig_title(font_size=self.font_size + 3, font_weight='bold')
        self.fig.show()

    def basemap_plot(self, datatype='all', tick_interval=None, save=False,
                     savepath=None, new_figure=True, mesh_rotation_angle=0.,
                     show_topography=False, **basemap_kwargs):
        """
        plot RMS misfit on a basemap using basemap modules in matplotlib

        :param datatype: type of data to plot misfit for, either 'z', 'tip', or
                         'all' to plot overall RMS
        :param tick_interval: tick interval on map in degrees, if None it is 
                              calculated from the data extent
        :param save: True/False, whether or not to save and close figure
        :param savepath: full path of file to save to, if None, saves to 
                         self.save_path
        :new_figure: True/False, whether or not to initiate a new figure for
                     the plot
        :param mesh_rotation_angle: rotation angle of mesh, in degrees 
                                    clockwise from north
        :param show_topography: True/False, option to show the topograpy in the
                                background
        :param **basemap_kwargs: provide any valid arguments to Basemap 
                                 instance (e.g. projection etc - see 
                                 https://basemaptutorial.readthedocs.io/en/latest/basemap.html)
                                 and these will be passed to the map.

        """
        if self.model_epsg is None:
            print("No projection information provided, please provide the model epsg code relevant to your model")
            return

        if new_figure:
            self.fig = plt.figure()

        # rotate stations
        if mesh_rotation_angle != 0:
            if hasattr(self, 'mesh_rotation_angle'):
                angle_to_rotate = self.mesh_rotation_angle - mesh_rotation_angle
            else:
                angle_to_rotate = -mesh_rotation_angle

            self.mesh_rotation_angle = mesh_rotation_angle

            self.residual.station_locations.rotate_stations(angle_to_rotate)

        # get relative locations
        seast, snorth = self.residual.station_locations.rel_east + self.residual.station_locations.center_point['east'],\
            self.residual.station_locations.rel_north + self.residual.station_locations.center_point['north']

        # project station location eastings and northings to lat/long
        slon, slat = epsg_project(seast, snorth, self.model_epsg, 4326)
        self.residual.station_locations.station_locations['lon'] = slon
        self.residual.station_locations.station_locations['lat'] = slat

        # initialise a basemap with extents, projection etc calculated from data
        # if not provided in basemap_kwargs # BM: todo? 
        self.bm = basemap_tools.initialise_basemap(self.residual.station_locations, **basemap_kwargs)
        basemap_tools.add_basemap_frame(self.bm, tick_interval=tick_interval)

        # project to basemap coordinates
        sx, sy = self.bm(slon, slat)

        # make scatter plot
        if datatype == 'all':
            if self.period_index == 'all':
                rms = self.residual.rms_array['rms']
            else:
                rms = self.residual.rms_array['rms_period'][:, self.period_index]
        elif datatype in ['z', 'tip']:
            if self.period_index == 'all':
                rms = self.residual.rms_array['rms_{}'.format(datatype)]
            else:
                rms = self.residual.rms_array['rms_{}_period'.format(datatype)][:, self.period_index]

        filt = np.nan_to_num(rms).astype(bool)

        self.bm.scatter(sx[filt], sy[filt],
                        c=rms[filt],
                        marker=self.marker,
                        edgecolors=(0, 0, 0),
                        cmap=self.rms_cmap,
                        norm=colors.Normalize(vmin=self.rms_min,
                                              vmax=self.rms_max)
                        )

        if not np.all(filt):
            filt2 = (1 - filt).astype(bool)
            self.bm.plot(sx[filt2], sy[filt2], 'k.')

        color_bar = plt.colorbar(cmap=self.rms_cmap,
                                 shrink=0.6,
                                 norm=colors.Normalize(vmin=self.rms_min,
                                                       vmax=self.rms_max),
                                 orientation='vertical')

        color_bar.set_label('RMS')

        title_dict = {'all': 'Z + Tipper', 'z': 'Z', 'tip': 'Tipper'}

        if self.period_index == 'all':
            plt.title('RMS misfit over all periods for ' + title_dict[datatype])
        else:
            plt.title('RMS misfit for period = {0:.5g} (s)'.format(self.residual.period_list[self.period_index]))

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

    def plot_loop(self, fig_format='png', style='point'):
        """
        loop over all periods and save figures accordingly

        :param: style [ 'point' | 'map' ]
        """

        for f_index in range(self.residual.period_list.size):
            self.period_index = f_index
            if style == 'point':
                self.plot()
                self.save_figure(fig_format=fig_format)
            elif style == 'map':
                self.plot_map()
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
