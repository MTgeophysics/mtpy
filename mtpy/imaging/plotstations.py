# -*- coding: utf-8 -*-
"""
===============
PlotStations
===============

Plots station locations in map view.


Created on Fri Jun 07 18:20:00 2013

@author: jpeacock-pr
"""

#==============================================================================

import matplotlib.pyplot as plt
import numpy as np
import os
import mtpy.imaging.mtplottools as mtpt
import mtpy.utils.exceptions as mtex

#==============================================================================


class PlotStations(object):
    """
    plot station locations in map view.

    Need to input one of the following lists:

    Arguments:
    ----------
        **fn_list** : list of strings
                     full paths to .edi files to plot. *default* is None

        **mt_object** : class mtpy.imaging.mtplot.MTplot
                        object of mtpy.imaging.mtplot.MTplot
                        *default* is None

    Optional Key Words:
    -------------------
        *fig_dpi*: float
                   dots per inch resolution of figure. *default* is 300.

        *fig_num*: int
                   number of figure instance. *default* is 1.

        *fig_size*: [x, y]
                    figure dimensions in inches. *default* is None.

        *font_size*: float
                     size of tick labels, axes labels will be +2.
                     *default* is 7

        *image_extent*: (xmin, xmax, ymin, ymax)
                        extent of image in map coordinates, must be input if
                        image_file is not None. *default* is None.

        *image_file*: string
                      full path to base image file, can be .jpg, .png, or .svg
                      *default* is None.

        *map_scale*: [ 'latlon' | 'eastnorth' | 'eastnorthkm' ]
                     scale of map, either in:
                         - 'latlon' --> latitude and longitude in decimal
                                        degrees
                         - 'eastnorth' --> easting and northing in meters
                         - 'eastnorthkm' --> easting and northing in kilometers

        *marker*: string or int
                  type of marker used to represent station location.  For all
                  marker options:

                ============================== ================================
                marker                         description
                ============================== ================================
                ``7``                          caretdown
                ``4``                          caretleft
                ``5``                          caretright
                ``6``                          caretup
                ``'o'``                        circle
                ``'D'``                        diamond
                ``'h'``                        hexagon1
                ``'H'``                        hexagon2
                ``'_'``                        hline
                ``''``                         nothing
                ``'None'``                     nothing
                ``None``                       nothing
                ``' '``                        nothing
                ``'8'``                        octagon
                ``'p'``                        pentagon
                ``','``                        pixel
                ``'+'``                        plus
                ``'.'``                        point
                ``'s'``                        square
                ``'*'``                        star
                ``'d'``                        thin_diamond
                ``3``                          tickdown
                ``0``                          tickleft
                ``1``                          tickright
                ``2``                          tickup
                ``'1'``                        tri_down
                ``'3'``                        tri_left
                ``'4'``                        tri_right
                ``'2'``                        tri_up
                ``'v'``                        triangle_down
                ``'<'``                        triangle_left
                ``'>'``                        triangle_right
                ``'^'``                        triangle_up
                ``'|'``                        vline
                ``'x'``                        x
                ``'$...$'``                    render the string using mathtext
                ============================== ================================

        *marker_color*: string or (red, green, blue) on a scale of 0 to 1.
                        color of station marker. *default* is black

        *marker_size*: float
                       size of station marker in points. *default* is 10.

        *plot_names*: [ True | False ]
                     plot station names next to marker. *default* is True.

        *plot_title*: string
                      title of plot

        *plot_yn*: [ 'y' | 'n' ]
                   plot on initialization.  *default* is 'y'.

        *ref_point*: (x, y)
                     reference point to center map on. *default* is (0, 0).

        *text_angle*: float
                      angle of station label in degrees. *default* is 0.

        *text_color*: string or (red, green, blue) on a scale [0,1]
                      color of station label. *default* is black.

        *text_ha*: [ 'center' | 'right' | 'left' ]
                   horizontal alignment of station label relative to the
                   station location. *default* is 'center'.

        *text_pad*: float
                    padding from station marker to station label in map
                    units.

        *text_size*: float
                     font size of station label.  *default* is 7.

        *text_va*: [ 'center' | 'top' | 'bottom' | 'baseline' ]
                   vertical alignment of station label. *default* is 'baseline'

        *text_weight*: [ 'ultralight' | 'light' | 'normal' | 'regular' |
                         'book' | 'medium' | 'roman' | 'semibold' |
                         'demibold' | 'demi' | 'bold' | 'heavy' |
                         'extra bold' | 'black' ]
                       weight of station label.  *default* is 'normal'.

        *xlimits*: (xmin, xmax)
                   limits of map in east-west direction.  *default* is None,
                   which computes limits from data.

        *ylimits*: (ymin, ymax)
                   limits of map in north-south direction.  *default* is None,
                   which computes limits from data.

    :Example: ::

        >>> import mtpy.imaging.mtplot as mtplot
        >>> import os
        >>> edipath = '/home/MT/edifiles'
        >>> edilist = [os.path.join(edipath, edi)
        >>> ...       for edi in os.listdir(edipath)
        >>> ...       if edi.find('.edi')>0]
        >>> ps1 = mtplot.plot_station_locations(fn_list=edilist)
        >>> # change station label padding and properties
        >>> ps1.text_pad = .001
        >>> ps1.text_angle = 60
        >>> ps1.text_size = 8
        >>> ps1.text_color = (.5, .5, 0) #orangeish
        >>> ps1.redraw_plot()
        >>> ps1.save_plot('/home/MT/figures', file_format='pdf')
        saved figure to '/home/MT/figures/station_map.pdf'

    =================== =======================================================
     Attributes          Description
    =================== =======================================================
        ax              matplotlib.axes instance of station map
        fig             matplotlib.figure instance of map figure
        fig_dpi         dots-per-inch resolution of figure
        fig_num         number of figure instance
        fig_size        size of figure in inches
        font_size       font size of tick labels
        image_extent    (xmin, xmax, ymin, ymax) extent of image if input
        image_file      full path to base image file, can be .jpg, .png or .svg
        map_scale       [ 'latlon' | 'eastnorth' | 'eastnorthkm' ] map scale
        marker          station marker, see above for options
        marker_color    color of marker
        marker_size     size of marker in points
        mt_list          list of mtpy.imaging.mtplottools.MTplot instances
        plot_names      [ True | False ] plot station names next to markers
        plot_title      title of plot
        plot_yn         [ 'y' | 'n' ] plot on initializing PlotStations
        ref_point       reference point to center map on.
        stationid       (index0, index1) to get station label from station name
        text_angle      angle of station label
        text_color      color of station label
        text_ha         horizontal alignment of station label, see above
        text_pad        padding of station label in y direction
        text_size       font size of station label
        text_va         vertical alignment of station label
        text_weight     font weight of station label
        xlimits         limits of map in east-west direction
        ylimits         limits of map in north-south direction
    =================== =======================================================

    Methods:
    ---------
        * *plot*: plots the pseudosection according to keywords
        * *redraw_plot*: redraws the plot, use if you change some of the
                         attributes.
        * *update_plot*: updates the plot, use if you change some of the
                         axes attributes, figure needs to be open to update.
        * *save_plot*: saves the plot to given filepath.

    """

    def __init__(self, **kwargs):

        fn_list = kwargs.pop('fn_list', None)
        mt_object_list = kwargs.pop('mt_object_list', None)

        #----set attributes for the class-------------------------
        self.mt_list = mtpt.MTplot_list(fn_list=fn_list,
                                        mt_object_list=mt_object_list)

        #--> set plot properties
        self.fig_num = kwargs.pop('fig_num', 1)
        self.plot_title = kwargs.pop('plot_title', None)
        self.fig_dpi = kwargs.pop('fig_dpi', 300)
        self.fig_size = kwargs.pop('fig_size', None)
        self.font_size = kwargs.pop('font_size', 7)
        self.stationid = kwargs.pop('stationid', [0, 4])
        self.xlimits = kwargs.pop('xlimits', None)
        self.ylimits = kwargs.pop('ylimits', None)
        self.ref_point = kwargs.pop('ref_point', (0, 0))

        self.map_scale = kwargs.pop('map_scale', 'latlon')
        self.marker = kwargs.pop('marker', 'v')
        self.marker_size = kwargs.pop('marker_size', 10)
        self.marker_color = kwargs.pop('marker_color', 'k')
        self.plot_names = kwargs.pop('plot_names', True)

        self.text_size = kwargs.pop('text_size', 7)
        self.text_weight = kwargs.pop('text_weight', 'normal')
        self.text_color = kwargs.pop('text_color', 'k')
        self.text_ha = kwargs.pop('text_ha', 'center')
        self.text_va = kwargs.pop('text_va', 'baseline')
        self.text_angle = kwargs.pop('text_angle', 0)
        self.text_x_pad = kwargs.pop('text_x_pad', None)
        self.text_y_pad = kwargs.pop('text_y_pad', None)

        self.image_file = kwargs.pop('image_file', None)
        self.image_extent = kwargs.pop('image_extent', None)

        if self.image_file is not None:
            if self.image_extent is None:
                raise mtex.MTpyError_inputarguments('Need to input extents ' +
                                                    'of the image as' +
                                                    '(x0, y0, x1, y1)')

        #--> plot if desired
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()

    def plot(self):
        """
        plots the station locations

        """
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = .09
        plt.rcParams['figure.subplot.right'] = .98
        plt.rcParams['figure.subplot.bottom'] = .09
        plt.rcParams['figure.subplot.top'] = .98

        # get station locations
        self.mt_list.get_station_locations(map_scale=self.map_scale,
                                           ref_point=self.ref_point)

        text_dict = {'size': self.text_size,
                     'weight': self.text_weight,
                     'rotation': self.text_angle,
                     'color': self.text_color}

        font_dict = {'size': self.font_size + 2, 'weight': 'bold'}

        if self.xlimits is None:
            if np.sign(self.mt_list.map_xarr.min()) == -1:
                self.xlimits = (self.mt_list.map_xarr.min() * 1.002,
                                self.mt_list.map_xarr.max() * .998)
            else:
                self.xlimits = (self.mt_list.map_xarr.min() * .998,
                                self.mt_list.map_xarr.max() * 1.002)

        if self.ylimits is None:
            if np.sign(self.mt_list.map_yarr.min()) == -1:
                self.ylimits = (self.mt_list.map_yarr.min() * 1.002,
                                self.mt_list.map_yarr.max() * .998)
            else:
                self.ylimits = (self.mt_list.map_yarr.min() * .998,
                                self.mt_list.map_yarr.max() * 1.002)

        if self.map_scale == 'latlon':
            xlabel = 'Longitude (deg)'
            ylabel = 'Latitude (deg)'

        elif self.map_scale == 'eastnorth':
            xlabel = 'Easting (m)'
            ylabel = 'Northing (m)'

        elif self.map_scale == 'eastnorthkm':
            xlabel = 'Easting (km)'
            ylabel = 'Northing (km)'

        # make a figure instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)

        #add and axes
        self.ax = self.fig.add_subplot(1, 1, 1, aspect='equal')

        #--> plot the background image if desired-----------------------
        if self.image_file is not None:
            im = plt.imread(self.image_file)
            self.ax.imshow(im, origin='lower', extent=self.image_extent,
                           aspect='auto')

        for key in list(self.mt_list.map_dict.keys()):
            self.ax.scatter(self.mt_list.map_dict[key][0],
                            self.mt_list.map_dict[key][1],
                            marker=self.marker,
                            c=self.marker_color,
                            s=self.marker_size)

            if self.plot_names == True:
                if self.text_x_pad is None:
                    self.text_x_pad = .0009 * self.mt_list.map_dict[key][0]
                if self.text_y_pad is None:
                    self.text_y_pad = .0009 * self.mt_list.map_dict[key][1]

                self.ax.text(self.mt_list.map_dict[key][0] + self.text_x_pad,
                             self.mt_list.map_dict[key][1] + self.text_y_pad *
                             np.sign(self.mt_list.map_dict[key][1]),
                             key[self.stationid[0]:self.stationid[1]],
                             verticalalignment=self.text_va,
                             horizontalalignment=self.text_ha,
                             fontdict=text_dict)

        # set axis properties
        self.ax.set_xlabel(xlabel, fontdict=font_dict)
        self.ax.set_ylabel(ylabel, fontdict=font_dict)
        self.ax.grid(alpha=.35, color=(.25, .25, .25))
        self.ax.set_xlim(self.xlimits)
        self.ax.set_ylim(self.ylimits)

        plt.show()

    def write_station_locations(self, save_path=None):
        """
        Write text file containing station locations in map coordinates and
        relative to ref_point.

        Arguments:
        ----------
            **save_path**: string
                           full path to folder to save file, or full path to
                           the file to save to. *default* is None, which uses
                           the directory path of files used to plot.

        Returns:
        ---------
            **fn_save_path**: string
                              full path to text file

        """

        if save_path is None:
            try:
                svpath = os.path.dirname(self.mt_list.mt_list[0].fn)
            except TypeError:
                raise IOError('Need to input save_path, could not find a path')
        else:
            svpath = save_path

        if self.map_scale == 'latlon':
            hdr_list = ['Station', 'Longitude(deg)', 'Latitude(deg)',
                        'Elevation(m)']
        elif self.map_scale == 'eastnorth':
            hdr_list = ['Station', 'Easting(m)', 'Northing(m)',
                        'Elevation(m)']
        elif self.map_scale == 'eastnorthkm':
            hdr_list = ['Station', 'Easting(km)', 'Northing(km)',
                        'Elevation(m)']

        self.mt_list.get_station_locations(map_scale=self.map_scale,
                                           ref_point=self.ref_point)

        fn_svpath = os.path.join(svpath, 'StationLocations_{0}.txt'.format(
            self.map_scale))
        tfid = file(fn_svpath, 'w')

        hdr_str = ['{0:<15}'.format(hdr_list[0])] +\
                  ['{0:^15}'.format(hh) for hh in hdr_list[1:]] + ['\n']

        tfid.write(''.join(hdr_str))
        for ss in list(self.mt_list.map_dict.keys()):
            x = self.mt_list.map_dict[ss][0]
            y = self.mt_list.map_dict[ss][1]
            z = self.mt_list.map_dict[ss][2]
            if self.map_scale == 'latlon':
                tline = '{0:<15}{1: ^15.3f}{2: ^15.3f}{3: ^15.1f}\n'.format(ss,
                                                                            x,
                                                                            y,
                                                                            z)
            else:
                tline = '{0:<15}{1: ^15.1f}{2: ^15.1f}{3: ^15.1f}\n'.format(ss,
                                                                            x,
                                                                            y,
                                                                            z)
            tfid.write(tline)

        tfid.close()

        print('Saved file to: ', fn_svpath)

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

        sf = '_{0:.6g}'.format(self.plot_freq)

        if fig_dpi is None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation)
            plt.clf()
            plt.close(self.fig)

        else:
            if not os.path.exists(save_fn):
                os.mkdir(save_fn)
            if not os.path.exists(os.path.join(save_fn, 'station_map')):
                os.mkdir(os.path.join(save_fn, 'station_map'))
                save_fn = os.path.join(save_fn, 'station_map')

            save_fn = os.path.join(save_fn, 'PTmap_' + self.ellipse_colorby + sf +
                                   'Hz.' + file_format)
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
        return "Plot station locations"
