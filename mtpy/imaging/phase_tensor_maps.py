# -*- coding: utf-8 -*-
"""
Plot phase tensor map in Lat-Lon Coordinate System

Revision History:
    Created by @author: jpeacock-pr on Thu May 30 18:20:04 2013

    Modified by Fei.Zhang@ga.gov.au 2017-03:

    brenainn.moushall 26-03-2020 15:07:14 AEDT:
        Add plotting of geotiff as basemap background.
"""
import os
import glob

import numpy as np
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.colorbar as mcb
import matplotlib.tri as tri

import mtpy.utils.gis_tools as gis_tools
import mtpy.imaging.mtcolors as mtcl
import mtpy.imaging.mtplottools as mtpl
import mtpy.analysis.pt as MTpt
from mtpy.utils.mtpylog import MtPyLog
from mtpy.utils.plot_geotiff_imshow import plot_geotiff_on_axes

_logger = MtPyLog.get_mtpy_logger(__name__)


# ==============================================================================

class PlotPhaseTensorMaps(mtpl.PlotSettings):
    """
    Plots phase tensor ellipses in map view from a list of edi files

    Arguments:
    -------------

        **fn_list** : list of strings
                          full paths to .edi files to plot

        **z_object** : class mtpy.core.z.Z
                      object of mtpy.core.z.  If this is input be sure the
                      attribute z.freq is filled.  *default* is None

        **mt_object** : class mtpy.imaging.mtplot.MTplot
                        object of mtpy.imaging.mtplot.MTplot
                        *default* is None

        **pt_object** : class mtpy.analysis.pt
                        phase tensor object of mtpy.analysis.pt.  If this is
                        input then the ._mt attribute is set to None cause
                        at the moment cannot tranform the phase tensor to z
                        *default* is None

        **plot_freq** : float
                             freq to plot in Hz
                             *default* is 1

        **ftol** : float
                   tolerance in freq range to look for in each file.
                   *default* is 0.1 (10 percent). ftol has no effect when
                   interpolate is True

        **interpolate** : bool
                   flag to indicate whether to interpolate values on to plot_freq.
                   *default* in True

        **ellipse_dict** : dictionary
                          dictionary of parameters for the phase tensor
                          ellipses with keys:
                              * 'size' -> size of ellipse in points
                                         *default* is 2

                              * 'colorby' : [ 'phimin' | 'phimax' | 'skew' |
                                              'skew_seg' | 'phidet' |
                                              'ellipticity' ]

                                        - 'phimin' -> colors by minimum phase
                                        - 'phimax' -> colors by maximum phase
                                        - 'skew' -> colors by skew
                                        - 'skew_seg' -> colors by skew in
                                                       discrete segments
                                                       defined by the range
                                        - 'normalized_skew' -> colors by skew
                                                see [Booker, 2014]
                                        - 'normalized_skew_seg' -> colors by
                                                       normalized skew in
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


        **cb_dict** : dictionary to control the color bar

                      * 'orientation' : [ 'vertical' | 'horizontal' ]
                                       orientation of the color bar
                                       *default* is vertical

                      * 'position' : tuple (x,y,dx,dy)
                                    - x -> lateral position of left hand corner
                                          of the color bar in figure between
                                          [0,1], 0 is left side

                                    - y -> vertical position of the bottom of
                                          the color bar in figure between
                                          [0,1], 0 is bottom side.

                                    - dx -> width of the color bar [0,1]

                                    - dy -> height of the color bar [0,1]

        **arrow_dict** : dictionary for arrow properties
                        * 'size' : float
                                  multiplier to scale the arrow. *default* is 5
                        * 'head_length' : float
                                         length of the arrow head *default* is
                                         1.5
                        * 'head_width' : float
                                        width of the arrow head *default* is
                                        1.5
                        * 'lw' : float
                                line width of the arrow *default* is .5

                        * 'color' : tuple (real, imaginary)
                                   color of the arrows for real and imaginary

                        * 'threshold': float
                                      threshold of which any arrow larger than
                                      this number will not be plotted, helps
                                      clean up if the data is not good.
                                      *default* is 1, note this is before
                                      scaling by 'size'

                        * 'direction : [ 0 | 1 ]
                                     - 0 for arrows to point toward a conductor
                                     - 1 for arrow to point away from conductor

        **xpad** : float
                   padding in the east-west direction of plot boundaries.  Note
                   this is by default set to lat and long units, so if you use
                   easting/northing put it the respective units.
                   *default* is 0.2

        **ypad** : float
                   padding in the north-south direction of plot boundaries.
                   Note this is by default set to lat and long units, so if you
                   use easting/northing put it the respective units.
                   *default* is 0.2

        **rotz** : float
                   angle in degrees to rotate the data clockwise positive.
                   *Default* is 0

        **fig_size** : tuple or list (x, y) in inches
                      dimensions of the figure box in inches, this is a default
                      unit of matplotlib.  You can use this so make the plot
                      fit the figure box to minimize spaces from the plot axes
                      to the figure box.  *default* is [8, 8]

        **station_dict** : dictionary
                            * 'id' --> for station id index.  Ex: If you want
                                      'S01' from 'S01dr' input as (0,3).

                            * 'pad' --> pad from the center of the ellipse to
                                       the station label

                            * 'font_dict'--> dictionary of font properties
                                           font dictionary for station name.
                                           Keys can be matplotlib.text
                                           properties, common ones are:
                                           * 'size'   -> for font size
                                           * 'weight' -> for font weight
                                           * 'color'  -> for color of font
                                           * 'angle'  -> for angle of text

        **tscale** : [ 'period' | 'freq' ]

                     * 'period'    -> plot vertical scale in period

                     * 'freq' -> plot vertical scale in freq

        **mapscale** : [ 'deg' | 'm' | 'km' ]
                       Scale of the map coordinates.

                       * 'deg' --> degrees in latitude and longitude

                       * 'm' --> meters for easting and northing

                       * 'km' --> kilometers for easting and northing

        **plot_yn** : [ 'y' | 'n' ]
                      *'y' to plot on creating an instance

                      *'n' to not plot on creating an instance

        **fig_num** : int
                     figure number.  *Default* is 1

        **title** : string
                    figure title

        **dpi** : int
                  dots per inch of the resolution. *default* is 300

        **plot_tipper** : [ 'yri' | 'yr' | 'yi' | 'n' ]
                        * 'yri' to plot induction both real and imaginary
                           induction arrows

                        * 'yr' to plot just the real induction arrows

                        * 'yi' to plot the imaginary induction arrows

                        * 'n' to not plot them

                        * *Default* is 'n'

                        **Note: convention is to point towards a conductor but
                        can be changed in arrow_dict['direction']**


        **font_size** : float
                        size of the font that labels the plot, 2 will be added
                        to this number for the axis labels.


        **station_dict** : dictionary
                           font dictionary for station name. Keys can be
                           matplotlib.text properties, common ones are:
                               * 'size'   -> for font size
                               * 'weight' -> for font weight
                               * 'color'  -> for color of font

        **arrow_legend_dict** : dictionary of properties for legend with keys:
                               * 'position' -> placement of arrow legend can be:
                                   - 'upper right'
                                   - 'lower right'
                                   - 'upper left'
                                   - 'lower left'
                               * 'xborderpad'-> padding from x axis
                               * 'yborderpad'-> padding from y axis
                               *'fontpad'   -> padding between arrow and
                                               legend text
                               * 'fontdict'  -> dictionary of font properties



        **reference_point** : tuple (x0,y0)
                              reference point estimate relative distance to.
                              This point will be (0,0) on the map and
                              everything else is referenced to this point

    :Example: ::

        >>> #change the axis label and grid color
        >>> ptmap.ax.set_xlabel('Latitude (deg)')
        >>> ptmap.ax.grid(which='major', color=(.5,1,0))
        >>> ptmap.update_plot()

    Attributes:
    --------------

        -arrow_color_imag         imaginary induction arrow color
        -arrow_color_real         real induction arrow color
        -arrow_direction          directional convention of arrows
        -arrow_head_length        head length of arrows in relative points
        -arrow_head_width         head width of arrows in relative points
        -arrow_legend_fontdict    font dictionary for arrow legend
        -arrow_legend_fontpad     padding between font of legend and arrow
        -arrow_legend_position    legend position see matplotlib.legend
        -arrow_legend_xborderpad  padding between legend and x axis
        -arrow_legend_yborderpad  padding between legend and y axis
        -arrow_lw                 arrow line width
        -arrow_size               scaling factor to make arrows visible
        -arrow_threshold          threshold for plotting arrows anything above
                                  this number will not be plotted

        -ax                   matplotlib.axes instance for the main plot
        -ax2                  matplotlib.axes instance for the color bar
        -cb                   matplotlib.colors.ColorBar instance for color bar
        -cb_orientation       color bar orientation ('vertical' | 'horizontal')
        -cb_position          color bar position (x, y, dx, dy)

        -dpi                  dots-per-inch resolution

        -ellipse_cmap         ellipse color map, see above for options
        -ellipse_colorby      parameter to color ellipse by
        -ellipse_range        (min, max, step) values to color ellipses
        -ellipse_size         scaling factor to make ellipses visible

        -fig                  matplotlib.figure instance for the figure
        -fig_num               number of figure being plotted
        -fig_size              size of figure in inches
        -font_size            font size of axes tick label, axes labels will be
                              font_size + 2

        -ftol                 tolerance to look for matching freq 10% default
        -jj                   index of plot freq

        -mapscale             scale of map

        -mt_list               list of mtplot.MTplot instances containing all
                              the important information for each station

        -plot_freq       freq in Hz to plot
        -plot_reference_point  reference point of map, everything will be
                               measured relative to this point
        -plot_tipper           string to indicate to plot induction arrows

        -plot_xarr            array of x-coordinates for stations
        -plot_yarr            array of y-coordinates for stations
        -plot_yn              plot on instance creation

        -rot_z                rotates the data by this angle assuming North is
                              0 and angle measures clockwise

        -tickstrfmt           format of tick strings
        -title                title of figure
        -tscale               temporal scale of y-axis ('freq' | 'period')
        -xpad                 padding between furthest station in x-direction
                              and the axes edge
        -ypad                 padding between furthes station in y-direction
                              and axes edge

    Methods:
    ----------

        -plot                 plots the pseudo section
        -redraw_plot          on call redraws the plot from scratch
        -save_figure          saves figure to a file of given format
        -update_plot          updates the plot while still active
        -export_params_to_file  writes parameters of the phase tensor and tipper to text files.

    """

    def __init__(self, **kwargs):
        """
        Initialise the object
        :param kwargs: keyword-value pairs
        """
        super(PlotPhaseTensorMaps, self).__init__(**kwargs)

        fn_list = kwargs.pop('fn_list', None)
        z_object_list = kwargs.pop('z_object_list', None)
        tipper_object_list = kwargs.pop('tipper_object_list', None)
        mt_object_list = kwargs.pop('mt_object_list', None)
        res_object_list = kwargs.pop('res_object_list', None)

        # ----set attributes for the class-------------------------
        self.mt_list = mtpl.get_mtlist(fn_list=fn_list,
                                       res_object_list=res_object_list,
                                       z_object_list=z_object_list,
                                       tipper_object_list=tipper_object_list,
                                       mt_object_list=mt_object_list)

        # set the freq to plot
        self.plot_freq = kwargs.pop('plot_freq', 1.0)
        self.ftol = kwargs.pop('ftol', 0.1)
        self.interpolate = kwargs.pop('interpolate', True)
        # read in map scale
        self.mapscale = kwargs.pop('mapscale', 'deg')
        # map background image
        self.background_image = kwargs.pop('background_image', None)
        self.bimg_band = kwargs.pop('bimg_band', None)
        self.bimg_cmap = kwargs.pop('bimg_cmap', 'viridis')

        # --> set the ellipse properties -------------------
        # set default size to 2
        if self.mapscale == 'deg':
            self.ellipse_size = kwargs.pop('ellipse_size', .05)
            self.arrow_size = kwargs.pop('arrow_size', .05)
            self.arrow_head_length = kwargs.pop('arrow_head_length', .005)
            self.arrow_head_width = kwargs.pop('arrow_head_width', .005)
            self.arrow_lw = kwargs.pop('arrow_lw', .0005)
            self.xpad = kwargs.pop('xpad', .05)
            self.ypad = kwargs.pop('xpad', .05)

        elif self.mapscale == 'm':
            self.ellipse_size = kwargs.pop('ellipse_size', 500)
            self.arrow_size = kwargs.pop('arrow_size', 500)
            self.arrow_head_length = kwargs.pop('arrow_head_length', 50)
            self.arrow_head_width = kwargs.pop('arrow_head_width', 50)
            self.arrow_lw = kwargs.pop('arrow_lw', 5)
            self.xpad = kwargs.pop('xpad', 500)
            self.ypad = kwargs.pop('xpad', 500)

        elif self.mapscale == 'km':
            self.ellipse_size = kwargs.pop('ellipse_size', .5)
            self.arrow_size = kwargs.pop('arrow_size', .5)
            self.arrow_head_length = kwargs.pop('arrow_head_length', .05)
            self.arrow_head_width = kwargs.pop('arrow_head_width', .05)
            self.arrow_lw = kwargs.pop('arrow_lw', .005)
            self.xpad = kwargs.pop('xpad', .5)
            self.ypad = kwargs.pop('xpad', .5)

        self.minorticks_on = kwargs.pop('minorticks_on', True)

        # --> set colorbar properties---------------------------------
        # set orientation to horizontal
        cb_dict = kwargs.pop('cb_dict', {})
        try:
            self.cb_orientation = cb_dict['orientation']
        except KeyError:
            self.cb_orientation = 'vertical'

        # set the position to middle outside the plot
        try:
            self.cb_position = cb_dict['position']
        except KeyError:
            self.cb_position = None

        # --> set plot properties ------------------------------
        # set some of the properties as attributes much to Lars' discontent
        self.fig_num = kwargs.pop('fig_num', 1)
        self.plot_num = kwargs.pop('plot_num', 1)
        self.plot_style = kwargs.pop('plot_style', '1')
        self.plot_title = kwargs.pop('plot_title', None)
        self.fig_dpi = kwargs.pop('fig_dpi', 300)

        self.tscale = kwargs.pop('tscale', 'period')
        self.fig_size = kwargs.pop('fig_size', [8, 8])
        self.tickstrfmt = kwargs.pop('tickstrfmt', '%.2f')

        self.font_size = kwargs.pop('font_size', 7)

        self.ref_ax_loc = kwargs.pop('ref_ax_loc', (.85, .1, .1, .1))
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

        # --> set induction arrow properties -------------------------------
        self.plot_tipper = kwargs.pop('plot_tipper', 'n')

        # --> set arrow legend properties -------------------------------
        arrow_legend_dict = kwargs.pop('arrow_legend_dict', {})

        self.arrow_legend_position = arrow_legend_dict.pop('position', 'lower right')

        # set x-border pad
        self.arrow_legend_xborderpad = arrow_legend_dict.pop('xborderpad', .2)

        # set y-border pad
        self.arrow_legend_yborderpad = arrow_legend_dict.pop('yborderpad', .2)

        # set font pad
        self.arrow_legend_fontpad = arrow_legend_dict.pop('fontpad', .05)

        # set font properties
        self.arrow_legend_fontdict = arrow_legend_dict.pop('fontdict', {'size': self.font_size,
                                                                        'weight': 'bold'})

        # --> set a central reference point
        self.plot_reference_point = kwargs.pop('reference_point', (0, 0))

        # --> set station name properties
        station_dict = kwargs.pop('station_dict', None)
        if station_dict is not None:
            self.station_id = station_dict.pop('id', (0, 2))

            # set spacing of station name and ellipse
            self.station_pad = station_dict.pop('pad', .0005)

            # set font properties of the station label
            self.station_font_dict = station_dict.pop('font_dict', {'size': self.font_size,
                                                                    'weight': 'bold'})

        self.plot_yn = kwargs.pop('plot_yn', 'y')
        self.save_fn = kwargs.pop('save_fn', "/c/tmp")

        # By this stage all keyword arguments meant to be set as class properties will have
        # been processed. Popping all class properties that still exist in kwargs
        self.kwargs = kwargs
        for key in vars(self):
            self.kwargs.pop(key, None)

        # --> plot if desired ------------------------
        if self.plot_yn == 'y':
            self.plot()
            # self.save_figure(self.save_fn, file_format='png')

    # ---need to rotate data on setting rotz
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
            mt.rot_z = self._rot_z[ii]

    def _get_rot_z(self):
        return self._rot_z

    rot_z = property(fget=_get_rot_z, fset=_set_rot_z,
                     doc="""rotation angle(s)""")

    # -----------------------------------------------
    # The main plot method for this module
    # -----------------------------------------------
    def plot(self, fig=None, save_path=None, show=True,
             raster_dict={'lons': [], 'lats': [],
                          'vals': [], 'levels': 50, 'cmap': 'rainbow',
                          'cbar_title': 'Arbitrary units',
                          'cbar_position': None}):
        """
        Plots the phase tensor map.
        :param fig: optional figure object
        :param save_path: path to folder for saving plots
        :param show: show plots if True
        :param raster_dict: Plotting of raster data is currently only supported when mapscale='deg'.
                            This parameter is a dictionary of parameters for plotting raster data,
                            on top of which phase tensor data are plotted. 'lons', 'lats' and 'vals'
                            are one dimensional lists (or numpy arrays) for longitudes, latitudes
                            and corresponding values, respectively. 'levels', 'cmap' and 'cbar_title'
                            are the number of levels to be used in the colormap, the colormap and its
                            title, respectively.
        """

        # set position properties for the plot
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = .1
        plt.rcParams['figure.subplot.right'] = .98
        plt.rcParams['figure.subplot.bottom'] = .1
        plt.rcParams['figure.subplot.top'] = .93
        plt.rcParams['figure.subplot.wspace'] = .55
        plt.rcParams['figure.subplot.hspace'] = .70
        # FZ: tweaks to make plot positioned better
        # plt.rcParams['font.size']=self.font_size
        # plt.rcParams['figure.subplot.left']=.1
        # plt.rcParams['figure.subplot.right']=.90
        # plt.rcParams['figure.subplot.bottom']=.2
        # plt.rcParams['figure.subplot.top']=.90
        # plt.rcParams['figure.subplot.wspace']=.70
        # plt.rcParams['figure.subplot.hspace']=.70

        lpfig = None
        lpax = None
        lpax2 = None
        # make figure instance
        if(fig is None):
            self.fig = plt.figure(self.fig_num, figsize=self.fig_size, dpi=self.fig_dpi)
            # self.fig = plt.figure(self.fig_num, dpi=self.fig_dpi)

            # clear the figure if there is already one up
            plt.clf()
            lpfig = self.fig
        else:
            lpfig = fig
        # end if

        lpax = lpfig.add_subplot(1, 1, 1, aspect='equal')

        # plt.locator_params(axis='x', nbins=3)  # control number of ticks in axis (nbins ticks)
        plt.xticks(rotation='vertical')  # FZ: control tick rotation=30 not that good

        # get the reference point
        refpoint = self.plot_reference_point

        # plot raster data if provided and if mapscale is 'deg'
        if(len(raster_dict['lons']) and self.mapscale == 'deg'):
            lons = np.array(raster_dict['lons'])
            lats = np.array(raster_dict['lats'])

            # retain masking if a masked array is passed in
            if type(raster_dict['vals']) == np.ma.core.MaskedArray:
                vals = np.ma.masked_array(raster_dict['vals'])
            else:
                vals = np.array(raster_dict['vals'])

            assert len(lons) == len(lats) == len(vals), 'Lons, Lats and Vals must all have the same length'

            lons -= refpoint[0]
            lats -= refpoint[1]
            levels = raster_dict.pop('levels', 50)
            cmap = raster_dict.pop('cmap', 'rainbow')
            cbar_title = raster_dict.pop('cbar_title', 'Arbitrary Units')
            triangulation = tri.Triangulation(lons, lats)
            cbinfo = lpax.tricontourf(triangulation, vals,
                                      levels=np.linspace(vals.min(), vals.max(), levels),
                                      cmap=cmap)
            if raster_dict['cbar_position'] is not None:
                cbax = self.fig.add_axes(raster_dict['cbar_position'])
            else:
                cbax, kw = mcb.make_axes(lpax,
                                         orientation=self.cb_orientation,
                                         shrink=.35)
            cbar = lpfig.colorbar(cbinfo, cbax)

            if(self.cb_orientation == 'horizontal'):
                cbar.ax.set_xlabel(cbar_title)
                cbar.ax.xaxis.set_label_position('top')
                cbar.ax.xaxis.set_label_coords(.5, 1.3)
            else:
                cbar.ax.set_ylabel(cbar_title, fontsize=self.font_size,
                                   fontweight='bold')
                cbar.ax.yaxis.set_label_position('right')
                cbar.ax.yaxis.set_label_coords(1.25, .5)
                cbar.ax.yaxis.tick_left()
                cbar.ax.tick_params(axis='y', direction='in')

        # end if

        # set some local parameters
        es = float(self.ellipse_size)
        cmap = self.ellipse_cmap
        ckmin = float(self.ellipse_range[0])
        ckmax = float(self.ellipse_range[1])
        try:
            ckstep = float(self.ellipse_range[2])
        except IndexError:
            if cmap == 'mt_seg_bl2wh2rd':
                raise ValueError('Need to input range as (min, max, step)')
            else:
                ckstep = 3
        nseg = float((ckmax - ckmin) / (2 * ckstep))
        ck = self.ellipse_colorby

        # --> set the bounds on the segmented colormap
        if cmap == 'mt_seg_bl2wh2rd':
            bounds = np.arange(ckmin, ckmax + ckstep, ckstep)

        # set tick parameters depending on the mapscale
        if self.mapscale == 'deg':
            self.tickstrfmt = '%.2f'

        elif self.mapscale == 'm' or self.mapscale == 'km':
            self.tickstrfmt = '%.0f'

        # make some empty arrays
        elliplist = []
        latlist = np.zeros(len(self.mt_list))
        lonlist = np.zeros(len(self.mt_list))
        self.plot_xarr = np.zeros(len(self.mt_list))
        self.plot_yarr = np.zeros(len(self.mt_list))

        for ii, mt in enumerate(self.mt_list):

            newZ = None
            newTipper = None
            fidx = 0
            if(self.interpolate):
                newZ, newTipper = mt.interpolate([self.plot_freq], bounds_error=False)
            else:
                fidx = np.argmin(np.fabs(mt.Z.freq - self.plot_freq))

            if((not self.interpolate and np.fabs(mt.Z.freq[fidx] - self.plot_freq) < self.ftol) or
                    (self.interpolate)):

                self.jj = fidx
                jj = fidx

                # get phase tensor
                if(not self.interpolate):
                    pt = mt.pt
                else:
                    newZ.compute_resistivity_phase()
                    pt = MTpt.PhaseTensor(z_object=newZ)

                # if map scale is lat lon set parameters
                if self.mapscale == 'deg':
                    latlist[ii] = mt.lat
                    lonlist[ii] = mt.lon
                    plotx = mt.lon - refpoint[0]
                    ploty = mt.lat - refpoint[1]

                # if map scale is in meters easting and northing
                elif self.mapscale == 'm':
                    east, north, zone = gis_tools.project_point_ll2utm(mt.lat,
                                                                       mt.lon)

                    # set the first point read in as a refernce other points
                    if ii == 0:
                        zone1 = zone
                        plotx = east - refpoint[0]
                        ploty = north - refpoint[1]

                    # read in all the other point
                    else:
                        # check to make sure the zone is the same this needs
                        # to be more rigorously done
                        if zone1 != zone:
                            print('Zone change at station ' + mt.station)
                            if zone1[0:2] == zone[0:2]:
                                pass
                            elif int(zone1[0:2]) < int(zone[0:2]):
                                east += 500000
                            else:
                                east -= -500000
                            latlist[ii] = north - refpoint[1]
                            lonlist[ii] = east - refpoint[0]
                            plotx = east - refpoint[0]
                            ploty = north - refpoint[1]
                        else:
                            latlist[ii] = north - refpoint[1]
                            lonlist[ii] = east - refpoint[0]
                            plotx = east - refpoint[0]
                            ploty = north - refpoint[1]

                # if mapscale is in km easting and northing
                elif self.mapscale == 'km':
                    east, north, zone = gis_tools.project_point_ll2utm(mt.lat,
                                                                       mt.lon)
                    if ii == 0:
                        zone1 = zone
                        plotx = (east - refpoint[0]) / 1000.
                        ploty = (north - refpoint[1]) / 1000.

                    else:
                        if zone1 != zone:
                            print('Zone change at station ' + mt.station)
                            if zone1[0:2] == zone[0:2]:
                                pass
                            elif int(zone1[0:2]) < int(zone[0:2]):
                                east += 500000
                            else:
                                east -= 500000
                            latlist[ii] = (north - refpoint[1]) / 1000.
                            lonlist[ii] = (east - refpoint[0]) / 1000.
                            plotx = (east - refpoint[0]) / 1000.
                            ploty = (north - refpoint[1]) / 1000.
                        else:
                            latlist[ii] = (north - refpoint[1]) / 1000.
                            lonlist[ii] = (east - refpoint[0]) / 1000.
                            plotx = (east - refpoint[0]) / 1000.
                            ploty = (north - refpoint[1]) / 1000.
                else:
                    raise NameError('mapscale not recognized')

                # put the location of each ellipse into an array in x and y
                self.plot_xarr[ii] = plotx
                self.plot_yarr[ii] = ploty

                # --> set local variables
                phimin = np.nan_to_num(pt.phimin[jj])
                phimax = np.nan_to_num(pt.phimax[jj])
                eangle = np.nan_to_num(pt.azimuth[jj])

                # output to csv file:
                # print('OUTCSV', mt.station, plotx, ploty, phimin, phimax, eangle)

                if cmap == 'mt_seg_bl2wh2rd':
                    bounds = np.arange(ckmin, ckmax + ckstep, ckstep)
                    nseg = float((ckmax - ckmin) / (2 * ckstep))

                # get the properties to color the ellipses by
                if self.ellipse_colorby == 'phiminang' or \
                        self.ellipse_colorby == 'phimin':
                    colorarray = pt.phimin[jj]

                elif self.ellipse_colorby == 'phimax':
                    colorarray = pt.phimax[jj]

                elif self.ellipse_colorby == 'phidet':
                    colorarray = np.sqrt(abs(pt.det[jj])) * (180 / np.pi)

                elif self.ellipse_colorby == 'skew' or \
                        self.ellipse_colorby == 'skew_seg':
                    colorarray = pt.beta[jj]

                elif self.ellipse_colorby == 'normalized_skew' or \
                        self.ellipse_colorby == 'normalized_skew_seg':
                    colorarray = 2 * pt.beta[jj]

                elif self.ellipse_colorby == 'ellipticity':
                    colorarray = pt.ellipticity[jj]

                else:
                    raise NameError(self.ellipse_colorby + ' is not supported')

                # --> get ellipse properties
                # if the ellipse size is not physically correct make it a dot
                if phimax == 0 or phimax > 100 or phimin == 0 or phimin > 100:
                    eheight = .0000001 * es
                    ewidth = .0000001 * es
                else:
                    scaling = es / phimax
                    eheight = phimin * scaling
                    ewidth = phimax * scaling

                # make an ellipse
                ellipd = patches.Ellipse((plotx, ploty),
                                         width=ewidth,
                                         height=eheight,
                                         angle=90 - eangle,
                                         lw=self.lw,
                                         **self.kwargs)

                # get ellipse color
                if cmap.find('seg') > 0:
                    ellipd.set_facecolor(mtcl.get_plot_color(colorarray,
                                                             self.ellipse_colorby,
                                                             cmap,
                                                             ckmin,
                                                             ckmax,
                                                             bounds=bounds))
                else:
                    ellipd.set_facecolor(mtcl.get_plot_color(colorarray,
                                                             self.ellipse_colorby,
                                                             cmap,
                                                             ckmin,
                                                             ckmax))

                # ==> add ellipse to the plot
                elliplist.append(ellipd)
                lpax.add_artist(ellipd)

                # -----------Plot Induction Arrows---------------------------
                if self.plot_tipper.find('y') == 0:

                    # get tipper

                    if mt.Tipper.tipper is None:
                        mt.Tipper.tipper = np.zeros((len(mt.period), 1, 2),
                                                    dtype='complex')

                    # make some local parameters for easier typing
                    ascale = self.arrow_size
                    adir = self.arrow_direction * np.pi

                    # plot real tipper
                    ti = None
                    if(not self.interpolate):
                        ti = mt.Tipper
                    else:
                        ti = newTipper
                    # end if

                    if self.plot_tipper == 'yri' or self.plot_tipper == 'yr':
                        if ti.mag_real[jj] <= self.arrow_threshold:
                            txr = ti.mag_real[jj] * ascale * \
                                np.sin((ti.angle_real[jj]) * np.pi / 180 + adir)
                            tyr = ti.mag_real[jj] * ascale * \
                                np.cos((ti.angle_real[jj]) * np.pi / 180 + adir)

                            lpax.arrow(plotx,
                                       ploty,
                                       txr,
                                       tyr,
                                       width=self.arrow_lw,
                                       facecolor=self.arrow_color_real,
                                       edgecolor=self.arrow_color_real,
                                       length_includes_head=False,
                                       head_width=self.arrow_head_width,
                                       head_length=self.arrow_head_length)
                        else:
                            pass

                    # plot imaginary tipper
                    if self.plot_tipper == 'yri' or self.plot_tipper == 'yi':
                        if ti.mag_imag[jj] <= self.arrow_threshold:
                            txi = ti.mag_imag[jj] * ascale * \
                                np.sin((ti.angle_imag[jj]) * np.pi / 180 + adir)
                            tyi = ti.mag_imag[jj] * ascale * \
                                np.cos((ti.angle_imag[jj]) * np.pi / 180 + adir)

                            lpax.arrow(plotx,
                                       ploty,
                                       txi,
                                       tyi,
                                       width=self.arrow_lw,
                                       facecolor=self.arrow_color_imag,
                                       edgecolor=self.arrow_color_imag,
                                       length_includes_head=False,
                                       head_width=self.arrow_head_width,
                                       head_length=self.arrow_head_length)

                # ------------Plot station name------------------------------
                try:
                    lpax.text(plotx,
                              ploty + self.station_pad,
                              mt.station[self.station_id[0]:self.station_id[1]],
                              horizontalalignment='center',
                              verticalalignment='baseline',
                              fontdict=self.station_font_dict)
                except AttributeError:
                    pass

            # ==> print a message if couldn't find the freq
            else:
                _logger.warn('Did not find {0:.5g} Hz for station {1}'.format(self.plot_freq, mt.station))

        # --> set axes properties depending on map scale------------------------
        if self.mapscale == 'deg':
            lpax.set_xlabel('Longitude',
                            fontsize=self.font_size,  # +2,
                            fontweight='bold')
            lpax.set_ylabel('Latitude',
                            fontsize=self.font_size,  # +2,
                            fontweight='bold')

        elif self.mapscale == 'm':
            lpax.set_xlabel('Easting (m)',
                            fontsize=self.font_size,  # +2,
                            fontweight='bold')
            lpax.set_ylabel('Northing (m)',
                            fontsize=self.font_size,  # +2,
                            fontweight='bold')

        elif self.mapscale == 'km':
            lpax.set_xlabel('Easting (km)',
                            fontsize=self.font_size,  # +2,
                            fontweight='bold')
            lpax.set_ylabel('Northing (km)',
                            fontsize=self.font_size,  # +2,
                            fontweight='bold')

        # --> set plot limits
        #    need to exclude zero values from the calculation of min/max!!!!
        lpax.set_xlim(self.plot_xarr[self.plot_xarr != 0.].min() - self.xpad,
                      self.plot_xarr[self.plot_xarr != 0.].max() + self.xpad)
        lpax.set_ylim(self.plot_yarr[self.plot_yarr != 0.].min() - self.xpad,
                      self.plot_yarr[self.plot_xarr != 0.].max() + self.xpad)

        # BM: Now that we have the bounds of the axis, we can plot a
        # background image on the map.
        if self.background_image:
            plot_geotiff_on_axes(self.background_image, lpax,
                                 epsg_code=4326, band_number=self.bimg_band,
                                 cmap=self.bimg_cmap)

        # --> set tick label format
        lpax.xaxis.set_major_formatter(FormatStrFormatter(self.tickstrfmt))
        lpax.yaxis.set_major_formatter(FormatStrFormatter(self.tickstrfmt))
#       lpax.set_xticklabels(np.round(self.plot_xarr, decimals=2),
#                                rotation=45)
        plt.setp(lpax.get_xticklabels(), rotation=45)

        # --> set title in period or freq
        if self.tscale == 'period':
            titlefreq = '{0:.5g} (s)'.format(1. / self.plot_freq)
        else:
            titlefreq = '{0:.5g} (Hz)'.format(self.plot_freq)

        if not self.plot_title:
            lpax.set_title('Phase Tensor Map for ' + titlefreq,
                           fontsize=self.font_size + 2, fontweight='bold')
        else:
            lpax.set_title(self.plot_title + titlefreq,
                           fontsize=self.font_size + 2, fontweight='bold')

        # --> plot induction arrow scale bar -----------------------------------
        if self.plot_tipper.find('y') == 0:
            parrx = lpax.get_xlim()
            parry = lpax.get_ylim()
            try:
                axpad = self.arrow_legend_xborderpad
            except AttributeError:
                axpad = self.xpad + self.arrow_size

            try:
                aypad = self.arrow_legend_yborderpad
            except AttributeError:
                aypad = self.ypad

            try:
                txtpad = self.arrow_legend_fontpad
            except AttributeError:
                txtpad = .25 * es

            # make arrow legend postion and arrows coordinates
            if self.arrow_legend_position == 'lower right':
                pax = parrx[1] - axpad
                pay = parry[0] + aypad
                # ptx = self.arrow_size
                # pty = 0
                txa = parrx[1] - axpad + self.arrow_size / 2.
                # txy = pay+txtpad

            elif self.arrow_legend_position == 'upper right':
                pax = parrx[1] - axpad
                pay = parry[1] - aypad
                # ptx = self.arrow_size
                # pty = 0
                txa = parrx[1] - axpad + self.arrow_size / 2.
                # txy = pay+txtpad

            elif self.arrow_legend_position == 'lower left':
                pax = parrx[0] + axpad
                pay = parry[0] + aypad
                # ptx = self.arrow_size
                # pty = 0
                txa = parrx[0] + axpad + self.arrow_size / 2.
                # txy = pay+txtpad

            elif self.arrow_legend_position == 'upper left':
                pax = parrx[0] + axpad
                pay = parry[1] - aypad
                # ptx = self.arrow_size
                # pty = 0
                txa = parrx[0] + axpad + self.arrow_size / 2.
                # txy = pay+txtpad
            else:
                pass  # raise NameError('arrowlegend not supported.')

            # FZ: Sudpip?
            ptx = self.arrow_size
            pty = 0
            # txy = pay + txtpad

            lpax.arrow(pax,
                       pay,
                       ptx,
                       pty,
                       width=self.arrow_lw,
                       facecolor=self.arrow_color_real,
                       edgecolor=self.arrow_color_real,
                       length_includes_head=False,
                       head_width=self.arrow_head_width,
                       head_length=self.arrow_head_length)

            # FZ: what is this '|T|=1'? and the horizontal line?
            # lpax.text(txa,
            #              txy,
            #              '|T|=1',
            #              horizontalalignment='center',
            #              verticalalignment='baseline',
            #              fontdict={'size':self.font_size,'weight':'bold'})
        # END: if self.plot_tipper.find('yes') == 0 ---------------------------

        # make a grid with color lines
        lpax.grid(True, alpha=.3, which='both', color=(0.5, 0.5, 0.5))
        if self.minorticks_on:
            plt.minorticks_on()  # turn on minor ticks automatically

        # ==> make a colorbar with appropriate colors
        if self.cb_position is None:
            lpax2, kw = mcb.make_axes(lpax,
                                      orientation=self.cb_orientation,
                                      shrink=.35)
            # FZ: try to fix colorbar h-position
            # from mpl_toolkits.axes_grid1 import make_axes_locatable
            #
            # # create an axes on the right side of ax. The width of cax will be 5%
            # # of ax and the padding between cax and ax will be fixed at 0.05 inch.
            # divider = make_axes_locatable(self.ax)
            # self.ax2 = divider.append_axes("right", size="5%", pad=0.05)

        else:
            lpax2 = lpfig.add_axes(self.cb_position)

        if cmap == 'mt_seg_bl2wh2rd':
            # make a color list
            self.clist = [(cc, cc, 1)
                          for cc in np.arange(0, 1 + 1. / (nseg), 1. / (nseg))] + \
                         [(1, cc, cc)
                          for cc in np.arange(1, -1. / (nseg), -1. / (nseg))]

            # make segmented colormap
            mt_seg_bl2wh2rd = colors.ListedColormap(self.clist)

            # make bounds so that the middle is white
            bounds = np.arange(ckmin - ckstep, ckmax + 2 * ckstep, ckstep)

            # normalize the colors
            norms = colors.BoundaryNorm(bounds, mt_seg_bl2wh2rd.N)

            # make the colorbar
            self.cb = mcb.ColorbarBase(lpax2,
                                       cmap=mt_seg_bl2wh2rd,
                                       norm=norms,
                                       orientation=self.cb_orientation,
                                       ticks=bounds[1:-1])
        else:
            if cmap in list(mtcl.cmapdict.keys()):
                cmap_input = mtcl.cmapdict[cmap]
            else:
                cmap_input = mtcl.cm.get_cmap(cmap)
            self.cb = mcb.ColorbarBase(lpax2,
                                       cmap=cmap_input,  # mtcl.cmapdict[cmap],
                                       norm=colors.Normalize(vmin=ckmin,
                                                             vmax=ckmax),
                                       orientation=self.cb_orientation)

        # label the color bar accordingly
        self.cb.set_label(mtpl.ckdict[ck],
                          fontdict={'size': self.font_size, 'weight': 'bold'})

        # place the label in the correct location
        if self.cb_orientation == 'horizontal':
            self.cb.ax.xaxis.set_label_position('top')
            self.cb.ax.xaxis.set_label_coords(.5, 1.3)

        elif self.cb_orientation == 'vertical':
            self.cb.ax.yaxis.set_label_position('right')
            self.cb.ax.yaxis.set_label_coords(1.25, .5)
            self.cb.ax.yaxis.tick_left()
            self.cb.ax.tick_params(axis='y', direction='in')

        # --> add reference ellipse:  (legend of ellipse size=1)
        # FZ: remove the following section if no show of Phi
        show_phi = False  # JingMingDuan does not want to show the black circle - it's not useful
        if show_phi is True:
            ref_ellip = patches.Ellipse((0, .0),
                                        width=es,
                                        height=es,
                                        angle=0)
            ref_ellip.set_facecolor((0, 0, 0))
            ref_ax_loc = list(lpax2.get_position().bounds)
            ref_ax_loc[0] *= .95
            ref_ax_loc[1] -= .17
            ref_ax_loc[2] = .1
            ref_ax_loc[3] = .1
            self.ref_ax = lpfig.add_axes(ref_ax_loc, aspect='equal')
            self.ref_ax.add_artist(ref_ellip)
            self.ref_ax.set_xlim(-es / 2. * 1.05, es / 2. * 1.05)
            self.ref_ax.set_ylim(-es / 2. * 1.05, es / 2. * 1.05)
            plt.setp(self.ref_ax.xaxis.get_ticklabels(), visible=False)
            plt.setp(self.ref_ax.yaxis.get_ticklabels(), visible=False)
            self.ref_ax.set_title(r'$\Phi$ = 1')

        if show:
            # always show, and adjust the figure before saving it below. The
            # figure size ratio are all different!!
            plt.show()

        # the figure need to be closed (X) then the following code save it to a
        # file.
        if save_path is not None:
            figfile = self.save_figure(save_path)  # , fig_dpi=300)
        else:
            figfile = None

        # self.export_params_to_file('E:/tmp')

        return figfile

    def save_figure(self, save_fn, file_format='pdf',
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

            **file_format** : [ jpg | png | pdf | eps | svg ]
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
        """

        sf = '_{0:.6g}'.format(self.plot_freq)

        if fig_dpi is None:
            fig_dpi = self.fig_dpi

        # FZ: fixed the following logic
        if os.path.isdir(save_fn):  # FZ: assume save-fn is a directory
            if not os.path.exists(save_fn):
                os.mkdir(save_fn)

            # make a file name
            fname = 'PTmap_DPI%s_%s_%sHz.%s' % (
                str(self.fig_dpi), self.ellipse_colorby, sf, file_format)
            path2savefile = os.path.join(save_fn, fname)
#            self.fig.savefig(path2savefile, dpi=fig_dpi, format=file_format, orientation=orientation,
#                             bbox_inches='tight')
        else:  # FZ: assume save-fn is a path2file= "path2/afile.fmt"
            file_format = save_fn.split('.')[-1]
            if file_format is None or file_format not in ['png', 'jpg']:
                _logger.error(
                    "Error: output file name is not correctly provided:",
                    save_fn)
                raise Exception(
                    "output file name is not correctly provided!!!")

            path2savefile = save_fn
            self.fig.savefig(path2savefile, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')
            # plt.clf()
            # plt.close(self.fig)

        if close_plot == 'y':
            plt.clf()
            plt.close(self.fig)
        else:
            pass

        self.fig_fn = save_fn
        _logger.debug('Saved figure to: %s', self.fig_fn)

    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.

        Use this if you change an of the .fig or axes properties
        """

        self.fig.canvas.draw()

    def redraw_plot(self):
        """
        use this function if you updated some attributes and want to re-plot.
        """

        plt.close(self.fig)
        self.plot()

    def export_params_to_file(self, save_path=None):
        """
        write text files for all the phase tensor parameters.
        :param save_path: string path to save files into.
        File naming pattern is like save_path/PhaseTensorTipper_Params_freq.csv/table
        **Files Content **
                *station
                *lon
                *lat
                *phi_min
                *phi_max
                *skew
                *ellipticity
                *azimuth
                *tipper_mag_real
                *tipper_ang_real
                *tipper_mag_imag
                *tipper_ang_imag

        :return: path2savedfile
        """

        # create a save path for files
        if save_path is None:
            _logger.info("No save_path provided. ")
            return None
            # try:
            #     svpath = os.path.join(os.path.dirname(self.mt_list[0].fn),
            #                           'PTMaps')
            # except TypeError:
            #     raise IOError('Need to input save_path, could not find path')
        else:
            svpath = save_path

        # if the folder doesn't exist make it
        if not os.path.exists(svpath):
            os.mkdir(svpath)

        # make sure the attributes are there if not get them
        try:
            self.plot_xarr
        except AttributeError:
            self.plot(show=False)

        # sort the x and y in ascending order to get placement in the file right
        xlist = np.sort(abs(self.plot_xarr))
        ylist = np.sort(abs(self.plot_yarr))

        # get the indicies of where the values should go in map view of the
        # text file
        nx = self.plot_xarr.shape[0]
        xyloc = np.zeros((nx, 2), np.uint)
        for jj, xx in enumerate(self.plot_xarr):
            xyloc[jj, 0] = np.where(xlist == abs(xx))[0][0]
            xyloc[jj, 1] = np.where(ylist == abs(self.plot_yarr[jj]))[0][0]

        # create arrays that simulate map view in a text file
        phiminmap = np.zeros((xlist.shape[0], ylist.shape[0]))
        phimaxmap = np.zeros((xlist.shape[0], ylist.shape[0]))
        azimuthmap = np.zeros((xlist.shape[0], ylist.shape[0]))
        ellipmap = np.zeros((xlist.shape[0], ylist.shape[0]))
        betamap = np.zeros((xlist.shape[0], ylist.shape[0]))
        trmap = np.zeros((xlist.shape[0], ylist.shape[0]))
        trazmap = np.zeros((xlist.shape[0], ylist.shape[0]))
        timap = np.zeros((xlist.shape[0], ylist.shape[0]))
        tiazmap = np.zeros((xlist.shape[0], ylist.shape[0]))
        stationmap = np.zeros((xlist.shape[0], ylist.shape[0]), dtype='|S8')

        station_location = {}  # a dict to store all MT stations (lan lat)

        # put the information into the zeroed arrays
        for ii in range(nx):
            mt1 = self.mt_list[ii]

            # try to find the freq in the freq list of each file
            freqfind = [ff for ff, f2 in enumerate(mt1.Z.freq)
                        if self.plot_freq * (1 - self.ftol) < f2 < self.plot_freq * (1 + self.ftol)]
            try:
                j2 = freqfind[0]
                # freq0 = mt1.Z.freq[j2]  # should use the closest freq from this list found within the tolerance range

                pt = mt1.pt
                tp = mt1.Tipper

                if pt.phimin[0] is not None:
                    phiminmap[xyloc[ii, 0], xyloc[ii, 1]] = pt.phimin[0][self.jj]
                    phimaxmap[xyloc[ii, 0], xyloc[ii, 1]] = pt.phimax[0][self.jj]
                    azimuthmap[xyloc[ii, 0], xyloc[ii, 1]] = pt.azimuth[0][self.jj]
                    ellipmap[xyloc[ii, 0], xyloc[ii, 1]] = pt.ellipticity[0][self.jj]
                    betamap[xyloc[ii, 0], xyloc[ii, 1]] = pt.beta[0][self.jj]

                if tp.mag_real is not None:
                    trmap[xyloc[ii, 0], xyloc[ii, 1]] = tp.mag_real[self.jj]
                    trazmap[xyloc[ii, 0], xyloc[ii, 1]] = tp.angle_real[self.jj]
                    timap[xyloc[ii, 0], xyloc[ii, 1]] = tp.mag_imag[self.jj]
                    tiazmap[xyloc[ii, 0], xyloc[ii, 1]] = tp.angle_imag[self.jj]
                try:
                    stationmap[xyloc[ii, 0], xyloc[ii, 1]] = \
                        mt1.station[self.station_id[0]:self.station_id[1]]
                except AttributeError:
                    stationmap[xyloc[ii, 0], xyloc[ii, 1]] = mt1.station

                station_location[stationmap[xyloc[ii, 0], xyloc[ii, 1]]] = (mt1.lon, mt1.lat, mt1.freq[j2])
            except IndexError:
                _logger.warn('Did not find {0:.5g} Hz for station {1}'.format(self.plot_freq, mt1.station))

        # ----------------------write files-------------------------------------
        svfn = 'Map_{0:.6g}Hz'.format(self.plot_freq)
        ptminfid = open(os.path.join(svpath, svfn + '.phimin'), 'w')
        ptmaxfid = open(os.path.join(svpath, svfn + '.phimax'), 'w')
        ptazmfid = open(os.path.join(svpath, svfn + '.azimuth'), 'w')
        ptskwfid = open(os.path.join(svpath, svfn + '.skew'), 'w')
        ptellfid = open(os.path.join(svpath, svfn + '.ellipticity'), 'w')
        tprmgfid = open(os.path.join(svpath, svfn + '.tipper_mag_real'), 'w')
        tprazfid = open(os.path.join(svpath, svfn + '.tipper_ang_real'), 'w')
        tpimgfid = open(os.path.join(svpath, svfn + '.tipper_mag_imag'), 'w')
        tpiazfid = open(os.path.join(svpath, svfn + '.tipper_ang_imag'), 'w')
        statnfid = open(os.path.join(svpath, svfn + '.station'), 'w')
        tablefid = open(os.path.join(svpath, svfn + '.table'), 'w')

        for ly in range(ylist.shape[0]):
            for lx in range(xlist.shape[0]):
                # if there is nothing there write some spaces
                if phiminmap[lx, ly] == 0.0:
                    ptminfid.write('{0:^8}'.format(' '))
                    ptmaxfid.write('{0:^8}'.format(' '))
                    ptazmfid.write('{0:^8}'.format(' '))
                    ptskwfid.write('{0:^8}'.format(' '))
                    ptellfid.write('{0:^8}'.format(' '))
                    tprmgfid.write('{0:^8}'.format(' '))
                    tprazfid.write('{0:^8}'.format(' '))
                    tpimgfid.write('{0:^8}'.format(' '))
                    tpiazfid.write('{0:^8}'.format(' '))
                    statnfid.write('{0:^8}'.format(' '))

                else:
                    ptminfid.write(mtpl.make_value_str(phiminmap[lx, ly]))
                    ptmaxfid.write(mtpl.make_value_str(phimaxmap[lx, ly]))
                    ptazmfid.write(mtpl.make_value_str(azimuthmap[lx, ly]))
                    ptskwfid.write(mtpl.make_value_str(betamap[lx, ly]))
                    ptellfid.write(mtpl.make_value_str(ellipmap[lx, ly]))
                    tprmgfid.write(mtpl.make_value_str(trmap[lx, ly]))
                    tprazfid.write(mtpl.make_value_str(trazmap[lx, ly]))
                    tpimgfid.write(mtpl.make_value_str(timap[lx, ly]))
                    tpiazfid.write(mtpl.make_value_str(tiazmap[lx, ly]))
                    statnfid.write('{0:^8}'.format(stationmap[lx, ly]))

            # make sure there is an end of line
            ptminfid.write('\n')
            ptmaxfid.write('\n')
            ptazmfid.write('\n')
            ptskwfid.write('\n')
            ptellfid.write('\n')
            tprmgfid.write('\n')
            tprazfid.write('\n')
            tpimgfid.write('\n')
            tpiazfid.write('\n')
            statnfid.write('\n')

        # close the files
        ptminfid.close()
        ptmaxfid.close()
        ptazmfid.close()
        ptskwfid.close()
        ptellfid.close()
        tprmgfid.close()
        tprazfid.close()
        tpimgfid.close()
        tpiazfid.close()
        statnfid.close()

        # --> write the table file
        # write header
        for ss in ['station', 'lon', 'lat', 'phi_min', 'phi_max', 'skew', 'ellipticity',
                   'azimuth', 'tip_mag_re', 'tip_ang_re', 'tip_mag_im',
                   'tip_ang_im', 'frequency']:
            tablefid.write('{0:^12}'.format(ss))
        tablefid.write('\n')

        for ii in range(nx):
            xx, yy = xyloc[ii, 0], xyloc[ii, 1]
            if stationmap[xx, yy] is None or len(stationmap[xx, yy]) < 1:
                pass  # station not having the freq
            else:  # only those stations with the given freq
                # Station
                tablefid.write('{0:^12}'.format(stationmap[xx, yy]))
                # lon
                tablefid.write(mtpl.make_value_str(station_location[stationmap[xx, yy]][0], spacing='{0:^12}'))
                # lat
                tablefid.write(mtpl.make_value_str(station_location[stationmap[xx, yy]][1], spacing='{0:^12}'))
                # phi_min
                tablefid.write(mtpl.make_value_str(phiminmap[xx, yy],
                                                   spacing='{0:^12}'))
                # phi_max
                tablefid.write(mtpl.make_value_str(phimaxmap[xx, yy],
                                                   spacing='{0:^12}'))
                # beta_skew
                tablefid.write(mtpl.make_value_str(betamap[xx, yy],
                                                   spacing='{0:^12}'))
                # ellip
                tablefid.write(mtpl.make_value_str(ellipmap[xx, yy],
                                                   spacing='{0:^12}'))
                # azimuth
                tablefid.write(mtpl.make_value_str(azimuthmap[xx, yy],
                                                   spacing='{0:^12}'))
                # tiprmag
                tablefid.write(mtpl.make_value_str(trmap[xx, yy],
                                                   spacing='{0:^12}'))
                # tiprang
                tablefid.write(mtpl.make_value_str(trazmap[xx, yy],
                                                   spacing='{0:^12}'))
                # tipimag
                tablefid.write(mtpl.make_value_str(timap[xx, yy],
                                                   spacing='{0:^12}'))
                # tipiang
                tablefid.write(mtpl.make_value_str(tiazmap[xx, yy],
                                                   spacing='{0:^12}'))
                # frequency
                tablefid.write(mtpl.make_value_str(station_location[stationmap[xx, yy]][2], spacing='{0:^12}'))
                tablefid.write('\n')

        tablefid.write('\n')

        _logger.info('Wrote files to {}'.format(svpath))

        return svpath

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return "Plots phase tensor maps for one freq"


# ====================================================================
# How to test use this module
# User modify scripts below to provide three things, 1) edidir the path to an edi folder, 2) plot_freq 3) output_dir
# For more advanced commandline script see: examples/cmdline
#  =======================================================================
if __name__ == "__main__":
    imaging = os.path.dirname(__file__)
    mtpy = os.path.dirname(imaging)
    base = os.path.dirname(mtpy)
    examples = os.path.join(base, 'examples')
    data = os.path.join(examples, 'data')
    edidir = os.path.join(data, 'edi_files_2')

    edi_file_list = glob.glob(edidir + '/*.edi')
    savedir = '/tmp'

    plot_freq = 1e-2
    ptm_obj = PlotPhaseTensorMaps(fn_list=edi_file_list,
                                  plot_freq=plot_freq,
                                  ftol=.2,
                                  interpolate=True,
                                  xpad=0.02,
                                  plot_tipper='yr',
                                  edgecolor='k',
                                  lw=0.1,
                                  alpha=1,
                                  minorticks_on=False,
                                  ellipse_size=.2,
                                  ellipse_range=[-10, 10, 2],
                                  ellipse_colorby='skew',
                                  arrow_size=0.5,
                                  ellipse_cmap='mt_seg_bl2wh2rd',
                                  plot_yn='n')

    # generate raster data
    lons, lats = np.meshgrid(np.linspace(136, 140, 50), np.linspace(-19, -22, 50))
    vals = np.exp(-np.sin(np.radians(lons) * 50)**2 - np.cos(np.radians(lats) * 70)**3)

    # plot with raster data
    f = plt.figure(figsize=(8, 6))
    ptm_obj.plot(fig=f, show=True,
                 raster_dict={'lons': lons.flatten(),
                              'lats': lats.flatten(),
                              'vals': vals.flatten(),
                              'levels': 50,
                              'cmap': 'rainbow',
                              'cbar_title': 'Arbitrary Units'})

    ptm_obj.export_params_to_file(save_path=savedir)
