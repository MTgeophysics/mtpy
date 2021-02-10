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
from mtpy.utils.mtpy_logger import get_mtpy_logger
from mtpy.utils.plot_geotiff_imshow import plot_geotiff_on_axes

_logger = get_mtpy_logger(__name__)


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
        
        self.fig = None
        self.ax_ptm = None
        self.ax_cb = None

        fn_list = kwargs.pop("fn_list", None)
        z_object_list = kwargs.pop("z_object_list", None)
        tipper_object_list = kwargs.pop("tipper_object_list", None)
        mt_object_list = kwargs.pop("mt_object_list", None)
        res_object_list = kwargs.pop("res_object_list", None)

        # ----set attributes for the class-------------------------
        self.mt_list = mtpl.get_mtlist(
            fn_list=fn_list,
            res_object_list=res_object_list,
            z_object_list=z_object_list,
            tipper_object_list=tipper_object_list,
            mt_object_list=mt_object_list,
        )

        # set the freq to plot
        self.plot_freq =  1.0
        self.ftol = 0.1
        self.interpolate = True
        # read in map scale
        self.mapscale = 'deg'
        # map background image
        self.background_image = None
        self.bimg_band = None
        self.bimg_cmap = 'viridis'

        # --> set the ellipse properties -------------------
        # set default size to 2
        if self.mapscale == 'deg':
            self.ellipse_size = .05
            self.arrow_size = .05
            self.arrow_head_length = .005
            self.arrow_head_width = .005
            self.arrow_lw = .0005
            self.xpad = .05
            self.ypad = .05

        elif self.mapscale == 'm':
            self.ellipse_size = 500
            self.arrow_size = 500
            self.arrow_head_length = 50
            self.arrow_head_width = 50
            self.arrow_lw = 5
            self.xpad = 500
            self.ypad = 500

        elif self.mapscale == 'km':
            self.ellipse_size = .5
            self.arrow_size = .5
            self.arrow_head_length = .05
            self.arrow_head_width = .05
            self.arrow_lw = .005
            self.xpad = .5
            self.ypad = .5

        self.minorticks_on = True

        # --> set colorbar properties---------------------------------
        # set orientation to horizontal
        self.cb_orientation = 'vertical'
        self.cb_position = None

        # --> set plot properties ------------------------------
        # set some of the properties as attributes much to Lars' discontent
        self.fig_num = 1
        self.plot_title = 'Phase Tensor Map for '
        self.fig_dpi = 300

        self.tscale = 'period'
        self.fig_size = None
        self.tickstrfmt = '%.2f'

        self.font_size = 7

        self.ref_ax_loc = (.85, .1, .1, .1)
        # if rotation angle is an int or float make an array the length of
        # mt_list for plotting purposes

        # --> set induction arrow properties -------------------------------
        self.plot_tipper = 'n'

        # --> set arrow legend properties -------------------------------

        self.arrow_legend_position = 'lower right'
        self.arrow_legend_xborderpad = .2
        self.arrow_legend_yborderpad = .2
        self.arrow_legend_fontpad = .05
        self.arrow_legend_fontdict = {'size': self.font_size, 
                                      'weight': 'bold'}

        # --> set a central reference point
        self.plot_reference_point = (0, 0)

        # --> set station name properties
        self.plot_station = False
        self.station_id = (0, 2)
        

        # set spacing of station name and ellipse
        self.station_pad = .0005

        # set font properties of the station label
        self.station_font_dict = {'size': self.font_size, 'weight': 'bold'}

        self.plot_yn = 'y'
        self.save_fn = "/c/tmp"

        for key, value in kwargs.items():
            setattr(self, key, value)

        # --> plot if desired ------------------------
        if self.plot_yn == "y":
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
        if value == 0:
            return
        for ii, mt_obj in enumerate(self.mt_list):
            mt_obj.rotation_angle = value
            
        self._rotation_angle = value
        
    def fold_strike(self, strike):
        """
        
        """
        strike = strike % 180
        strike[np.where(strike > 90)] -= 180
        
        return strike
    
    def add_raster(self, raster_dict):
        """
        """
        
        # plot raster data if provided and if mapscale is 'deg'
        if len(raster_dict["lons"]) and self.mapscale == "deg":
            lons = np.array(raster_dict["lons"])
            lats = np.array(raster_dict["lats"])

            # retain masking if a masked array is passed in
            if type(raster_dict["vals"]) == np.ma.core.MaskedArray:
                vals = np.ma.masked_array(raster_dict["vals"])
            else:
                vals = np.array(raster_dict["vals"])

            if not len(lons) == len(lats) == len(vals):
                msg= 'Raster Lons, Lats and Vals must all have the same length'
                _logger.error(msg)

            lons -= self.plot_reference_point[0]
            lats -= self.plot_reference_point[1]
            levels = raster_dict.pop('levels', 50)
            cmap = raster_dict.pop('cmap', 'rainbow')
            cbar_title = raster_dict.pop('cbar_title', 'Arbitrary Units')
            triangulation = tri.Triangulation(lons, lats)
            cbinfo = self.ax_ptm.tricontourf(triangulation, vals,
                                      levels=np.linspace(vals.min(), 
                                                         vals.max(), 
                                                         levels),
                                      cmap=cmap)
            if raster_dict['cbar_position'] is not None:
                cbax = self.fig.add_axes(raster_dict['cbar_position'])
            else:
                cbax, kw = mcb.make_axes(self.ax_ptm,
                                         orientation=self.cb_orientation,
                                         shrink=.35)
            cbar = self.fig.colorbar(cbinfo, cbax)

            if self.cb_orientation == "horizontal":
                cbar.ax.set_xlabel(cbar_title)
                cbar.ax.xaxis.set_label_position("top")
                cbar.ax.xaxis.set_label_coords(0.5, 1.3)
            else:
                cbar.ax.set_ylabel(
                    cbar_title, fontsize=self.font_size, fontweight="bold"
                )
                cbar.ax.yaxis.set_label_position("right")
                cbar.ax.yaxis.set_label_coords(1.25, 0.5)
                cbar.ax.yaxis.tick_left()
                cbar.ax.tick_params(axis="y", direction="in")

        # end if
        
    def add_tipper(self, tipper_obj, plot_x, plot_y):
        """
        
        :param tipper: DESCRIPTION
        :type tipper: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        # make some local parameters for easier typing
        ascale = self.arrow_size
        adir = self.arrow_direction * np.pi

        # plot real tipper
        ti = tipper_obj
        if ti is None:
            return

        if self.plot_tipper == 'yri' or self.plot_tipper == 'yr':
            if ti.mag_real[self.jj] <= self.arrow_threshold:
                txr = ti.mag_real[self.jj] * ascale * \
                    np.sin((ti.angle_real[self.jj]) * np.pi / 180 + adir)
                tyr = ti.mag_real[self.jj] * ascale * \
                    np.cos((ti.angle_real[self.jj]) * np.pi / 180 + adir)

                self.ax_ptm.arrow(plot_x,
                                  plot_y,
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
            if ti.mag_imag[self.jj] <= self.arrow_threshold:
                txi = ti.mag_imag[self.jj] * ascale * \
                    np.sin((ti.angle_imag[self.jj]) * np.pi / 180 + adir)
                tyi = ti.mag_imag[self.jj] * ascale * \
                    np.cos((ti.angle_imag[self.jj]) * np.pi / 180 + adir)

                self.ax_ptm.arrow(plot_x,
                                  plot_y,
                                  txi,
                                  tyi,
                                  width=self.arrow_lw,
                                  facecolor=self.arrow_color_imag,
                                  edgecolor=self.arrow_color_imag,
                                  length_includes_head=False,
                                  head_width=self.arrow_head_width,
                                  head_length=self.arrow_head_length)
    
    def get_xy(self, mt_obj):
        """
        """
        # if map scale is lat lon set parameters
        if self.mapscale == 'deg':
            plot_x = mt_obj.lon - self.plot_reference_point[0]
            plot_y = mt_obj.lat - self.plot_reference_point[1]

        # if map scale is in meters easting and northing
        elif self.mapscale in ['km', 'm']:
            if self.mapscale == 'km':
                scale = 1./1000
            else:
                scale = 1.
            east, north, zone = gis_tools.project_point_ll2utm(mt_obj.lat,
                                                               mt_obj.lon)

                
            plot_x *= scale
            plot_y *= scale

        else:
            raise NameError('mapscale not recognized')
            
        return plot_x, plot_y
    
    def get_pt_ellipse(self, pt, plot_x, plot_y):
        """
        """

        # --> set local variables
        phimin = np.nan_to_num(pt.phimin[self.jj])
        phimax = np.nan_to_num(pt.phimax[self.jj])
        eangle = np.nan_to_num(pt.azimuth[self.jj])

        if self.ellipse_cmap == 'mt_seg_bl2wh2rd':
            bounds = np.arange(self.ellipse_range[0], 
                               self.ellipse_range[1] + self.ellipse_range[2],
                               self.ellipse_range[2])
            nseg = float((self.ellipse_range[1] - self.ellipse_range[0]) / \
                         (2 * self.ellipse_range[2]))

        # get the properties to color the ellipses by
        if self.ellipse_colorby == 'phiminang' or \
                self.ellipse_colorby == 'phimin':
            colorarray = pt.phimin[self.jj]

        elif self.ellipse_colorby == 'phimax':
            colorarray = pt.phimax[self.jj]

        elif self.ellipse_colorby == 'phidet':
            colorarray = np.sqrt(abs(pt.det[self.jj])) * (180 / np.pi)

        elif self.ellipse_colorby == 'skew' or \
                self.ellipse_colorby == 'skew_seg':
            colorarray = pt.beta[self.jj]

        elif self.ellipse_colorby == 'normalized_skew' or \
                self.ellipse_colorby == 'normalized_skew_seg':
            colorarray = 2 * pt.beta[self.jj]

        elif self.ellipse_colorby == 'ellipticity':
            colorarray = pt.ellipticity[self.jj]
            
        elif self.ellipse_colorby in ['strike', 'azimuth']:
            colorarray = pt.azimuth[self.jj] % 180
            if colorarray > 90:
                colorarray -= 180
            self.ellipse_range = (-90, 90)

        else:
            raise NameError(self.ellipse_colorby + ' is not supported')

        # --> get ellipse properties
        # if the ellipse size is not physically correct make it a dot
        if phimax == 0 or phimax > 100 or phimin == 0 or phimin > 100:
            eheight = .0000001 * self.ellipse_size
            ewidth = .0000001 * self.ellipse_size
        else:
            scaling = self.ellipse_size / phimax
            eheight = phimin * scaling
            ewidth = phimax * scaling

        # make an ellipse
        ellipd = patches.Ellipse((plot_x, plot_y),
                                 width=ewidth,
                                 height=eheight,
                                 angle=90 - eangle,
                                 lw=self.lw)

        # get ellipse color
        if self.ellipse_cmap.find('seg') > 0:
            ellipd.set_facecolor(mtcl.get_plot_color(colorarray,
                                                     self.ellipse_colorby,
                                                     self.ellipse_cmap,
                                                     self.ellipse_range[0],
                                                     self.ellipse_range[1],
                                                     bounds=bounds))
        else:
            ellipd.set_facecolor(mtcl.get_plot_color(colorarray,
                                                     self.ellipse_colorby,
                                                     self.ellipse_cmap,
                                                     self.ellipse_range[0],
                                                     self.ellipse_range[1]))
            
        return ellipd
        

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
        :param raster_dict: Plotting of raster data is currently only 
                            supported when mapscale='deg'.
                            This parameter is a dictionary of parameters for 
                            plotting raster data,
                            on top of which phase tensor data are plotted.
                            'lons', 'lats' and 'vals'  are one dimensional 
                            lists (or numpy arrays) for longitudes, latitudes
                            and corresponding values, respectively. 'levels', 
                            'cmap' and 'cbar_title' are the number of levels
                            to be used in the colormap, the colormap and its
                            title, respectively.
        """

        # set position properties for the plot
        plt.rcParams['font.size'] = self.font_size

        # make figure instance
        if fig is not None:
            self.fig = fig
        else:
            self.fig = plt.figure(self.fig_num, figsize=self.fig_size,
                                  dpi=self.fig_dpi)
        plt.clf()

        self.ax_ptm = self.fig.add_subplot(1, 1, 1, aspect='equal')

        # plt.locator_params(axis='x', nbins=3)  
        # control number of ticks in axis (nbins ticks)
        # FZ: control tick rotation=30 not that good
        plt.xticks(rotation='vertical')  

        self.add_raster(raster_dict)

        try:
            step = float(self.ellipse_range[2])
        except IndexError:
            if self.ellipse_cmap == 'mt_seg_bl2wh2rd':
                raise ValueError('Need to input range as (min, max, step)')
            else:
                step = 3
        nseg = float((self.ellipse_range[1] - self.ellipse_range[0]) / \
                     (2 * step))
        ck = self.ellipse_colorby

        # --> set the bounds on the segmented colormap
        if self.ellipse_cmap == 'mt_seg_bl2wh2rd':
            bounds = np.arange(self.ellipse_range[0], 
                               self.ellipse_range[1] + step,
                               step)

        # set tick parameters depending on the mapscale
        if self.mapscale == "deg":
            self.tickstrfmt = "%.2f"

        elif self.mapscale == "m" or self.mapscale == "km":
            self.tickstrfmt = "%.0f"

        # make some empty arrays
        elliplist = []
        self.plot_xarr = np.zeros(len(self.mt_list))
        self.plot_yarr = np.zeros(len(self.mt_list))

        for ii, mt_obj in enumerate(self.mt_list):

            new_z_obj = None
            new_tipper_obj = None
            fidx = 0
            if self.interpolate:
                new_z_obj, new_tipper_obj = mt_obj.interpolate([self.plot_freq],
                                                 bounds_error=False)
            else:
                fidx = np.argmin(np.fabs(mt_obj.Z.freq - self.plot_freq))

            if((not self.interpolate and \
                np.fabs(mt_obj.Z.freq[fidx] - self.plot_freq) < self.ftol) or \
               (self.interpolate)):

                self.jj = fidx

                # get phase tensor
                if(not self.interpolate):
                    pt = mt_obj.pt
                else:
                    new_z_obj.compute_resistivity_phase()
                    pt = MTpt.PhaseTensor(z_object=new_z_obj)

                plot_x, plot_y = self.get_xy(mt_obj)
                # put the location of each ellipse into an array in x and y
                self.plot_xarr[ii] = plot_x
                self.plot_yarr[ii] = plot_y
                
                ellipd = self.get_pt_ellipse(pt, plot_x, plot_y)

                # ==> add ellipse to the plot
                elliplist.append(ellipd)
                self.ax_ptm.add_artist(ellipd)

                # -----------Plot Induction Arrows---------------------------
                if self.plot_tipper.find('y') == 0:
                    # get phase tensor
                    if not self.interpolate:
                        new_tipper_obj = mt_obj.Tipper
                    else:
                        new_tipper_obj.compute_mag_direction()
                    
                    if mt_obj.Tipper.tipper is None:
                        continue
                    self.add_tipper(new_tipper_obj, plot_x, plot_y)

                # ------------Plot station name------------------------------
                if self.plot_station:
                    try:
                        name = mt_obj.station[self.station_id[0]:self.station_id[1]]
                        self.ax_ptm.text(plot_x,
                                         plot_y + self.station_pad,
                                         name,
                                         horizontalalignment='center',
                                         verticalalignment='baseline',
                                         fontdict=self.station_font_dict)
                    except AttributeError:
                        pass

            # ==> print a message if couldn't find the freq
            else:
                _logger.warn('Did not find {0:.5g} Hz for station {1}'.format(self.plot_freq, mt_obj.station))

        # --> set axes properties depending on map scale------------------------
        if self.mapscale == 'deg':
            self.ax_ptm.set_xlabel('Longitude',
                            fontsize=self.font_size,  # +2,
                            fontweight='bold')
            self.ax_ptm.set_ylabel('Latitude',
                            fontsize=self.font_size,  # +2,
                            fontweight='bold')

        elif self.mapscale == 'm':
            self.ax_ptm.set_xlabel('Easting (m)',
                            fontsize=self.font_size,  # +2,
                            fontweight='bold')
            self.ax_ptm.set_ylabel('Northing (m)',
                            fontsize=self.font_size,  # +2,
                            fontweight='bold')

        elif self.mapscale == 'km':
            self.ax_ptm.set_xlabel('Easting (km)',
                            fontsize=self.font_size,  # +2,
                            fontweight='bold')
            self.ax_ptm.set_ylabel('Northing (km)',
                            fontsize=self.font_size,  # +2,
                            fontweight='bold')

        # --> set plot limits
        #    need to exclude zero values from the calculation of min/max!!!!
        self.ax_ptm.set_xlim(self.plot_xarr[self.plot_xarr != 0.].min() - self.xpad,
                      self.plot_xarr[self.plot_xarr != 0.].max() + self.xpad)
        self.ax_ptm.set_ylim(self.plot_yarr[self.plot_yarr != 0.].min() - self.xpad,
                      self.plot_yarr[self.plot_xarr != 0.].max() + self.xpad)

        # BM: Now that we have the bounds of the axis, we can plot a
        # background image on the map.
        if self.background_image:
            plot_geotiff_on_axes(self.background_image, self.ax_ptm,
                                 epsg_code=4326, band_number=self.bimg_band,
                                 cmap=self.bimg_cmap)

        # --> set tick label format
        self.ax_ptm.xaxis.set_major_formatter(FormatStrFormatter(self.tickstrfmt))
        self.ax_ptm.yaxis.set_major_formatter(FormatStrFormatter(self.tickstrfmt))
        plt.setp(self.ax_ptm.get_xticklabels(), rotation=45)

        # --> set title in period or freq
        if self.tscale == "period":
            titlefreq = "{0:.5g} (s)".format(1.0 / self.plot_freq)
        else:
            titlefreq = "{0:.5g} (Hz)".format(self.plot_freq)

        self.ax_ptm.set_title('{0} {1}'.format(self.plot_title, titlefreq),
                              fontsize=self.font_size + 2, fontweight='bold')

        # --> plot induction arrow scale bar -----------------------------------
        if self.plot_tipper.find('y') == 0:
            parrx = self.ax_ptm.get_xlim()
            parry = self.ax_ptm.get_ylim()
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
                txtpad = .25 * self.ellipse_size

            # make arrow legend postion and arrows coordinates
            if self.arrow_legend_position == "lower right":
                pax = parrx[1] - axpad
                pay = parry[0] + aypad
                # ptx = self.arrow_size
                # pty = 0
                txa = parrx[1] - axpad + self.arrow_size / 2.0
                # txy = pay+txtpad

            elif self.arrow_legend_position == "upper right":
                pax = parrx[1] - axpad
                pay = parry[1] - aypad
                # ptx = self.arrow_size
                # pty = 0
                txa = parrx[1] - axpad + self.arrow_size / 2.0
                # txy = pay+txtpad

            elif self.arrow_legend_position == "lower left":
                pax = parrx[0] + axpad
                pay = parry[0] + aypad
                # ptx = self.arrow_size
                # pty = 0
                txa = parrx[0] + axpad + self.arrow_size / 2.0
                # txy = pay+txtpad

            elif self.arrow_legend_position == "upper left":
                pax = parrx[0] + axpad
                pay = parry[1] - aypad
                # ptx = self.arrow_size
                # pty = 0
                txa = parrx[0] + axpad + self.arrow_size / 2.0
                # txy = pay+txtpad
            else:
                pass  # raise NameError('arrowlegend not supported.')

            # FZ: Sudpip?
            ptx = self.arrow_size
            pty = 0
            # txy = pay + txtpad

            self.ax_ptm.arrow(pax,
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
            # self.ax_ptm.text(txa,
            #              txy,
            #              '|T|=1',
            #              horizontalalignment='center',
            #              verticalalignment='baseline',
            #              fontdict={'size':self.font_size,'weight':'bold'})
        # END: if self.plot_tipper.find('yes') == 0 ---------------------------

        # make a grid with color lines
        self.ax_ptm.grid(True, alpha=.3, which='both', color=(0.5, 0.5, 0.5))
        self.ax_ptm.set_axisbelow(True)
        if self.minorticks_on:
            plt.minorticks_on()  # turn on minor ticks automatically

        # ==> make a colorbar with appropriate colors
        if self.cb_position is None:
            self.ax_cb, kw = mcb.make_axes(self.ax_ptm,
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
            self.ax_cb = self.fig.add_axes(self.cb_position)

        if self.ellipse_cmap == 'mt_seg_bl2wh2rd':
            # make a color list
            self.clist = [
                (cc, cc, 1) for cc in np.arange(0, 1 + 1.0 / (nseg), 1.0 / (nseg))
            ] + [(1, cc, cc) for cc in np.arange(1, -1.0 / (nseg), -1.0 / (nseg))]

            # make segmented colormap
            mt_seg_bl2wh2rd = colors.ListedColormap(self.clist)

            # make bounds so that the middle is white
            bounds = np.arange(self.ellipse_range[0] - self.ellipse_range[2], self.ellipse_range[1] + 2 * self.ellipse_range[2], self.ellipse_range[2])

            # normalize the colors
            norms = colors.BoundaryNorm(bounds, mt_seg_bl2wh2rd.N)

            # make the colorbar
            self.cb = mcb.ColorbarBase(self.ax_cb,
                                       cmap=mt_seg_bl2wh2rd,
                                       norm=norms,
                                       orientation=self.cb_orientation,
                                       ticks=bounds[1:-1])
        else:
            if self.ellipse_cmap in list(mtcl.cmapdict.keys()):
                cmap_input = mtcl.cmapdict[self.ellipse_cmap]
            else:
                cmap_input = mtcl.cm.get_cmap(self.ellipse_cmap)
            self.cb = mcb.ColorbarBase(self.ax_cb,
                                       cmap=cmap_input,  # mtcl.cmapdict[cmap],
                                       norm=colors.Normalize(vmin=self.ellipse_range[0],
                                                             vmax=self.ellipse_range[1]),
                                       orientation=self.cb_orientation)

        # label the color bar accordingly
        self.cb.set_label(
            mtpl.ckdict[ck], fontdict={"size": self.font_size, "weight": "bold"}
        )

        # place the label in the correct location
        if self.cb_orientation == "horizontal":
            self.cb.ax.xaxis.set_label_position("top")
            self.cb.ax.xaxis.set_label_coords(0.5, 1.3)

        elif self.cb_orientation == "vertical":
            self.cb.ax.yaxis.set_label_position("right")
            self.cb.ax.yaxis.set_label_coords(1.25, 0.5)
            self.cb.ax.yaxis.tick_left()
            self.cb.ax.tick_params(axis="y", direction="in")

        # --> add reference ellipse:  (legend of ellipse size=1)
        # FZ: remove the following section if no show of Phi
        show_phi = False  # JingMingDuan does not want to show the black circle - it's not useful
        if show_phi is True:
            ref_ellip = patches.Ellipse((0, .0),
                                        width=self.ellipse_size,
                                        height=self.ellipse_size,
                                        angle=0)
            ref_ellip.set_facecolor((0, 0, 0))
            ref_ax_loc = list(self.ax_cb.get_position().bounds)
            ref_ax_loc[0] *= .95
            ref_ax_loc[1] -= .17
            ref_ax_loc[2] = .1
            ref_ax_loc[3] = .1
            self.ref_ax = self.fig.add_axes(ref_ax_loc, aspect='equal')
            self.ref_ax.add_artist(ref_ellip)
            self.ref_ax.set_xlim(-self.ellipse_size / 2. * 1.05, self.ellipse_size / 2. * 1.05)
            self.ref_ax.set_ylim(-self.ellipse_size / 2. * 1.05, self.ellipse_size / 2. * 1.05)
            plt.setp(self.ref_ax.xaxis.get_ticklabels(), visible=False)
            plt.setp(self.ref_ax.yaxis.get_ticklabels(), visible=False)
            self.ref_ax.set_title(r"$\Phi$ = 1")
            

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

    def save_figure(
        self,
        save_fn,
        file_format="pdf",
        orientation="portrait",
        fig_dpi=None,
        close_plot="y",
    ):
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

        sf = "_{0:.6g}".format(self.plot_freq)

        if fig_dpi is None:
            fig_dpi = self.fig_dpi

        # FZ: fixed the following logic
        if os.path.isdir(save_fn):  # FZ: assume save-fn is a directory
            if not os.path.exists(save_fn):
                os.mkdir(save_fn)

            # make a file name
            fname = "PTmap_DPI%s_%s_%sHz.%s" % (
                str(self.fig_dpi),
                self.ellipse_colorby,
                sf,
                file_format,
            )
            path2savefile = os.path.join(save_fn, fname)
        #            self.fig.savefig(path2savefile, dpi=fig_dpi, format=file_format, orientation=orientation,
        #                             bbox_inches='tight')
        else:  # FZ: assume save-fn is a path2file= "path2/afile.fmt"
            file_format = save_fn.split(".")[-1]
            if file_format is None or file_format not in ["png", "jpg"]:
                _logger.error(
                    "Error: output file name is not correctly provided:", save_fn
                )
                raise Exception("output file name is not correctly provided!!!")

            path2savefile = save_fn
            self.fig.savefig(
                path2savefile,
                dpi=fig_dpi,
                format=file_format,
                orientation=orientation,
                bbox_inches="tight",
            )
            # plt.clf()
            # plt.close(self.fig)

        if close_plot == "y":
            plt.clf()
            plt.close(self.fig)
        else:
            pass

        self.fig_fn = save_fn
        _logger.debug("Saved figure to: %s", self.fig_fn)

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

        self.fig.clf()
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
        stationmap = np.zeros((xlist.shape[0], ylist.shape[0]), dtype="|S8")

        station_location = {}  # a dict to store all mt_obj stations (lan lat)

        # put the information into the zeroed arrays
        for ii in range(nx):
            mt1 = self.mt_list[ii]

            # try to find the freq in the freq list of each file
            freqfind = [
                ff
                for ff, f2 in enumerate(mt1.Z.freq)
                if self.plot_freq * (1 - self.ftol)
                < f2
                < self.plot_freq * (1 + self.ftol)
            ]
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
                    stationmap[xyloc[ii, 0], xyloc[ii, 1]] = mt1.station[
                        self.station_id[0] : self.station_id[1]
                    ]
                except AttributeError:
                    stationmap[xyloc[ii, 0], xyloc[ii, 1]] = mt1.station

                station_location[stationmap[xyloc[ii, 0], xyloc[ii, 1]]] = (
                    mt1.longitude,
                    mt1.latitude,
                    mt1.freq[j2],
                )
            except IndexError:
                _logger.warn(
                    "Did not find {0:.5g} Hz for station {1}".format(
                        self.plot_freq, mt1.station
                    )
                )

        # ----------------------write files-------------------------------------
        svfn = "Map_{0:.6g}Hz".format(self.plot_freq)
        ptminfid = open(os.path.join(svpath, svfn + ".phimin"), "w")
        ptmaxfid = open(os.path.join(svpath, svfn + ".phimax"), "w")
        ptazmfid = open(os.path.join(svpath, svfn + ".azimuth"), "w")
        ptskwfid = open(os.path.join(svpath, svfn + ".skew"), "w")
        ptellfid = open(os.path.join(svpath, svfn + ".ellipticity"), "w")
        tprmgfid = open(os.path.join(svpath, svfn + ".tipper_mag_real"), "w")
        tprazfid = open(os.path.join(svpath, svfn + ".tipper_ang_real"), "w")
        tpimgfid = open(os.path.join(svpath, svfn + ".tipper_mag_imag"), "w")
        tpiazfid = open(os.path.join(svpath, svfn + ".tipper_ang_imag"), "w")
        statnfid = open(os.path.join(svpath, svfn + ".station"), "w")
        tablefid = open(os.path.join(svpath, svfn + ".table"), "w")

        for ly in range(ylist.shape[0]):
            for lx in range(xlist.shape[0]):
                # if there is nothing there write some spaces
                if phiminmap[lx, ly] == 0.0:
                    ptminfid.write("{0:^8}".format(" "))
                    ptmaxfid.write("{0:^8}".format(" "))
                    ptazmfid.write("{0:^8}".format(" "))
                    ptskwfid.write("{0:^8}".format(" "))
                    ptellfid.write("{0:^8}".format(" "))
                    tprmgfid.write("{0:^8}".format(" "))
                    tprazfid.write("{0:^8}".format(" "))
                    tpimgfid.write("{0:^8}".format(" "))
                    tpiazfid.write("{0:^8}".format(" "))
                    statnfid.write("{0:^8}".format(" "))

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
                    statnfid.write("{0:^8}".format(stationmap[lx, ly]))

            # make sure there is an end of line
            ptminfid.write("\n")
            ptmaxfid.write("\n")
            ptazmfid.write("\n")
            ptskwfid.write("\n")
            ptellfid.write("\n")
            tprmgfid.write("\n")
            tprazfid.write("\n")
            tpimgfid.write("\n")
            tpiazfid.write("\n")
            statnfid.write("\n")

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
        for ss in [
            "station",
            "lon",
            "lat",
            "phi_min",
            "phi_max",
            "skew",
            "ellipticity",
            "azimuth",
            "tip_mag_re",
            "tip_ang_re",
            "tip_mag_im",
            "tip_ang_im",
            "frequency",
        ]:
            tablefid.write("{0:^12}".format(ss))
        tablefid.write("\n")

        for ii in range(nx):
            xx, yy = xyloc[ii, 0], xyloc[ii, 1]
            if stationmap[xx, yy] is None or len(stationmap[xx, yy]) < 1:
                pass  # station not having the freq
            else:  # only those stations with the given freq
                # Station
                tablefid.write("{0:^12}".format(stationmap[xx, yy]))
                # lon
                tablefid.write(
                    mtpl.make_value_str(
                        station_location[stationmap[xx, yy]][0], spacing="{0:^12}"
                    )
                )
                # lat
                tablefid.write(
                    mtpl.make_value_str(
                        station_location[stationmap[xx, yy]][1], spacing="{0:^12}"
                    )
                )
                # phi_min
                tablefid.write(
                    mtpl.make_value_str(phiminmap[xx, yy], spacing="{0:^12}")
                )
                # phi_max
                tablefid.write(
                    mtpl.make_value_str(phimaxmap[xx, yy], spacing="{0:^12}")
                )
                # beta_skew
                tablefid.write(mtpl.make_value_str(betamap[xx, yy], spacing="{0:^12}"))
                # ellip
                tablefid.write(mtpl.make_value_str(ellipmap[xx, yy], spacing="{0:^12}"))
                # azimuth
                tablefid.write(
                    mtpl.make_value_str(azimuthmap[xx, yy], spacing="{0:^12}")
                )
                # tiprmag
                tablefid.write(mtpl.make_value_str(trmap[xx, yy], spacing="{0:^12}"))
                # tiprang
                tablefid.write(mtpl.make_value_str(trazmap[xx, yy], spacing="{0:^12}"))
                # tipimag
                tablefid.write(mtpl.make_value_str(timap[xx, yy], spacing="{0:^12}"))
                # tipiang
                tablefid.write(mtpl.make_value_str(tiazmap[xx, yy], spacing="{0:^12}"))
                # frequency
                tablefid.write(
                    mtpl.make_value_str(
                        station_location[stationmap[xx, yy]][2], spacing="{0:^12}"
                    )
                )
                tablefid.write("\n")

        tablefid.write("\n")

        _logger.info("Wrote files to {}".format(svpath))

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
    examples = os.path.join(base, "examples")
    data = os.path.join(examples, "data")
    edidir = os.path.join(data, "edi_files_2")

    edi_file_list = glob.glob(edidir + "/*.edi")
    savedir = "/tmp"

    plot_freq = 1e-2
    ptm_obj = PlotPhaseTensorMaps(
        fn_list=edi_file_list,
        plot_freq=plot_freq,
        ftol=0.2,
        interpolate=True,
        xpad=0.02,
        plot_tipper="yr",
        edgecolor="k",
        lw=0.1,
        alpha=1,
        minorticks_on=False,
        ellipse_size=0.2,
        ellipse_range=[-10, 10, 2],
        ellipse_colorby="skew",
        arrow_size=0.5,
        ellipse_cmap="mt_seg_bl2wh2rd",
        plot_yn="n",
    )

    # generate raster data
    lons, lats = np.meshgrid(np.linspace(136, 140, 50), np.linspace(-19, -22, 50))
    vals = np.exp(
        -np.sin(np.radians(lons) * 50) ** 2 - np.cos(np.radians(lats) * 70) ** 3
    )

    # plot with raster data
    f = plt.figure(figsize=(8, 6))
    ptm_obj.plot(
        fig=f,
        show=True,
        raster_dict={
            "lons": lons.flatten(),
            "lats": lats.flatten(),
            "vals": vals.flatten(),
            "levels": 50,
            "cmap": "rainbow",
            "cbar_title": "Arbitrary Units",
        },
    )

    ptm_obj.export_params_to_file(save_path=savedir)
