# -*- coding: utf-8 -*-
"""
mtplot_tools contains helper functions and classes for plotting



@author: jpeacock
"""
# ==============================================================================
from pathlib import Path
import numpy as np

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.colorbar as mcb
from matplotlib.lines import Line2D
import matplotlib.tri as tri

import mtpy.imaging.mtcolors as mtcl
from mtpy.utils.mtpy_logger import get_mtpy_logger

# =============================================================================
# Global Parameters
# =============================================================================

period_label_dict = dict(
    [(ii, "$10^{" + str(ii) + "}$") for ii in range(-20, 21)]
)


def get_period_limits(period):
    return (
        10 ** (np.floor(np.log10(period.min()))),
        10 ** (np.ceil(np.log10(period.max()))),
    )


# ==============================================================================
# Arrows properties for induction vectors
# ==============================================================================
class MTArrows:
    """
    Helper class to read a dictionary of arrow properties

    Arguments:
    -----------
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

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.arrow_size = 2.5
        self.arrow_head_length = 0.15 * self.arrow_size
        self.arrow_head_width = 0.1 * self.arrow_size
        self.arrow_lw = 0.5 * self.arrow_size
        self.arrow_threshold = 2
        self.arrow_color_imag = "c"
        self.arrow_color_real = "k"
        self.arrow_direction = 0

        # Set class property values from kwargs and pop them
        for v in vars(self):
            if v in list(kwargs.keys()):
                setattr(self, v, kwargs.pop(v, None))


# ==============================================================================
#  ellipse properties
# ==============================================================================
class MTEllipse:
    """
    helper class for getting ellipse properties from an input dictionary

    Arguments:
    -------------

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
              - 'normalized_skew' -> colors by
                              normalized_skew
                              see Booker, 2014
              - 'normalized_skew_seg' -> colors by
                             normalized_skew
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
                 'mt_rd2gr2bl']

             - 'mt_yl2rd'       --> yellow to red
             - 'mt_bl2yl2rd'    --> blue to yellow to red
             - 'mt_wh2bl'       --> white to blue
             - 'mt_rd2bl'       --> red to blue
             - 'mt_bl2wh2rd'    --> blue to white to red
             - 'mt_bl2gr2rd'    --> blue to green to red
             - 'mt_rd2gr2bl'    --> red to green to blue
             - 'mt_seg_bl2wh2rd' --> discrete blue to
                                     white to red

    """

    def __init__(self, **kwargs):
        self.ellipse_size = 2
        self.ellipse_colorby = "phimin"
        self.ellipse_range = (0, 90, 10)
        self.ellipse_cmap = "mt_bl2gr2rd"
        self.ellipse_spacing = 1

        # Set class property values from kwargs and pop them
        for v in vars(self):
            if v in list(kwargs.keys()):
                setattr(self, v, kwargs.pop(v, None))
        self.get_range()
        self.get_color_map()

    def get_color_map(self):
        """
        get a color map
        """
        if self.ellipse_colorby in ["skew_seg", "normalized_skew_seg"]:
            self.ellipse_cmap = "mt_seg_bl2wh2rd"

    def get_range(self):
        """
        get an appropriate range for the colorby
        """
        # set color ranges
        if self.ellipse_range[0] == self.ellipse_range[1]:
            if self.ellipse_colorby in [
                "skew",
                "skew_seg",
                "normalized_skew",
                "normalized_skew_seg",
            ]:

                self.ellipse_range = (-9, 9, 3)
            elif self.ellipse_colorby == "ellipticity":
                self.ellipse_range = (0, 1, 0.1)
            else:
                self.ellipse_range = (0, 90, 5)

    def get_pt_color_array(self, pt_object):
        """
        Get the appropriat color by array
        """

        # get the properties to color the ellipses by
        if (
            self.ellipse_colorby == "phiminang"
            or self.ellipse_colorby == "phimin"
        ):
            color_array = pt_object.phimin
        elif (
            self.ellipse_colorby == "phimaxang"
            or self.ellipse_colorby == "phimax"
        ):
            color_array = pt_object.phimax
        elif self.ellipse_colorby == "phidet":
            color_array = np.sqrt(abs(pt_object.det)) * (180 / np.pi)
        elif (
            self.ellipse_colorby == "skew"
            or self.ellipse_colorby == "skew_seg"
        ):
            color_array = pt_object.beta
        elif self.ellipse_colorby == "ellipticity":
            color_array = pt_object.ellipticity
        elif self.ellipse_colorby in ["strike", "azimuth"]:
            color_array = pt_object.azimuth % 180
            color_array[np.where(color_array > 90)] -= 180
        else:
            raise NameError(self.ellipse_colorby + " is not supported")
        return color_array

    @property
    def ellipse_properties(self):
        return {
            "size": self.ellipse_size,
            "range": self.ellipse_range,
            "cmap": self.ellipse_cmap,
            "colorby": self.ellipse_colorby,
            "spacing": self.ellipse_spacing,
        }


# ==============================================================================
# Plot settings
# ==============================================================================
class PlotSettings(MTArrows, MTEllipse):
    """
    Hold all the plot settings that one might need
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # figure properties:
        self.fig_num = 1
        self.fig_dpi = 150
        self.fig_size = None
        self.show_plot = True

        self.font_size = 7
        self.marker_size = 2.5
        self.marker_lw = 0.75
        self.marker_color = "b"
        self.marker = "v"
        self.lw = 1
        self.plot_title = None

        # line styles:
        self.xy_ls = ":"
        self.yx_ls = ":"
        self.det_ls = ":"
        self.skew_ls = ":"
        self.strike_ls = ":"

        # marker styles:
        self.xy_marker = "s"
        self.yx_marker = "o"
        self.det_marker = "v"
        self.skew_marker = "d"
        self.strike_inv_marker = "v"
        self.strike_pt_marker = "^"
        self.strike_tip_marker = ">"

        # marker color styles:
        self.xy_color = (0.25, 0.35, 0.75)
        self.yx_color = (0.75, 0.25, 0.25)
        self.det_color = (0.25, 0.75, 0.25)
        self.skew_color = (0.85, 0.35, 0)
        self.strike_inv_color = (0.2, 0.2, 0.7)
        self.strike_pt_color = (0.7, 0.2, 0.2)
        self.strike_tip_color = (0.2, 0.7, 0.2)

        # marker face color styles:
        self.xy_mfc = (0.25, 0.35, 0.75)
        self.yx_mfc = (0.75, 0.25, 0.25)
        self.det_mfc = (0.25, 0.75, 0.25)
        self.skew_mfc = (0.85, 0.35, 0)
        self.strike_inv_mfc = (0.2, 0.2, 0.7)
        self.strike_pt_mfc = (0.7, 0.2, 0.2)
        self.strike_tip_mfc = (0.2, 0.7, 0.2)

        # plot limits
        self.x_limits = None
        self.y_limits = None
        self.res_limits = None
        self.phase_limits = None
        self.tipper_limits = None
        self.strike_limits = None
        self.skew_limits = None
        self.pt_limits = None

        # Show Plot
        self.show_plot = True
        self.plot_tipper = "n"
        self.plot_pt = False
        self.plot_strike = False
        self.plot_skew = False

        self.text_size = 7
        self.text_weight = "normal"
        self.text_color = "k"
        self.text_ha = "center"
        self.text_va = "baseline"
        self.text_angle = 0
        self.text_x_pad = 0
        self.text_y_pad = 0

        self.subplot_left = 0.09
        self.subplot_right = 0.9
        self.subplot_bottom = 0.09
        self.subplot_top = 0.98
        self.subplot_wspace = None
        self.subplot_hspace = None

        # Set class property values from kwargs and pop them
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.cb_label_dict = {
            "phiminang": r"$\Phi_{min}$ (deg)",
            "phimin": r"$\Phi_{min}$ (deg)",
            "phimaxang": r"$\Phi_{max}$ (deg)",
            "phimax": r"$\Phi_{max}$ (deg)",
            "phidet": r"Det{$\Phi$} (deg)",
            "skew": r"Skew (deg)",
            "normalized_skew": r"Normalized Skew (deg)",
            "ellipticity": r"Ellipticity",
            "skew_seg": r"Skew (deg)",
            "normalized_skew_seg": r"Normalized Skew (deg)",
            "geometric_mean": r"$\sqrt{\Phi_{min} \cdot \Phi_{max}}$",
            "strike": r"Azimuth (deg)",
            "azimuth": r"Azimuth (deg)",
        }

        self.period_label_dict = period_label_dict

    def set_period_limits(self, period):
        """
        set period limits

        :return: DESCRIPTION
        :rtype: TYPE

        """

        return (
            10 ** (np.floor(np.log10(period.min()))),
            10 ** (np.ceil(np.log10(period.max()))),
        )

    def set_resistivity_limits(self, resistivity, mode="od", scale="log"):
        """
        set resistivity limits

        :param resistivity: DESCRIPTION
        :type resistivity: TYPE
        :param mode: DESCRIPTION, defaults to "od"
        :type mode: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if mode in ["od", "det", "det_only"]:
            nz_xy = np.nonzero(resistivity[:, 0, 1])
            nz_yx = np.nonzero(resistivity[:, 1, 0])
            limits = [
                10
                ** (
                    np.floor(
                        np.log10(
                            min(
                                [
                                    np.nanmin(resistivity[nz_xy, 0, 1]),
                                    np.nanmin(resistivity[nz_yx, 1, 0]),
                                ]
                            )
                        )
                    )
                ),
                10
                ** (
                    np.ceil(
                        np.log10(
                            max(
                                [
                                    np.nanmax(resistivity[nz_xy, 0, 1]),
                                    np.nanmax(resistivity[nz_yx, 1, 0]),
                                ]
                            )
                        )
                    )
                ),
            ]
        if mode == "d":
            nz_xx = np.nonzero(resistivity[:, 0, 1])
            nz_yy = np.nonzero(resistivity[:, 1, 0])
            limits = [
                10
                ** (
                    np.floor(
                        np.log10(
                            min(
                                [
                                    np.nanmin(resistivity[nz_xx, 0, 0]),
                                    np.nanmin(resistivity[nz_yy, 1, 1]),
                                ]
                            )
                        )
                    )
                ),
                10
                ** (
                    np.ceil(
                        np.log10(
                            max(
                                [
                                    np.nanmax(resistivity[nz_xx, 0, 0]),
                                    np.nanmax(resistivity[nz_yy, 1, 1]),
                                ]
                            )
                        )
                    )
                ),
            ]
        if scale == "log":
            if limits[0] == 0:
                limits[0] = 0.1
        return limits

    def set_phase_limits(self, phase, mode="od"):
        if mode in ["od", "det", "det_only"]:
            nz_xy = np.nonzero(phase[:, 0, 1])
            nz_yx = np.nonzero(phase[:, 1, 0])

            ph_min = min(
                [
                    np.nanmin(phase[nz_xy, 0, 1]),
                    np.nanmin(phase[nz_yx, 1, 0] + 180),
                ]
            )
            if ph_min > 0:
                ph_min = 0
            else:
                ph_min = round(ph_min / 5) * 5
            ph_max = max(
                [np.nanmax(phase[:, 0, 1]), np.nanmax(phase[:, 1, 0] + 180)]
            )
            if ph_max < 91:
                ph_max = 89.9
            else:
                ph_max = round(ph_max / 5) * 5
            return (ph_min, ph_max)
        elif mode == "d":
            ph_min = -180
            ph_max = 180

            return (ph_min, ph_max)

    @property
    def xy_error_bar_properties(self):
        """
        xy error bar properties for xy mode
        :return: DESCRIPTION
        :rtype: TYPE

        """
        return {
            "marker": self.xy_marker,
            "ms": self.marker_size,
            "mew": self.lw,
            "mec": self.xy_color,
            "color": self.xy_color,
            "ecolor": self.xy_color,
            "ls": self.xy_ls,
            "lw": self.lw,
            "capsize": self.marker_size,
            "capthick": self.lw,
        }

    @property
    def yx_error_bar_properties(self):
        """
        xy error bar properties for xy mode
        :return: DESCRIPTION
        :rtype: TYPE

        """
        return {
            "marker": self.yx_marker,
            "ms": self.marker_size,
            "mew": self.lw,
            "mec": self.yx_color,
            "color": self.yx_color,
            "ecolor": self.yx_color,
            "ls": self.yx_ls,
            "lw": self.lw,
            "capsize": self.marker_size,
            "capthick": self.lw,
        }

    @property
    def det_error_bar_properties(self):
        """
        xy error bar properties for xy mode
        :return: DESCRIPTION
        :rtype: TYPE

        """
        return {
            "marker": self.det_marker,
            "ms": self.marker_size,
            "mew": self.lw,
            "mec": self.det_color,
            "color": self.det_color,
            "ecolor": self.det_color,
            "ls": self.det_ls,
            "lw": self.lw,
            "capsize": self.marker_size,
            "capthick": self.lw,
        }

    @property
    def font_dict(self):
        return {"size": self.font_size + 2, "weight": "bold"}

    @property
    def arrow_real_properties(self):
        return {
            "lw": self.arrow_lw,
            "facecolor": self.arrow_color_real,
            "edgecolor": self.arrow_color_real,
            "head_width": self.arrow_head_width,
            "head_length": self.arrow_head_length,
            "length_includes_head": False,
        }

    @property
    def arrow_imag_properties(self):
        return {
            "lw": self.arrow_lw,
            "facecolor": self.arrow_color_imag,
            "edgecolor": self.arrow_color_imag,
            "head_width": self.arrow_head_width,
            "head_length": self.arrow_head_length,
            "length_includes_head": False,
        }

    @property
    def text_dict(self):
        return {
            "size": self.text_size,
            "weight": self.text_weight,
            "rotation": self.text_angle,
            "color": self.text_color,
        }


# =============================================================================
# Base
# =============================================================================
class PlotBase(PlotSettings):
    """
    base class for plotting objects

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.logger = get_mtpy_logger(
            f"{self.__class__.__module__}.{self.__class__.__name__}"
        )

        self._basename = self.__class__.__name__.lower()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return f"Plotting {self._basename}"

    def __repr__(self):
        return self.__str__()

    def _set_subplot_params(self):
        # set some parameters of the figure and subplot spacing
        plt.rcParams["font.size"] = self.font_size
        plt.rcParams["figure.subplot.bottom"] = self.subplot_bottom
        plt.rcParams["figure.subplot.top"] = self.subplot_top
        plt.rcParams["figure.subplot.left"] = self.subplot_left
        plt.rcParams["figure.subplot.right"] = self.subplot_right

        if self.subplot_wspace is not None:
            plt.rcParams["figure.subplot.wspace"] = self.subplot_wspace
        if self.subplot_hspace is not None:
            plt.rcParams["figure.subplot.hspace"] = self.subplot_hspace

    def plot(self):
        pass

    def save_plot(
        self,
        save_fn,
        file_format="pdf",
        orientation="portrait",
        fig_dpi=None,
        close_plot=True,
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

            **file_format** : [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps...

            **orientation** : [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait

            **fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the fig_dpi will be that at
                          which the figure was made.  I don't think that
                          it can be larger than fig_dpi of the figure.

            **close_plot** : [ true | false ]
                             * True will close the plot after saving.
                             * False will leave plot open

        :Example: ::

            >>> # to save plot as jpg
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotResPhase(r'/home/MT/mt01.edi')
            >>> p1.save_plot(r'/home/MT/figures', file_format='jpg')

        """

        if fig_dpi is None:
            fig_dpi = self.fig_dpi
        save_fn = Path(save_fn)
        if not save_fn.is_dir():
            file_format = save_fn.suffix[1:]
        else:
            save_fn = save_fn.joinpath(f"{self._basename}.{file_format}")
        self.fig.savefig(
            save_fn, dpi=fig_dpi, format=file_format, orientation=orientation
        )

        if close_plot:
            plt.close(self.fig)
        else:
            pass
        self.fig_fn = save_fn
        self.logger.info(f"Saved figure to: {self.fig_fn}")

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


# =============================================================================
#  plotting functions
# =============================================================================
def plot_resistivity(ax, period, resistivity, error, **properties):
    """
    plot apparent resistivity to the given axis with given properties

    :param ax: DESCRIPTION
    :type ax: TYPE
    :param resistivity: DESCRIPTION
    :type resistivity: TYPE
    :param period: DESCRIPTION
    :type period: TYPE
    :param **properties: DESCRIPTION
    :type **properties: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    nz = np.nonzero(resistivity)

    return plot_errorbar(
        ax,
        period[nz],
        resistivity[nz],
        y_error=error[nz],
        **properties,
    )


def plot_phase(ax, period, phase, error, yx=False, **properties):
    """
    plot apparent resistivity to the given axis with given properties

    :param ax: DESCRIPTION
    :type ax: TYPE
    :param resistivity: DESCRIPTION
    :type resistivity: TYPE
    :param period: DESCRIPTION
    :type period: TYPE
    :param **properties: DESCRIPTION
    :type **properties: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    # need this for the yx component
    nz = np.nonzero(phase)
    if yx:
        return plot_errorbar(
            ax,
            period[nz],
            phase[nz] + 180,
            y_error=error[nz],
            **properties,
        )
    return plot_errorbar(
        ax,
        period[nz],
        phase[nz],
        y_error=error[nz],
        **properties,
    )


# ==============================================================================
# grid data onto a map view
# ==============================================================================
def grid_data(data_array, x, y, nx=None, ny=None):
    """
    Project data onto a regular grid for plotting.


    Arguments:
    -----------
        **data_array**: np.ndarray (len(x), len(y))
                        array of data values to be gridded

        **x**: np.ndarray(len(x))
               array of values that coorespond

        **nx**: int
                number of cells in the x-direction.  If none, 2 times the
                number of x components

        **ny**: int
                number of cells in the x-direction.  If none, 2 times the
                number of y components

    Returns:
    ---------
        **grid_array**: np.ndarray(nx, ny)
                        array of data set on a regular grid

        **xg**: np.ndarray(nx, ny)
                array of x-grid values

        **yg**: np.ndarray(nx, ny)
                array of y-grid values


    """

    if nx is None:
        nx = 2 * len(x)
    if ny is None:
        ny = 2 * len(y)
    # create evenly spaced intervals to grid over
    xi = np.linspace(x.min(), x.max(), num=nx, endpoint=True)
    yi = np.linspace(y.min(), y.max(), num=ny, endpoint=True)

    xg, yg = np.meshgrid(xi, yi)

    grid_array = mlab.griddata(x, y, data_array, xg, yg)

    return grid_array, xg, yg


# ==============================================================================
# function for writing values to file
# ==============================================================================
def make_value_str(
    value,
    value_list=None,
    spacing="{0:^8}",
    value_format="{0: .2f}",
    append=False,
    add=False,
):
    """
    helper function for writing values to a file, takes in a value and either
    appends or adds value to value_list according to the spacing and format of
    the string.

    Arguments:
    ----------
        **value** : float

        **value_list** : list of values converted to strings

        **spacing** : spacing of the string that the value will be converted
                      to.

        **value_format** : format of the string that the value is being
                            coverted to.

        **append** : [ True | False]
                     if True then appends the value to value list

        **add** : [ True | False ]
                  if True adds value string to the other value strings in
                  value_list

    Returns:
    --------
        **value_list** : the input value_list with the new value either
                        added or appended.
        or

        **value_str** : value string if add and append are false
    """

    value_str = spacing.format(value_format.format(value))

    if append is True:
        value_list.append(value_str)
        return value_list
    if add is True:
        value_list += value_str
        return value_list
    if append == False and add == False:
        return value_str
    return value_list


# ==============================================================================
# function for error bar plots
# ==============================================================================
def plot_errorbar(ax, x_array, y_array, y_error=None, x_error=None, **kwargs):
    """
    convinience function to make an error bar instance

    Arguments:
    ------------
        **ax** : matplotlib.axes instance
                 axes to put error bar plot on

        **x_array** : np.ndarray(nx)
                      array of x values to plot

        **y_array** : np.ndarray(nx)
                      array of y values to plot

        **y_error** : np.ndarray(nx)
                      array of errors in y-direction to plot

        **x_error** : np.ndarray(ns)
                      array of error in x-direction to plot

        **color** : string or (r, g, b)
                    color of marker, line and error bar

        **marker** : string
                     marker type to plot data as

        **mew** : string
                     marker edgewidth

        **ms** : float
                 size of marker

        **ls** : string
                 line style between markers

        **lw** : float
                 width of line between markers

        **e_capsize** : float
                        size of error bar cap

        **e_capthick** : float
                         thickness of error bar cap

        **picker** : float
                     radius in points to be able to pick a point.


    Returns:
    ---------
        **errorbar_object** : matplotlib.Axes.errorbar
                              error bar object containing line data,
                              errorbars, etc.
    """
    # this is to make sure error bars plot in full and not just a dashed line
    if x_error is not None:
        x_err = x_error
    else:
        x_err = None
    if y_error is not None:

        y_err = y_error
    else:
        y_err = None
    plt_settings = {
        "color": "k",
        "marker": "x",
        "mew": 1,
        "mec": "k",
        "ms": 2,
        "ls": ":",
        "lw": 1,
        "capsize": 2,
        "capthick": 0.5,
        "ecolor": "k",
        "elinewidth": 1,
        "picker": None,
    }

    for key, value in kwargs.items():
        plt_settings[key] = value
    errorbar_object = ax.errorbar(
        x_array, y_array, xerr=x_err, yerr=y_err, **plt_settings
    )
    return errorbar_object


def add_colorbar_axis(ax, fig):
    # add colorbar for PT
    axpos = ax.get_position()
    cb_position = (
        axpos.bounds[0] - 0.0575,
        axpos.bounds[1] + 0.02,
        0.01,
        axpos.bounds[3] * 0.75,
    )

    cbax = fig.add_axes(cb_position)
    return cbax


def plot_pt_lateral(
    ax,
    pt_obj,
    color_array,
    ellipse_properties,
    y_shift=0,
    fig=None,
    edge_color=None,
    n_index=0,
):
    """

    :param ax: DESCRIPTION
    :type ax: TYPE
    :param pt_obj: DESCRIPTION
    :type pt_obj: TYPE
    :param color_array: DESCRIPTION
    :type color_array: TYPE
    :param ellipse_properties: DESCRIPTION
    :type ellipse_properties: TYPE
    :param bounds: DESCRIPTION, defaults to None
    :type bounds: TYPE, optional
    :return: DESCRIPTION
    :rtype: TYPE

    """
    bounds = None
    try:
        ellipse_properties["range"][2]
    except IndexError:
        ellipse_properties["range"][2] = 3
    if ellipse_properties["cmap"] == "mt_seg_bl2wh2rd":
        bounds = np.arange(
            ellipse_properties["range"][0],
            ellipse_properties["range"][1] + ellipse_properties["range"][2],
            ellipse_properties["range"][2],
        )
        nseg = float(
            (ellipse_properties["range"][1] - ellipse_properties["range"][0])
            / (2 * ellipse_properties["range"][2])
        )
    # -------------plot ellipses-----------------------------------
    for ii, ff in enumerate(1.0 / pt_obj.freq):
        # make sure the ellipses will be visable
        if pt_obj.phimax[ii] == 0:
            continue
        eheight = (
            pt_obj.phimin[ii] / pt_obj.phimax[ii] * ellipse_properties["size"]
        )
        ewidth = (
            pt_obj.phimax[ii] / pt_obj.phimax[ii] * ellipse_properties["size"]
        )

        # create an ellipse scaled by phimin and phimax and oriented
        # along the azimuth which is calculated as clockwise but needs
        # to be plotted counter-clockwise hence the negative sign.
        ellipd = patches.Ellipse(
            (np.log10(ff) * ellipse_properties["spacing"], y_shift),
            width=ewidth,
            height=eheight,
            angle=90 - pt_obj.azimuth[ii],
        )

        ax.add_patch(ellipd)

        # get ellipse color
        ellipd.set_facecolor(
            mtcl.get_plot_color(
                color_array[ii],
                ellipse_properties["colorby"],
                ellipse_properties["cmap"],
                ellipse_properties["range"][0],
                ellipse_properties["range"][1],
                bounds=bounds,
            )
        )
        if edge_color is not None:
            ellipd.set_edgecolor(edge_color)
    # set axis properties
    ax.set_ylim(
        ymin=-1.5 * ellipse_properties["size"],
        ymax=y_shift + 1.5 * ellipse_properties["size"],
    )
    cbax = None
    cbpt = None
    if n_index == 0:
        if fig is not None:
            cbax = add_colorbar_axis(ax, fig)
        if ellipse_properties["cmap"] == "mt_seg_bl2wh2rd":
            # make the colorbar
            nseg = float(
                (
                    ellipse_properties["range"][1]
                    - ellipse_properties["range"][0]
                )
                / (2 * ellipse_properties["range"][2])
            )
            cbpt = make_color_list(
                cbax,
                nseg,
                ellipse_properties["range"][0],
                ellipse_properties["range"][1],
                ellipse_properties["range"][2],
            )
        else:
            cbpt = mcb.ColorbarBase(
                cbax,
                cmap=mtcl.cmapdict[ellipse_properties["cmap"]],
                norm=colors.Normalize(
                    vmin=ellipse_properties["range"][0],
                    vmax=ellipse_properties["range"][1],
                ),
                orientation="vertical",
            )
        cbpt.set_ticks(
            [
                ellipse_properties["range"][0],
                (
                    ellipse_properties["range"][1]
                    - ellipse_properties["range"][0]
                )
                / 2,
                ellipse_properties["range"][1],
            ]
        )
        cbpt.set_ticklabels(
            [
                f"{ellipse_properties['range'][0]:.0f}",
                f"{(ellipse_properties['range'][1] - ellipse_properties['range'][0]) / 2:.0f}",
                f"{ellipse_properties['range'][1]:.0f}",
            ]
        )

        cbpt.ax.yaxis.set_label_position("left")
        cbpt.ax.yaxis.set_label_coords(-1.05, 0.5)
        cbpt.ax.yaxis.tick_right()
        cbpt.ax.tick_params(axis="y", direction="in")
    return cbax, cbpt


def plot_tipper_lateral(
    axt,
    t_obj,
    plot_tipper,
    real_properties,
    imag_properties,
    font_size=6,
    legend=True,
    zero_reference=False,
):
    """

    :param axt: DESCRIPTION
    :type axt: TYPE
    :param t_obj: DESCRIPTION
    :type t_obj: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    if plot_tipper.find("y") == 0 or plot_tipper:
        txr = t_obj.mag_real * np.cos(np.deg2rad(-t_obj.angle_real))
        tyr = t_obj.mag_real * np.sin(np.deg2rad(-t_obj.angle_real))

        txi = t_obj.mag_imag * np.cos(np.deg2rad(-t_obj.angle_imag))
        tyi = t_obj.mag_imag * np.sin(np.deg2rad(-t_obj.angle_imag))

        nt = len(txr)
        period = 1.0 / t_obj.freq
        x_limits = get_period_limits(period)

        tiplist = []
        tiplabel = []

        if plot_tipper.find("r") > 0:
            line = Line2D([0], [0], color=real_properties["facecolor"], lw=1)
            tiplist.append(line)
            tiplabel.append("real")
        if plot_tipper.find("i") > 0:
            line = Line2D([0], [0], color=imag_properties["facecolor"], lw=1)
            tiplist.append(line)
            tiplabel.append("imag")
        for aa in range(nt):
            xlenr = txr[aa] * np.log10(period[aa])
            xleni = txi[aa] * np.log10(period[aa])

            if xlenr == 0 and xleni == 0:
                continue
            # --> plot real arrows
            if plot_tipper.find("r") > 0:
                axt.arrow(
                    np.log10(period[aa]),
                    0,
                    xlenr,
                    tyr[aa],
                    **real_properties,
                )
            # --> plot imaginary arrows
            if plot_tipper.find("i") > 0:
                axt.arrow(
                    np.log10(period[aa]),
                    0,
                    xleni,
                    tyi[aa],
                    **imag_properties,
                )
        # make a line at 0 for reference
        if zero_reference:
            axt.plot(np.log10(period), [0] * nt, "k", lw=0.5)
        if legend:
            axt.legend(
                tiplist,
                tiplabel,
                loc="upper left",
                markerscale=1,
                borderaxespad=0.01,
                labelspacing=0.07,
                handletextpad=0.2,
                borderpad=0.1,
                prop={"size": 6},
            )
        # set axis properties

        axt.set_xlim(np.log10(x_limits[0]), np.log10(x_limits[1]))

        tklabels = []
        xticks = []

        for tk in axt.get_xticks():
            try:
                tklabels.append(period_label_dict[tk])
                xticks.append(tk)
            except KeyError:
                pass
        axt.set_xticks(xticks)
        axt.set_xticklabels(tklabels, fontdict={"size": font_size})
        # need to reset the x_limits caouse they get reset when calling
        # set_ticks for some reason
        axt.set_xlim(np.log10(x_limits[0]), np.log10(x_limits[1]))

        # axt.set_xscale('log', nonpositive='clip')
        tmax = max([np.nanmax(tyr), np.nanmax(tyi)])
        if tmax > 1:
            tmax = 0.899
        tmin = min([np.nanmin(tyr), np.nanmin(tyi)])
        if tmin < -1:
            tmin = -0.899
        tipper_limits = (tmin - 0.1, tmax + 0.1)
        axt.set_ylim(tipper_limits)
        axt.grid(
            True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25
        )
    return axt, tiplist, tiplabel


def get_log_tick_labels(ax, spacing=1):
    """

    :param ax: DESCRIPTION
    :type ax: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    tklabels = []
    xticks = []
    for tk in ax.get_xticks():
        try:
            tklabels.append(period_label_dict[tk / spacing])
            xticks.append(tk)
        except KeyError:
            pass
    return tklabels, xticks


def make_color_list(cbax, nseg, ckmin, ckmax, ckstep):
    """ """

    # make a color list
    clist = [
        (cc, cc, 1) for cc in np.arange(0, 1 + 1.0 / (nseg), 1.0 / (nseg))
    ] + [(1, cc, cc) for cc in np.arange(1, -1.0 / (nseg), -1.0 / (nseg))]

    # make segmented colormap
    mt_seg_bl2wh2rd = colors.ListedColormap(clist)

    # make bounds so that the middle is white
    bounds = np.arange(ckmin - ckstep, ckmax + 2 * ckstep, ckstep)

    # normalize the colors
    norms = colors.BoundaryNorm(bounds, mt_seg_bl2wh2rd.N)

    # make the colorbar
    return mcb.ColorbarBase(
        cbax,
        cmap=mt_seg_bl2wh2rd,
        norm=norms,
        orientation="vertical",
        ticks=bounds[1:-1],
    )


def round_to_step(num, base=5):
    return base * round(num / base)


def add_raster(
    ax,
    raster_dict={
        "lons": [],
        "lats": [],
        "vals": [],
        "levels": 50,
        "cmap": "rainbow",
        "cbar_title": "Arbitrary units",
    },
    refpoint=(0, 0),
    add_colorbar=True,
    cb_orientation="horizontal",
):
    lons = np.array(raster_dict["lons"])
    lats = np.array(raster_dict["lats"])

    # retain masking if a masked array is passed in
    if type(raster_dict["vals"]) == np.ma.core.MaskedArray:
        vals = np.ma.masked_array(raster_dict["vals"])
    else:
        vals = np.array(raster_dict["vals"])
    lons -= refpoint[0]
    lats -= refpoint[1]
    levels = raster_dict.pop("levels", 50)
    cmap = raster_dict.pop("cmap", "rainbow")
    cbar_title = raster_dict.pop("cbar_title", "Arbitrary Units")

    # if a 2D array provided, can use contourf and no need to triangulate
    triangulate = True
    if len(vals.shape) > 1:
        if lons.shape == lats.shape == vals.shape:
            triangulate = False
        elif (lons.shape == vals.shape[1]) and (lats.shape == vals.shape[0]):
            triangulate = False
    if triangulate:
        assert (
            len(lons) == len(lats) == len(vals)
        ), "Lons, Lats and Vals must all have the same length"
        triangulation = tri.Triangulation(lons, lats)
        ax.tricontourf(
            triangulation,
            vals,
            levels=np.linspace(vals.min(), vals.max(), levels),
            cmap=cmap,
        )
    else:
        ax.contourf(
            lons,
            lats,
            vals,
            levels=np.linspace(vals.min(), vals.max(), levels),
            cmap=cmap,
        )
    if add_colorbar:
        cbax, kw = mcb.make_axes(ax, orientation=cb_orientation, shrink=0.35)

        if cb_orientation == "horizontal":
            cbax.ax.set_xlabel(cbar_title)
            cbax.ax.xaxis.set_label_position("top")
            cbax.ax.xaxis.set_label_coords(0.5, 1.3)
        else:
            cbax.ax.set_ylabel(cbar_title, fontweight="bold")
            cbax.ax.yaxis.set_label_position("right")
            cbax.ax.yaxis.set_label_coords(1.25, 0.5)
            cbax.ax.yaxis.tick_left()
            cbax.ax.tick_params(axis="y", direction="in")
