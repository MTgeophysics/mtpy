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

from mtpy.utils.mtpy_logger import get_mtpy_logger

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
        self.plot_tipper = False
        self.plot_pt = False

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

        # Set class property values from kwargs and pop them
        for v in vars(self):
            if v in list(kwargs.keys()):
                setattr(self, v, kwargs.pop(v, None))
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

        self.period_label_dict = dict(
            [(ii, "$10^{" + str(ii) + "}$") for ii in range(-20, 21)]
        )

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

    def set_resistivity_limits(self, resistivity, mode="od"):
        """
        set resistivity limits
        
        :param resistivity: DESCRIPTION
        :type resistivity: TYPE
        :param mode: DESCRIPTION, defaults to "od"
        :type mode: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if mode == "od":
            return (
                10
                ** (
                    np.floor(
                        np.log10(
                            min(
                                [
                                    np.nanmin(resistivity[:, 0, 1]),
                                    np.nanmin(resistivity[:, 1, 0]),
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
                                    np.nanmax(resistivity[:, 0, 1]),
                                    np.nanmax(resistivity[:, 1, 0]),
                                ]
                            )
                        )
                    )
                ),
            )
        if mode == "d":
            return (
                10
                ** (
                    np.floor(
                        np.log10(
                            min(
                                [
                                    np.nanmin(resistivity[:, 0, 0]),
                                    np.nanmin(resistivity[:, 1, 1]),
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
                                    np.nanmax(resistivity[:, 0, 0]),
                                    np.nanmax(resistivity[:, 1, 1]),
                                ]
                            )
                        )
                    )
                ),
            )

    def set_phase_limits(self, phase, mode="od"):
        if mode == "od":
            ph_min = 0
            ph_max = 90
            if min(phase[:, 0, 1]) < 0 or min(phase[:, 1, 0] + 180) < 0:
                ph_min = min([min(phase[:, 0, 1]), min(phase[:, 1, 0])])
                if ph_min > 0:
                    ph_min = 0
            else:
                ph_min = 0
            if max(phase[:, 0, 1]) > 90 or max(phase[:, 1, 0] + 180) > 90:
                ph_max = min([max(phase[:, 0, 1]), max(phase[:, 1, 0] + 180)])
                if ph_max < 91:
                    ph_max = 89.9
            else:
                ph_max = 89.9
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
            file_format = save_fn.suffix
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
