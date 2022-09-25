# -*- coding: utf-8 -*-
"""
Base classes for plotting classes

:author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================
from pathlib import Path
import numpy as np

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.colorbar as mcb
from matplotlib.lines import Line2D

import mtpy.imaging.mtcolors as mtcl

from .plot_settings import PlotSettings
from .map_interpolation_tools import interpolate_to_map
from mtpy.utils.mtpy_logger import get_mtpy_logger

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

        plt.close(self.fig)
        self.plot()


class PlotBaseMaps(PlotBase):
    """
    Base object for plot classes that use map views, includes methods for
    interpolation.

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def interpolate_to_map(
        self,
        plot_array,
        component,
        cell_size=0.002,
        n_padding_cells=10,
        interpolation_method="delaunay",
    ):
        """
        interpolate points onto a 2d map.

        :param plot_array: DESCRIPTION
        :type plot_array: TYPE
        :param component: DESCRIPTION
        :type component: TYPE
        :param cell_size: DESCRIPTION, defaults to 0.002
        :type cell_size: TYPE, optional
        :param n_padding_cells: DESCRIPTION, defaults to 10
        :type n_padding_cells: TYPE, optional
        :param interpolation_method: DESCRIPTION, defaults to "delaunay"
        :type interpolation_method: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        return interpolate_to_map(
            plot_array,
            component,
            cell_size=cell_size,
            n_padding_cells=n_padding_cells,
            interpolation_method=interpolation_method,
        )

    def _get_interpolated_z(self, tf):
        return np.array(
            [
                [
                    tf.z_interp_dict["zxx"]["real"](1 / self.plot_period)[0]
                    + 1j
                    * tf.z_interp_dict["zxx"]["imag"](1.0 / self.plot_period)[
                        0
                    ],
                    tf.z_interp_dict["zxy"]["real"](1.0 / self.plot_period)[0]
                    + 1j
                    * tf.z_interp_dict["zxy"]["imag"](1.0 / self.plot_period)[
                        0
                    ],
                ],
                [
                    tf.z_interp_dict["zyx"]["real"](1.0 / self.plot_period)[0]
                    + 1j
                    * tf.z_interp_dict["zyx"]["imag"](1.0 / self.plot_period)[
                        0
                    ],
                    tf.z_interp_dict["zyy"]["real"](1.0 / self.plot_period)[0]
                    + 1j
                    * tf.z_interp_dict["zyy"]["imag"](1.0 / self.plot_period)[
                        0
                    ],
                ],
            ]
        )

    def _get_interpolated_z_err(self, tf):
        return np.array(
            [
                [
                    tf.z_interp_dict["zxx"]["err"](self.plot_freq)[0],
                    tf.z_interp_dict["zxy"]["err"](self.plot_freq)[0],
                ],
                [
                    tf.z_interp_dict["zyx"]["err"](self.plot_freq)[0],
                    tf.z_interp_dict["zyy"]["err"](self.plot_freq)[0],
                ],
            ]
        )

    def _get_interpolated_t(self, tf):
        return np.array(
            [
                [
                    [
                        tf.t_interp_dict["tzx"]["real"](self.plot_freq)[0]
                        + 1j
                        * tf.t_interp_dict["tzx"]["imag"](self.plot_freq)[0],
                    ],
                    [
                        tf.z_interp_dict["tzy"]["real"](self.plot_freq)[0]
                        + 1j
                        * tf.t_interp_dict["tzy"]["imag"](self.plot_freq)[0],
                    ],
                ]
            ]
        )

    def _get_interpolated_t_err(self, tf):
        return np.array(
            [
                [
                    [
                        tf.t_interp_dict["tzx"]["err"](self.plot_freq)[0],
                    ],
                    [
                        tf.t_interp_dict["tzy"]["err"](self.plot_freq)[0],
                    ],
                ]
            ]
        )
