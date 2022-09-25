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

import mtpy.imaging.mtcolors as mtcl

from . import PlotSettings
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


def add_raster(ax, raster_fn, add_colorbar=True, **kwargs):
    """
    Add a raster to an axis using rasterio

    :param raster_fn: DESCRIPTION
    :type raster_fn: TYPE
    :param **kwargs: DESCRIPTION
    :type **kwargs: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    import rasterio
    from rasterio.plot import show

    tif = rasterio.open(raster_fn)
    ax2 = show(tif, ax=ax, **kwargs)
    cb = None
    if add_colorbar:
        im = ax2.get_images()[0]
        fig = ax2.get_figure()
        cb = fig.colorbar(im, ax=ax)

    return ax2, cb
