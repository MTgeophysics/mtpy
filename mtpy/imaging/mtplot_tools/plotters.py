# -*- coding: utf-8 -*-
"""
Simple plotters elements that can be assembled in various plotting classes

Created on Sun Sep 25 15:27:28 2022

:author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np

import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.colorbar as mcb
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt

from mtpy.imaging.mtcolors import get_plot_color
from .utils import (
    period_label_dict,
    get_period_limits,
    add_colorbar_axis,
    make_color_list,
)

# =============================================================================


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
    for ii, ff in enumerate(1.0 / pt_obj.frequency):
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
            get_plot_color(
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
                cmap=plt.get_cmap(ellipse_properties["cmap"]),
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
    if t_obj is None:
        return None, None, None

    if axt is None:
        return None, None, None

    if plot_tipper.find("y") == 0 or plot_tipper:
        txr = t_obj.mag_real * np.cos(np.deg2rad(-t_obj.angle_real))
        tyr = t_obj.mag_real * np.sin(np.deg2rad(-t_obj.angle_real))

        txi = t_obj.mag_imag * np.cos(np.deg2rad(-t_obj.angle_imag))
        tyi = t_obj.mag_imag * np.sin(np.deg2rad(-t_obj.angle_imag))

        nt = len(txr)
        period = 1.0 / t_obj.frequency
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
