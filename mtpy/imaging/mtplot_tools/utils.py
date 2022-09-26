# -*- coding: utf-8 -*-
"""
Utility functions for plotting

Created on Sun Sep 25 15:49:01 2022

:author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np

import matplotlib.colors as colors
import matplotlib.colorbar as mcb

# =============================================================================

period_label_dict = dict(
    [(ii, "$10^{" + str(ii) + "}$") for ii in range(-20, 21)]
)


def get_period_limits(period):
    return (
        10 ** (np.floor(np.log10(period.min()))),
        10 ** (np.ceil(np.log10(period.max()))),
    )


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
