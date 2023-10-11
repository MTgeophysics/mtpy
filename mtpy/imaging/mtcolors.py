# -*- coding: utf-8 -*-
"""
Created on Tue May 14 18:05:59 2013

@author: jpeacock-pr
"""
# =============================================================================
# Imports
# =============================================================================
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib import cm

import numpy as np

# ==============================================================================
# Make some color maps for plotting
# ==============================================================================
# yellow to red
ptcmapdict = {
    "red": ((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
    "green": ((0.0, 0.0, 1.0), (1.0, 0.0, 1.0)),
    "blue": ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
}

mt_yl2rd = colors.LinearSegmentedColormap("mt_yl2rd", ptcmapdict, 256)

# blue to yellow to red
skcmapdict = {
    "red": (
        (0.0, 0.0, 0.0),
        (0.5, 1.0, 1.0),
        (0.5, 0.0, 1.0),
        (1.0, 1.0, 1.0),
    ),
    "green": (
        (0.0, 1.0, 0.0),
        (0.5, 1.0, 0.0),
        (0.5, 0.0, 1.0),
        (1.0, 0.0, 1.0),
    ),
    "blue": (
        (0.0, 0.0, 1.0),
        (0.5, 0.0, 1.0),
        (0.5, 0.1, 0.1),
        (1.0, 0.1, 0.1),
    ),
}

mt_bl2yl2rd = colors.LinearSegmentedColormap("mt_bl2yl2rd", skcmapdict, 256)

# blue to white to red
skcmapdict2 = {
    "red": (
        (0.0, 0.0, 0.0),
        (0.25, 0.0, 0.0),
        (0.5, 0.8, 1.0),
        (0.75, 1.0, 1.0),
        (1.0, 0.4, 1.0),
    ),
    "green": (
        (0.0, 0.0, 0.0),
        (0.25, 0.0, 0.0),
        (0.5, 0.9, 0.9),
        (0.75, 0.0, 0.0),
        (1.0, 0.0, 0.0),
    ),
    "blue": (
        (0.0, 0.0, 0.4),
        (0.25, 1.0, 1.0),
        (0.5, 1.0, 0.8),
        (0.75, 0.0, 0.0),
        (1.0, 0.0, 0.0),
    ),
}

mt_bl2wh2rd = colors.LinearSegmentedColormap("mt_bl2wh2rd", skcmapdict2, 256)


# blue to white to red in segmented colors
mt_seg_bl2wh2rd = colors.ListedColormap(
    (
        (0, 0, 1),
        (0.5, 0.5, 1),
        (0.75, 0.75, 1),
        (0.9, 0.9, 1),
        (1, 1, 1),
        (1.0, 0.9, 0.9),
        (1, 0.75, 0.75),
        (1, 0.5, 0.5),
        (1, 0, 0),
    )
)

# white to blue
ptcmapdict3 = {
    "red": ((0.0, 1.0, 1.0), (1.0, 0.0, 0.0)),
    "green": ((0.0, 1.0, 1.0), (1.0, 0.0, 0.0)),
    "blue": ((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
}
mt_wh2bl = colors.LinearSegmentedColormap("mt_wh2bl", ptcmapdict3, 256)

# white to orange
cmapdict_wh2or = {
    "red": ((0.0, 1.0, 1.0), (1.0, 0.95, 0.0)),
    "green": ((0.0, 1.0, 1.0), (1.0, 0.45, 0.95)),
    "blue": ((0.0, 1.0, 1.0), (1.0, 0, 0)),
}
mt_wh2or = colors.LinearSegmentedColormap("mt_wh2or", cmapdict_wh2or, 256)

# red to blue
rtcmapdict = {
    "red": ((0.0, 0.0, 1.0), (1.0, 0.0, 1.0)),
    "green": ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
    "blue": ((0.0, 1.0, 0.0), (1.0, 1.0, 0.0)),
}
mt_rd2bl = colors.LinearSegmentedColormap("mt_rd2bl", rtcmapdict, 256)

# blue to green to red
ptcmapdict4 = {
    "red": (
        (0.0, 0.0, 0.0),
        (0.25, 0.0, 0.0),
        (0.5, 0.9, 1.0),
        (0.75, 1.0, 1.0),
        (1.0, 0.45, 1.0),
    ),
    "green": (
        (0.0, 0.0, 0.0),
        (0.25, 0.5, 0.5),
        (0.5, 1.0, 1.0),
        (0.75, 0.5, 0.5),
        (1.0, 0.0, 0.0),
    ),
    "blue": (
        (0.0, 0.0, 0.45),
        (0.25, 1.0, 1.0),
        (0.5, 1.0, 0.9),
        (0.75, 0.0, 0.0),
        (1.0, 0.0, 0.0),
    ),
}
mt_bl2gr2rd = colors.LinearSegmentedColormap("mt_bl2gr2rd", ptcmapdict4, 256)

# red to green to blue
ptcmapdict4 = {
    "red": (
        (0.0, 0.0, 0.45),
        (0.25, 1.0, 1.0),
        (0.5, 1.0, 0.9),
        (0.75, 0.0, 0.0),
        (1.0, 0.0, 0.0),
    ),
    "green": (
        (0.0, 0.0, 0.0),
        (0.25, 0.5, 0.5),
        (0.5, 1.0, 1.0),
        (0.75, 0.5, 0.5),
        (1.0, 0.0, 0.0),
    ),
    "blue": (
        (0.0, 0.0, 0.0),
        (0.25, 0.0, 0.0),
        (0.5, 0.9, 1.0),
        (0.75, 1.0, 1.0),
        (1.0, 0.45, 1.0),
    ),
}
mt_rd2gr2bl = colors.LinearSegmentedColormap("mt_rd2gr2bl", ptcmapdict4, 256)

mtcmapdict2 = {
    "red": (
        (0.0, 0.4, 0.2),
        (0.20, 0.992, 0.992),
        (0.50, 0.953, 0.953),
        (0.62, 0.384, 0.384),
        (1.0, 0.12, 0.12),
    ),
    "green": (
        (0.0, 0.0392, 0.0392),
        (0.20, 0.3098, 0.3098),
        (0.50, 0.953, 0.953),
        (0.62, 0.529, 0.529),
        (1.0, 0.094, 0.094),
    ),
    "blue": (
        (0.0, 0.00, 0.00),
        (0.20, 0.0, 0.0),
        (0.50, 0.953, 0.953),
        (0.62, 1.0, 1.0),
        (1.0, 0.45, 0.45),
    ),
}
mt_rd2wh2bl = colors.LinearSegmentedColormap("mt_rd2wh2bl", mtcmapdict2, 256)

mtcmapdict3 = {
    "red": (
        (0.0, 0.2, 0.2),
        (0.25, 0.2, 0.2),
        (0.5, 1.0, 1.0),
        (0.75, 1.0, 1.0),
        (1.0, 0.2, 0.5),
    ),
    "green": (
        (0.0, 0.3, 0.3),
        (0.25, 0.3, 0.3),
        (0.5, 1.0, 1.0),
        (0.75, 0.3, 0.3),
        (1.0, 0.0, 0.0),
    ),
    "blue": (
        (0.0, 0.75, 0.45),
        (0.25, 0.85, 0.85),
        (0.5, 1.0, 1.0),
        (0.75, 0.0, 0.0),
        (1.0, 0.0, 0.0),
    ),
}

mt_rd2wh2bl_r = colors.LinearSegmentedColormap(
    "mt_rd2wh2bl_r", mtcmapdict3, 256
)

rdylbu_data = {
    "blue": [
        [0.0, 0.000, 0.000],
        [0.1, 0.000, 0.000],
        [0.2, 0.000, 0.000],
        [0.3, 0.000, 0.000],
        [0.4, 0.000, 0.000],
        [0.5, 1.000, 1.000],
        [0.6, 0.990, 0.990],
        [0.7, 0.950, 0.950],
        [0.8, 0.900, 0.900],
        [0.9, 0.550, 0.550],
        [1.0, 0.250, 0.250],
    ],
    "green": [
        [0.0, 0.000, 0.000],
        [0.1, 0.100, 0.100],
        [0.2, 0.400, 0.400],
        [0.3, 0.800, 0.800],
        [0.4, 0.900, 0.900],
        [0.5, 1.000, 1.000],
        [0.6, 0.900, 0.900],
        [0.7, 0.800, 0.800],
        [0.8, 0.400, 0.400],
        [0.9, 0.100, 0.100],
        [1.0, 0.000, 0.000],
    ],
    "red": [
        [0.0, 0.250, 0.250],
        [0.1, 0.550, 0.550],
        [0.2, 0.900, 0.900],
        [0.3, 0.950, 0.950],
        [0.4, 0.990, 0.990],
        [0.5, 1.000, 1.000],
        [0.6, 0.000, 0.000],
        [0.7, 0.000, 0.000],
        [0.8, 0.000, 0.000],
        [0.9, 0.000, 0.000],
        [1.0, 0.000, 0.000],
    ],
}

mt_rdylbu = colors.LinearSegmentedColormap("mt_rdylbu", rdylbu_data, 256)

# Combine the lower and upper range of the terrain colormap with a gap in the middle
# to let the coastline appear more prominently.
# inspired by https://stackoverflow.com/questions/31051488/combining-two-matplotlib-colormaps
colors_undersea = cm.terrain(np.linspace(0, 0.17, 56))
colors_land = cm.terrain(np.linspace(0.25, 1, 200))


# combine them and build a new colormap
color_list = np.vstack((colors_undersea, colors_land))
cut_terrain_map = colors.LinearSegmentedColormap.from_list(
    "cut_terrain", color_list
)

MT_CMAP_DICT = {
    "mt_yl2rd": mt_yl2rd,
    "mt_bl2yl2rd": mt_bl2yl2rd,
    "mt_wh2bl": mt_wh2bl,
    "mt_rd2bl": mt_rd2bl,
    "mt_bl2wh2rd": mt_bl2wh2rd,
    "mt_seg_bl2wh2rd": mt_seg_bl2wh2rd,
    "mt_bl2gr2rd": mt_bl2gr2rd,
    "mt_rd2gr2bl": mt_rd2gr2bl,
    "mt_wh2or": mt_wh2or,
    "mt_rd2wh2bl": mt_rd2wh2bl,
    "mt_rd2wh2bl_r": mt_rd2wh2bl_r,
    "mt_rdylbu": mt_rdylbu,
    "cut_terrain": cut_terrain_map,
}


def register_cmaps(cmap_dict):
    for key, value in cmap_dict.items():
        plt.register_cmap(key, value)


def get_color(cvar, cmap):
    """
    gets the color to plot for the given color map

    """
    return cm.get_cmap(cmap)(cvar)


def get_plot_color(colorx, comp, cmap, ckmin=None, ckmax=None, bounds=None):
    """
    gets the color for the given compnent, color array and cmap

    Note: we now use the linearSegmentedColorMap objects, instead of the get_color function
    """

    # get face color info
    if comp in [
        "phimin",
        "phimax",
        "phidet",
        "ellipticity",
        "geometric_mean",
        "azimuth",
        "strike",
    ]:
        if ckmin is None or ckmax is None:
            raise IOError("Need to input min and max values for plotting")
        """
        cvar = (colorx-ckmin)/(ckmax-ckmin)
        if cmap == 'mt_bl2wh2rd' or cmap == 'mt_bl2yl2rd' or \
           cmap == 'mt_bl2gr2rd' or cmap == 'mt_rd2gr2bl' or \
           cmap == 'mt_rd2wh2bl' or cmap == 'mt_rd2wh2bl_r':
            cvar = 2*cvar-1

        return get_color(cvar, cmap)
        """
        norm = colors.Normalize(ckmin, ckmax)
        return cm.get_cmap(cmap)(norm(colorx))

    elif comp == "skew" or comp == "normalized_skew":
        """
        cvar = 2*colorx/(ckmax-ckmin)
        return get_color(cvar, cmap)
        """

        norm = colors.Normalize(ckmin, ckmax)
        return cm.get_cmap(cmap)(norm(colorx))

    elif comp == "skew_seg" or comp == "normalized_skew_seg":
        if bounds is None:
            raise IOError("Need to input bounds for segmented colormap")
        """
        for bb in range(bounds.shape[0]):
            if colorx >= bounds[bb] and colorx < bounds[bb+1]:
                cvar = float(bounds[bb])/bounds.max()
                return get_color(cvar, cmap)

            #if the skew is extremely negative make it blue
            elif colorx < bounds[0]:
                cvar = -1.0
                return get_color(cvar, cmap)

            #if skew is extremely positive make it red
            elif colorx > bounds[-1]:
                cvar = 1.0
                return get_color(cvar, cmap)
        """
        norm = colors.Normalize(bounds[0], bounds[-1])
        step = abs(bounds[1] - bounds[0])
        ### need to get the color into a bin so as to not smear the colors.

        if colorx > max(bounds):
            colorx = max(bounds)
        elif colorx < min(bounds):
            colorx = min(bounds)
        elif abs(colorx) <= step:
            colorx = 0
        else:
            colorx = int(
                step
                * round(
                    float(colorx - np.sign(colorx) * (abs(colorx) % step))
                    / step
                )
            )
        return cm.get_cmap(cmap)(norm(colorx))
    else:
        raise NameError("color key " + comp + " not supported")


def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    colors_i = np.concatenate((np.linspace(0, 1.0, N), (0.0, 0.0, 0.0, 0.0)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1.0, N + 1)
    cdict = {}
    for ki, key in enumerate(("red", "green", "blue")):
        cdict[key] = [
            (indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki])
            for i in range(N + 1)
        ]
    # Return colormap object.
    return colors.LinearSegmentedColormap(cmap.name + "_%d" % N, cdict, 1024)


class FixPointNormalize(colors.Normalize):
    """
    Inspired by https://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib
    Subclassing Normalize to obtain a colormap with a fixpoint
    somewhere in the middle of the colormap.
    This may be useful for a `terrain` map, to set the "sea level"
    to a color in the blue/turquise range.
    """

    def __init__(
        self, vmin=None, vmax=None, sealevel=0, col_val=0.21875, clip=False
    ):
        # sealevel is the fix point of the colormap (in data units)
        self.sealevel = sealevel
        # col_val is the color value in the range [0,1] that should represent the sealevel.
        self.col_val = col_val
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.sealevel, self.vmax], [0, self.col_val, 1]
        return np.ma.masked_array(np.interp(value, x, y))
