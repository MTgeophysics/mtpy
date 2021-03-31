import six
from matplotlib import colors as mcolors
from .plot_parameter import PlotParameter

COLORS = list(six.iteritems(mcolors.cnames))
SIMPLE_COLORS = [
    "b",  # blue
    "g",  # green
    "r",  # red
    "c",  # cyan
    "m",  # magenta
    "y",  # yellow
    "k",  # black
    "w",  # white
]
# # add the single letter colors
# for name, rgb in six.iteritems(mcolors.ColorConverter.colors):
#     hex_ = mcolors.rgb2hex(rgb)
#     COLORS.append((name, hex_))
# sort by name
COLORS.sort(key=lambda c: c[0])
