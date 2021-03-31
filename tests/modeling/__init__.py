from __future__ import print_function

import shutil
from difflib import unified_diff

import matplotlib
import os

import sys

from matplotlib import pyplot as plt

if os.name == "posix" and "DISPLAY" not in os.environ:
    print(
        "MATPLOTLIB: No Display found, using non-interactive svg backend",
        file=sys.stderr,
    )
    matplotlib.use("svg")
    import matplotlib.pyplot as plt
else:
    # matplotlib.use('svg')
    import matplotlib.pyplot as plt

    plt.ion()


def diff_files(after, before, ignores=None):
    """
    compare two files using diff
    :param ignores:
    :param before:
    :param after:
    :return: the number count of different lines
    """

    with open(before) as f2p:
        before_lines = f2p.readlines()
    with open(after) as f1p:
        after_lines = f1p.readlines()

    before_lines = [line.strip() for line in before_lines]
    after_lines = [line.strip() for line in after_lines]

    if ignores:
        for ignored_term in ignores:
            before_lines = [line for line in before_lines if ignored_term not in line]
            after_lines = [line for line in before_lines if ignored_term not in line]

    msg = "Comparing {} and {}:\n".format(before, after)

    lines = [
        line
        for line in unified_diff(
            before_lines,
            after_lines,
            fromfile="baseline ({})".format(before),
            tofile="test ({})".format(after),
            n=0,
        )
    ]

    if lines:
        msg += "  Found differences:\n\t" + "\n\t".join(lines)
        is_identical = False
    else:
        msg += " NO differences found."
        is_identical = True

    return is_identical, msg


def _clean_recreate(adir):
    if os.path.exists(adir):
        # clear dir if it already exist
        shutil.rmtree(adir)

    os.mkdir(adir)


# def show_patcher(show_func):
#     """
#     patch the plt.show() if interactive is enabled to display and then close the plot after 1 second
#     so plt.show() will not block the script and the figure is still visible to the user
#     :param show_func:
#     :return:
#     """
#
#     def new_show_func(*args, **kwargs):
#         stuff = show_func(*args, **kwargs)
#         # wait 1 second for the image to show on screen
#         figManager = plt.gcf()
#         if figManager is not None:
#             canvas = figManager.canvas
#             # if canvas.figure.stale:
#             #     canvas.draw()
#             # show(block=False)
#             try:
#                 canvas.start_event_loop(1)  # wait time = 1
#             except NotImplementedError:
#                 pass
#             finally:
#                 pass
#         plt.close()
#         return stuff
#
#     return new_show_func if plt.isinteractive() else show_func
