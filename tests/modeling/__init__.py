import shutil
from difflib import unified_diff

import matplotlib
import os

if os.name == "posix" and 'DISPLAY' not in os.environ:
    print("MATPLOTLIB: No Display found, using non-interactive svg backend")
    matplotlib.use('svg')
    import matplotlib.pyplot as plt
else:
    # matplotlib.use('svg')
    import matplotlib.pyplot as plt
    plt.ion()


def _diff_files(after, before, ignores=None):
    """
    compare two files using diff
    :param before:
    :param after:
    :return: the number count of different lines
    """

    with open(before) as f2p:
        before_lines = f2p.readlines()
    with open(after) as f1p:
        after_lines = f1p.readlines()

    if ignores:
        for ignored_term in ignores:
            before_lines = filter(lambda line: ignored_term not in line, before_lines)
            after_lines = filter(lambda line: ignored_term not in line, before_lines)

    msg = "Comparing {} and {}:\n".format(before, after)
    is_identical = False

    lines = [line for line in unified_diff(
        before_lines,
        after_lines,
        fromfile="baseline ({})".format(before),
        tofile="test ({})".format(after),
        n=0)]


    if lines:
        msg += "  Found differences:\n    " + "    ".join(lines)
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
