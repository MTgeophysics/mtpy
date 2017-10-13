from __future__ import print_function

import os
import sys

import matplotlib
import sip
from PyQt4.QtGui import QApplication

from mtpy.utils.decorator import ImageCompare
from mtpy.utils.mtpylog import MtPyLog

sip.setdestroyonexit(False)

app = QApplication(sys.argv)

if os.name == "posix" and 'DISPLAY' not in os.environ:
    print("MATPLOTLIB: No Display found, using non-interactive svg backend", sys.stderr)
    matplotlib.use('svg')
    import matplotlib.pyplot as plt
else:
    # matplotlib.use('svg')
    import matplotlib.pyplot as plt
    plt.ion()

MtPyLog().get_mtpy_logger(__name__).info("Testing using matplotlib backend {}".format(matplotlib.rcParams['backend']))


def reset_matplotlib():
    # save some important params
    interactive = matplotlib.rcParams['interactive']
    backend = matplotlib.rcParams['backend']
    # reset
    matplotlib.rcdefaults()  # reset the rcparams to default
    # recover
    matplotlib.rcParams['backend'] = backend
    matplotlib.rcParams['interactive'] = interactive
    logger = MtPyLog().get_mtpy_logger(__name__)
    logger.info("Testing using matplotlib backend {}".format(matplotlib.rcParams['backend']))

