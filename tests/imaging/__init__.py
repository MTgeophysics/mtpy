from __future__ import print_function

import os
import sys
import shutil
from unittest import TestCase
import matplotlib

from mtpy.utils.decorator import ImageCompare
from mtpy.utils.mtpylog import MtPyLog

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


ImageCompare.print_image_testing_note(file=sys.stderr)


class ImageTestCase(TestCase):
    @classmethod
    def setUpClass(cls):
        reset_matplotlib()
        cls._temp_dir = "tests/temp/{}".format(cls.__name__.split('.')[-1])
        if os.path.isdir(cls._temp_dir):
            shutil.rmtree(cls._temp_dir)
        os.mkdir(cls._temp_dir)

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def setUp(self):
        if plt.get_fignums():
            plt.clf()
        reset_matplotlib()
        if plt.isinteractive():
            plt.show(block=False)  # show an empty window first for drawing

    def tearDown(self):
        if plt.isinteractive() and plt.get_fignums():
            plt.pause(1)
        if plt.get_fignums():
            plt.close("all")
