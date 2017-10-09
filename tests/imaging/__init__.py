from __future__ import print_function

import os
import shutil
from unittest import TestCase

import matplotlib

import sys

from matplotlib import pyplot as plt

from mtpy.utils.decorator import ImageCompare

if os.name == "posix" and 'DISPLAY' not in os.environ:
    print("MATPLOTLIB: No Display found, using non-interactive svg backend", sys.stderr)
    matplotlib.use('svg')
    import matplotlib.pyplot as plt
else:
    # matplotlib.use('svg')
    import matplotlib.pyplot as plt

    plt.ion()


def reset_matplotlib():
    # save some important params
    interactive = matplotlib.rcParams['interactive']
    backend = matplotlib.rcParams['backend']
    # reset
    matplotlib.rcdefaults()  # reset the rcparams to default
    # recover
    matplotlib.rcParams['backend'] = backend
    matplotlib.rcParams['interactive'] = interactive


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
        plt.clf()

    def tearDown(self):
        if plt.isinteractive():
            plt.pause(1)
        plt.close("all")
