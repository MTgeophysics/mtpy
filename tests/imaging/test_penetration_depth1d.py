from unittest import TestCase

import os.path

# configure matplotlib for testing
import matplotlib

matplotlib.use('Agg')  # remove this line if you want to see the plots
import matplotlib.pyplot as plt
# plt.ioff()

from mtpy.imaging.penetration_depth1d import plot_edi_dir
from mtpy.imaging.penetration_depth1d import plot_edi_file


class TestPenetration_depth1d(TestCase):
    def setUp(self):
        self._temp_dir = "tests/temp"
        if not os.path.isdir(self._temp_dir):
            os.mkdir(self._temp_dir)


    def test_plot_edi_dir(self):
        """
        testing ploting all edi files in a given dir
        :return:
        """
        # all plots
        plot_edi_dir("tests/data/edifiles")
        # plt.close()

    def test_plot_edi_file(self):
        """
        testing ploting a single edi file
        :return:
        """
        plot_edi_file("tests/data/edifiles/15125A.edi")
        #['zxy', 'zyx', 'det']
        # zxy
        plot_edi_file("tests/data/edifiles/15125A.edi", ['zxy'])
        # zyx
        plot_edi_file("tests/data/edifiles/15125A.edi", ['zyx'])
        # det
        plot_edi_file("tests/data/edifiles/15125A.edi", ['det'])
        # zxy & zyx
        plot_edi_file("tests/data/edifiles/15125A.edi", ['zxy', 'zyx'])
        # zxy & det
        plot_edi_file("tests/data/edifiles/15125A.edi", ['zxy', 'det'])
        # zyx & det
        plot_edi_file("tests/data/edifiles/15125A.edi", ['zyx', 'det'])
        # plot unknown type
        plot_edi_file("tests/data/edifiles/15125A.edi", ['zyx', 'dat'])
        # plot empty set of tpyes
        plot_edi_file("tests/data/edifiles/15125A.edi", [])
        # save image
        fname = os.path.join(self._temp_dir,"TestPenetration_depth1d.jpg")
        if os.path.isfile(fname):
            os.remove(fname)    # remove test file if already exist
        plot_edi_file("tests/data/edifiles/15125A.edi", savefile=fname)
        assert (os.path.isfile(fname))
        pass
