from unittest import TestCase

import os.path

# configure matplotlib for testing
import matplotlib

matplotlib.use('Agg')  # comment out this line if you want to see the plots
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
        testing plotting all edi files in a given dir
        :return:
        """
        # all plots
        plot_edi_dir("tests/data/edifiles")
        # plt.close()

    def test_plot_edi_file_all(self):
        """
        testing plotting a single edi file
        :return:
        """
        plot_edi_file("tests/data/edifiles/15125A.edi")

    def test_plot_edi_file_zxy(self):
        """
        testing ploting zxy of a single edi file
        :return:
        """
        plot_edi_file("tests/data/edifiles/15125A.edi", ['zxy'])

    def test_plot_edi_file_zyx(self):
        """
        testing plotting zyx of a single edi file
        :return:
        """
        plot_edi_file("tests/data/edifiles/15125A.edi", ['zyx'])

    def test_plot_edi_file_det(self):
        """
        testing plotting det of a single edi file
        :return:
        """
        plot_edi_file("tests/data/edifiles/15125A.edi", ['det'])

    def test_plot_edi_file_zxy_zyx(self):
        """
        testing plotting zxy & zyx of a single edi file
        :return:
        """
        plot_edi_file("tests/data/edifiles/15125A.edi", ['zxy', 'zyx'])

    def test_plot_edi_file_zxy_det(self):
        """
        testing plotting zxy & det of a single edi file
        :return:
        """
        plot_edi_file("tests/data/edifiles/15125A.edi", ['zxy', 'det'])

    def test_plot_edi_file_zyx_det(self):
        """
        testing plotting zyx & det of a single edi file
        :return:
        """
        plot_edi_file("tests/data/edifiles/15125A.edi", ['zyx', 'det'])

    def test_plot_edi_file_unknown_type(self):
        """
        testing plotting zyx and an unknown of a single edi file
        :return:
        """
        plot_edi_file("tests/data/edifiles/15125A.edi", ['zyx', 'dat'])

    def test_plot_edi_file_empty_rholist(self):
        """
        testing plotting an empty rholist of a single edi file
        :return:
        """
        plot_edi_file("tests/data/edifiles/15125A.edi", [])

    def test_plot_edi_file_save_inage(self):
        """
        testing saving plot of a single edi file
        :return:
        """
        fname = os.path.join(self._temp_dir,"TestPenetration_depth1d.jpg")
        if os.path.isfile(fname):
            os.remove(fname)    # remove test file if already exist
        plot_edi_file("tests/data/edifiles/15125A.edi", savefile=fname)
        assert (os.path.isfile(fname))
        pass
