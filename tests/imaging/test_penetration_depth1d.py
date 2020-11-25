"""
Run these unit test at commandline.
$   cd u25656@PC /e/Githubz/mtpy2 (develop)
$   nosetests tests/imaging/
"""
import os.path

import matplotlib.pyplot as plt

from mtpy.imaging.penetration import ZComponentError
from mtpy.imaging.penetration_depth1d import plot_edi_dir
from mtpy.imaging.penetration_depth1d import plot_edi_file
from tests.imaging import ImageTestCase, ImageCompare


class TestPenetration_depth1d(ImageTestCase):
    def test_plot_edi_dir(self):
        """
        testing plotting all edi files in a given dir
        :return:
        """
        # all plots
        plot_edi_dir("data/edifiles")
        # plt.close()

    @ImageCompare(fig_size=(8, 6))
    def test_plot_edi_file(self):
        """
        testing plotting a single edi file
        :return:
        """
        plot_edi_file("data/edifiles/15125A.edi")

    @ImageCompare(fig_size=(8, 6))
    def test_plot_edi_file_zxy(self):
        """
        testing ploting zxy of a single edi file
        :return:
        """
        plot_edi_file("data/edifiles/15125A.edi", ["zxy"])

    @ImageCompare(fig_size=(8, 6))
    def test_plot_edi_file_zyx(self):
        """
        testing plotting zyx of a single edi file
        :return:
        """
        plot_edi_file("data/edifiles/15125A.edi", ["zyx"])

    @ImageCompare(fig_size=(8, 6))
    def test_plot_edi_file_det(self):
        """
        testing plotting det of a single edi file
        :return:
        """
        plot_edi_file("data/edifiles/15125A.edi", ["det"])

    @ImageCompare(fig_size=(8, 6))
    def test_plot_edi_file_zxy_zyx(self):
        """
        testing plotting zxy & zyx of a single edi file
        :return:
        """
        plot_edi_file("data/edifiles/15125A.edi", ["zxy", "zyx"])

    @ImageCompare(fig_size=(8, 6))
    def test_plot_edi_file_zxy_det(self):
        """
        testing plotting zxy & det of a single edi file
        :return:
        """
        plot_edi_file("data/edifiles/15125A.edi", ["zxy", "det"])

    @ImageCompare(fig_size=(8, 6))
    def test_plot_edi_file_zyx_det(self):
        """
        testing plotting zyx & det of a single edi file
        :return:
        """
        plot_edi_file("data/edifiles/15125A.edi", ["zyx", "det"])

    @ImageCompare(fig_size=(8, 6))
    def test_plot_edi_file_unknown_type(self):
        """
        testing plotting zyx and an unknown of a single edi file
        :return:
        """
        plot_edi_file("data/edifiles/15125A.edi", ["zyx", "dat"])

    def test_plot_edi_file_empty_rholist(self):
        """
        testing plotting an empty rholist of a single edi file
        :return:
        """
        with self.assertRaises(ZComponentError):
            plot_edi_file("data/edifiles/15125A.edi", [])

    def test_plot_edi_file_save_image(self):
        """
        testing saving plot of a single edi file
        :return:
        """
        fname = os.path.join(self._temp_dir, "TestPenetration_depth1d.png")
        if os.path.isfile(fname):
            os.remove(fname)  # remove test file if already exist
        plot_edi_file("data/edifiles/15125A.edi", savefile=fname)
        assert os.path.isfile(fname)
