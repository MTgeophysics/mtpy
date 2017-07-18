import os
from unittest import TestCase

# configure matplotlib for testing
import matplotlib.pyplot as plt

plt.ion()

from mtpy.imaging.penetration_depth3d import plot_latlon_depth_profile

class TestPenetration_depth3d(TestCase):
    @classmethod
    def setUpClass(cls):
        cls._temp_dir = "tests/temp"
        if not os.path.isdir(cls._temp_dir):
            os.mkdir(cls._temp_dir)
        cls._edifiles_small = "tests/data/edifiles"

    # def test_plot_many_periods(self):
    #     plot_many_periods(self._edifiles_small)

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def tearDown(self):
        plt.pause(1)

    def test_plot_latlon_depth_profile_period_index(self):
        plot_latlon_depth_profile(self._edifiles_small, 10, 'det', savefig=False)

    def test_plot_latlon_depth_profile_no_period(self):
        try:
            plot_latlon_depth_profile(self._edifiles_small, savefig=False)
            assert(False)  # if this statement reached, it is wrong
        except Exception as ex:
            print (ex)
            assert(True)

    def test_plot_latlon_depth_profile_period(self):
        plot_latlon_depth_profile(self._edifiles_small, 2.857, savefig=False)
