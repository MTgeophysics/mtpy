from unittest import TestCase
import os

from mtpy.imaging.penetration_depth3d import plot4all_periods
from mtpy.imaging.penetration_depth3d import plot_latlon_depth_profile

class TestPenetration_depth3d(TestCase):
    def setUp(self):
        self._temp_dir = "tests/temp"
        if not os.path.isdir(self._temp_dir):
            os.mkdir(self._temp_dir)
        self._edifiles_small = "tests/data/edifiles"

    # def test_plot4all_periods(self):
    #     plot4all_periods(self._edifiles_small)

    def test_plot_latlon_depth_profile_period_index(self):
        plot_latlon_depth_profile(self._edifiles_small, 10, 'det', savefig=False)

    def test_plot_latlon_depth_profile_no_period(self):
        try:
            plot_latlon_depth_profile(self._edifiles_small, savefig=False)
        except Exception:
            pass

    def test_plot_latlon_depth_profile_period(self):
        plot_latlon_depth_profile(self._edifiles_small, 2.857, savefig=False)
